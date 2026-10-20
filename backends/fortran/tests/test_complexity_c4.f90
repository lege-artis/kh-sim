! test_complexity_c4.f90 - C4 complexity-claim audit (LEGE-ARTIS language doctrine
! v0.1 s4.4 C4) for the hand-rolled FFT of the Fortran reference (src/kh_fft.f90).
!
! Claim under audit: the radix-2 DIT FFT is O(N log N) (Cooley-Tukey; the same
! algorithm as the Pascal FFT1D, after Press et al. (1986) S12.2).
!
! Procedure (doctrine C4): execution time t(N) of one transform, 5 runs per size,
! the median is used; least-squares slope of ln t over ln(N log2 N);
! PASS when 0.9 <= slope <= 1.1. An O(N^2) defect gives a slope near 2.
! The gate applies to the 1D transform, the routine that carries the claim
! (_audit/textbook-citations-allowlist.yaml: fft1d_inplace, O(N log N)).
! The 2D wrappers carry no complexity claim of their own; their slope is printed
! for information. It sits lower (about 0.88..0.94 measured), because the fixed
! cost per row/column call (copies, twiddle set-up) is large against N log2 N at
! 8 x 8 and decays over the range - a lower-order term, not a complexity defect.
!
!   1D: kh_fft_forward_2d / kh_fft_inverse_2d with nx = 1, ny = N,
!       N = 2^6 .. 2^16 (fft1d_inplace is private; with nx = 1 the second
!       pass is N transforms of length 1, O(N)).
!   2D: square nx = ny = 2^3 .. 2^8, so N = nx*ny = 2^6 .. 2^16 (even
!       exponents); cost nx*ny*(log2 nx + log2 ny) = N log2 N.
!
! One unit of work is a forward transform followed by the inverse one (which
! normalises), so the data stay bounded however often the pair is repeated. Each
! run repeats the pair until at least MIN_RUN_S seconds of process CPU time
! (cpu_time) have elapsed, so the clock resolution is irrelevant and preemption by
! other processes on a shared runner does not enter the fit.
!
! Build (release - C4 measures run time; doctrine: release-no-LTO):
!   gfortran -O2 -Jbuild/c4 -o build/c4/test_complexity_c4 \
!       src/kh_constants.f90 src/kh_fft.f90 tests/test_complexity_c4.f90
! Run:   ./build/c4/test_complexity_c4      exit 0 = PASS (1D gate)
! ASCII-only source.

program test_complexity_c4
  use kh_fft, only: kh_fft_forward_2d, kh_fft_inverse_2d
  implicit none

  integer,  parameter :: NRUNS = 5
  real(8),  parameter :: MIN_RUN_S = 0.10_8
  real(8),  parameter :: SLOPE_LO = 0.9_8, SLOPE_HI = 1.1_8

  real(8) :: xs(11), ys(11), s1, s2
  integer :: e, k, cnt
  logical :: ok1

  write(*,'(a)') 'C4 complexity audit - kh_fft.f90 radix-2 FFT, claim O(N log N)'
  write(*,'(a,i0,a,f4.2,a)') '  ', NRUNS, ' runs per size (median), >= ', &
                             MIN_RUN_S, ' s CPU per run'

  write(*,'(a)') '  1D (nx = 1, ny = N), N = 2^6..2^16'
  cnt = 0
  do e = 6, 16
    cnt = cnt + 1
    call measure(1, 2**e, xs(cnt), ys(cnt))
  end do
  s1 = slope(xs, ys, cnt)
  ok1 = (s1 >= SLOPE_LO) .and. (s1 <= SLOPE_HI)
  write(*,'(a,f6.3,a,f3.1,a,f3.1,a)') 'C4 FFT 1D slope = ', s1, '  (gate [', &
        SLOPE_LO, ', ', SLOPE_HI, '])'

  write(*,'(a)') '  2D (nx = ny), N = nx*ny = 2^6..2^16'
  cnt = 0
  do k = 3, 8
    cnt = cnt + 1
    call measure(2**k, 2**k, xs(cnt), ys(cnt))
  end do
  s2 = slope(xs, ys, cnt)
  write(*,'(a,f6.3,a)') 'C4 FFT 2D slope = ', s2, &
        '  (information only - no complexity claim on the 2D wrappers)'

  if (ok1) then
    write(*,'(a)') 'RESULT: PASS'
    stop
  else
    write(*,'(a)') 'RESULT: FAIL'
    error stop 1
  end if

contains

  subroutine measure(nx, ny, x, y)
    integer, intent(in)  :: nx, ny
    real(8), intent(out) :: x, y
    complex(8), allocatable :: a(:,:)
    real(8) :: samples(NRUNS), t, lg
    integer :: i, j, r, n

    n = nx * ny
    allocate(a(0:nx-1, 0:ny-1))
    do j = 0, ny-1
      do i = 0, nx-1
        a(i, j) = cmplx(sin(0.37_8*(i + nx*j)), cos(0.11_8*(i + nx*j)), kind=8)
      end do
    end do
    call kh_fft_forward_2d(a, nx, ny)      ! warm-up
    call kh_fft_inverse_2d(a, nx, ny)
    do r = 1, NRUNS
      samples(r) = time_pair(a, nx, ny)
    end do
    t = median(samples)
    lg = log(real(n, 8)) / log(2.0_8)
    x = log(real(n, 8) * lg)
    y = log(t)
    write(*,'(a,i4,a,i6,a,i6,a,f11.3,a,f8.3,a)') '    nx = ', nx, '  ny = ', ny, &
          '  N = ', n, '  t = ', t*1.0e6_8, ' us  t/(N log2 N) = ', &
          t*1.0e9_8/(real(n, 8)*lg), ' ns'
    deallocate(a)
  end subroutine measure

  ! Seconds per forward+inverse pair, one run.
  real(8) function time_pair(a, nx, ny)
    integer,    intent(in)    :: nx, ny
    complex(8), intent(inout) :: a(0:nx-1, 0:ny-1)
    integer(8) :: reps
    real(8)    :: c0, c1, el
    reps = 0
    call cpu_time(c0)
    do
      call kh_fft_forward_2d(a, nx, ny)
      call kh_fft_inverse_2d(a, nx, ny)
      reps = reps + 1
      call cpu_time(c1)
      el = c1 - c0
      if (el >= MIN_RUN_S) exit
    end do
    time_pair = el / real(reps, 8)
  end function time_pair

  real(8) function median(v)
    real(8), intent(in) :: v(:)
    real(8) :: w(size(v)), t
    integer :: i, j
    w = v
    do i = 2, size(w)
      t = w(i); j = i - 1
      do while (j >= 1)
        if (w(j) <= t) exit
        w(j+1) = w(j); j = j - 1
      end do
      w(j+1) = t
    end do
    median = w(size(w)/2 + 1)
  end function median

  real(8) function slope(x, y, n)
    real(8), intent(in) :: x(:), y(:)
    integer, intent(in) :: n
    real(8) :: sx, sy, sxx, sxy
    sx = sum(x(1:n)); sy = sum(y(1:n))
    sxx = sum(x(1:n)**2); sxy = sum(x(1:n)*y(1:n))
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
  end function slope

end program test_complexity_c4
