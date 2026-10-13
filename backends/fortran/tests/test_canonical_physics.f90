! test_canonical_physics.f90 -- phase 7f canonical-physics testbeds (Fortran)
! Mode 2: in-process, no HTTP server needed.
!
! Tests: CP-7F-E (Numerical Recipes in Pascal: E2 s17.5 Dirichlet-Poisson, E1 s12.4 convlv)
!          -- asserts here (see run_cp7fe), as in the Pascal suite (MASTER s3.5 E8).
!        CP-7F-A, B, C, D run from the shared harnesses in tests/integration (ruling 2):
!          A test_canonical_physics_cp7fa.py (E9, moved 2026-10-12), B and D
!          test_canonical_physics_cp7fbd.py (E7), C test_canonical_physics_cp7fc.py (E6);
!          the modes A, Afinal, B, C, Cfinal, D print SKIP with that pointer.
!   Invoke with: ./test_canonical_physics [all|A|B|C|D|E]
!
! Compilation:
!   gfortran -O2 -J/tmp/kh_mods -I/tmp/kh_mods -o test_canonical_physics \
!       src/kh_constants.f90 src/kh_fft.f90 src/kh_grid.f90 \
!       src/kh_initial_conditions.f90 src/kh_nonlinear.f90 \
!       src/kh_poisson.f90 src/kh_velocity.f90 \
!       src/kh_diagnostics.f90 src/kh_forcing.f90 src/kh_physics.f90 \
!       tests/test_canonical_physics.f90

program test_canonical_physics
  use kh_constants
  use kh_fft,                only: kh_fft_forward_2d, kh_fft_inverse_2d
  implicit none

  character(len=16) :: arg1, arg2
  integer           :: argc
  logical           :: g_fail = .false.   ! set by a failed CP-7F-E check -> exit 1

  ! ---- parse arguments ----
  argc = command_argument_count()
  arg1 = ''
  arg2 = ''
  if (argc >= 1) call get_command_argument(1, arg1)
  if (argc >= 2) call get_command_argument(2, arg2)

  if (arg1 == '' .or. arg1 == 'all') then
    ! Run everything in sequence (useful for short CI runs)
    call cp7fa_moved()
    call cp7fbd_moved('B')
    call cp7fc_moved()
    call cp7fbd_moved('D')
    call run_cp7fe()
  else if (arg1 == 'A' .or. arg1 == 'Afinal') then
    call cp7fa_moved()
  else if (arg1 == 'C' .or. arg1 == 'Cfinal') then
    call cp7fc_moved()
  else if (arg1 == 'B' .or. arg1 == 'D') then
    call cp7fbd_moved(trim(arg1))
  else if (arg1 == 'E') then
    call run_cp7fe()
  else
    write(*,'(A,A)') 'Unknown mode: ', trim(arg1)
    stop 1
  end if
  if (g_fail) stop 1

contains

  ! ===========================================================================
  ! CP-7F-B and CP-7F-D: MOVED 2026-10-10 to
  ! tests/integration/test_canonical_physics_cp7fbd.py (MASTER-BRIEF s3.2/s3.4,
  ! E7), which runs one specification each against the oracle and every backend
  ! over HTTP (ruling 8: the Fortran server computes forcing). The B stub
  ! (BLOCKED-pending-oracle) and the in-process D at Re 500..5000 that stood here
  ! are in the git history of this file.
  ! ===========================================================================
  subroutine cp7fbd_moved(which)
    character(len=*), intent(in) :: which
    if (which == 'B') then
      write(*,'(A)') 'CP-7F-B: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.2 E7)'
    else
      write(*,'(A)') 'CP-7F-D: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.4 E7)'
    end if
  end subroutine cp7fbd_moved

  ! ===========================================================================
  ! CP-7F-E: Numerical Recipes in Pascal -- the Pascal-binding strict-C1 anchor
  !          (MASTER s3.5), cross-checked here as for Pascal. Answer key: the
  !          oracle's test_cp7fe_nr1986_ch19_poisson / test_cp7fe_nr1986_convlv.
  !
  ! E2  Press, Flannery, Teukolsky, Vetterling, Numerical Recipes in Pascal
  !     (CUP 1989) s17.5 + Example Book (Pascal) d17r1, pp.217-219: SOR
  !     Dirichlet-Poisson, 11x11 lattice, u = 0 on the boundary, source 2.0 at
  !     the centre node, printed grid (2 decimals). Solved by the Fourier method
  !     for Dirichlet boundaries (NR s17.4): orthonormal sine transform on both
  !     axes, division by the 5-point eigenvalues, inverse transform. Gates: the
  !     printed points at tolerance_rel (2e-2, OQ-7F-PAS-01), the residual at
  !     the source = 2 (1e-9), the 5-point equation on every node (1e-12).
  ! E1  Same edition s12.4 + Example Book d12r6, pp.154-155: convlv, data[6..10]
  !     = 1 (n = 16) with respns[3..6] = 1 (m = 9, NR wrap-around lags), through
  !     the production kh_fft (IFFT(FFT(data) * FFT(h))). Gate 1e-6 per channel
  !     (exact integers), sum = 20.
  !
  ! The namelists are generated from the golden vectors by
  ! tools/nr1986_json_to_namelist.py (nr1986_pascal_ch19.nml, nr1986_convlv.nml)
  ! and found through KH_FORTRAN_TESTS_DIR, else tests/, ../tests/ or ./ . A
  ! missing namelist is a FAIL (the transcription exists since 2026-06-17).
  ! ===========================================================================
  subroutine run_cp7fe()
    write(*,'(A)') ''
    write(*,'(A)') 'CP-7F-E: Numerical Recipes in Pascal, Pascal-binding strict-C1 anchor (MASTER s3.5)'
    call cp7fe_e2_poisson()
    call cp7fe_e1_convlv()
    if (g_fail) then
      write(*,'(A)') 'CP-7F-E: FAIL'
    else
      write(*,'(A)') 'CP-7F-E: PASS (E2 + E1)'
    end if
  end subroutine run_cp7fe

  subroutine fe_check(what, pass, detail)
    character(len=*), intent(in) :: what, detail
    logical,          intent(in) :: pass
    if (pass) then
      write(*,'(A,A,A,A)') '  PASS  ', what, '  ', detail
    else
      write(*,'(A,A,A,A)') '  FAIL  ', what, '  ', detail
      g_fail = .true.
    end if
  end subroutine fe_check

  ! Locate a namelist; '' when absent.
  function fe_find(fname) result(path)
    character(len=*), intent(in) :: fname
    character(len=512) :: path
    character(len=512) :: d
    character(len=12), parameter :: cands(3) = [character(len=12) :: 'tests/', '../tests/', './']
    integer :: k, l
    logical :: there
    call get_environment_variable('KH_FORTRAN_TESTS_DIR', d, length=l)
    if (l > 0) then
      path = trim(d) // '/' // fname
      inquire(file=trim(path), exist=there)
      if (there) return
    end if
    do k = 1, size(cands)
      path = trim(cands(k)) // fname
      inquire(file=trim(path), exist=there)
      if (there) return
    end do
    path = ''
  end function fe_find

  ! Read /NR1986_CP7FE/ from one namelist file; ok = .false. on any problem.
  subroutine fe_read(fname, ok, gnx, gny, ns, si, sj, sv, iname, ival, tol, src_i, src_j, src_val, m_resp)
    integer, parameter :: MAXS = 64   ! keep in sync with generator MAX_SAMPLES
    character(len=*), intent(in)  :: fname
    logical,          intent(out) :: ok
    integer,          intent(out) :: gnx, gny, ns, src_i, src_j, m_resp
    integer,          intent(out) :: si(MAXS), sj(MAXS)
    real(c_double),   intent(out) :: sv(MAXS), ival, tol, src_val
    character(len=64), intent(out) :: iname

    ! namelist /NR1986_CP7FE/ -- names MUST match tools/nr1986_json_to_namelist.py
    logical            :: data_present
    character(len=64)  :: chapter, section, prob_type, boundary, integ_name
    character(len=128) :: example_name
    integer            :: grid_nx, grid_ny, n_samples
    integer            :: sample_i(MAXS), sample_j(MAXS)
    real(c_double)     :: sample_val(MAXS), integ_val, tolerance_rel
    namelist /NR1986_CP7FE/ data_present, chapter, section, example_name,      &
         prob_type, boundary, grid_nx, grid_ny, n_samples,                     &
         sample_i, sample_j, sample_val, integ_name, integ_val, tolerance_rel, &
         src_i, src_j, src_val, m_resp
    character(len=512) :: path
    integer :: u, ios

    ok = .false.
    path = fe_find(fname)
    if (len_trim(path) == 0) then
      call fe_check(fname, .false., 'namelist not found (KH_FORTRAN_TESTS_DIR, tests/, ../tests/, ./);' // &
                    ' make it with tools/nr1986_json_to_namelist.py')
      return
    end if
    data_present = .false.; n_samples = 0; tolerance_rel = 1.0e-3_c_double
    integ_name = ''; integ_val = 0.0_c_double
    sample_i = 0; sample_j = 0; sample_val = 0.0_c_double
    src_i = 0; src_j = 0; src_val = 0.0_c_double; m_resp = 0
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios == 0) then
      read(u, nml=NR1986_CP7FE, iostat=ios)
      close(u)
    end if
    if (ios /= 0 .or. .not. data_present .or. n_samples < 1 .or. n_samples > MAXS) then
      call fe_check(trim(path), .false., 'namelist unreadable, data_present=.false. or bad n_samples')
      return
    end if
    write(*,'(A,A,A,A,A,A,A,A)') '  ', trim(path), ': NR in Pascal s', trim(section), '  [', &
         trim(example_name), ']  ', trim(prob_type)
    gnx = grid_nx; gny = grid_ny; ns = n_samples
    si = sample_i; sj = sample_j; sv = sample_val
    iname = integ_name; ival = integ_val; tol = tolerance_rel
    ok = .true.
  end subroutine fe_read

  subroutine cp7fe_e2_poisson()
    integer, parameter :: MAXS = 64
    logical :: ok
    integer :: gnx, gny, ns, si(MAXS), sj(MAXS), src_i, src_j, m_resp, n, i, j, k, s
    real(c_double) :: sv(MAXS), ival, tol, src_val, c, rel, worst, resid, max_res, f, book, got
    character(len=64)  :: iname
    character(len=160) :: line
    real(c_double), allocatable :: sm(:,:), fm(:,:), tm(:,:), g(:,:)

    write(*,'(A)') 'CP-7F-E/E2: NR in Pascal s17.5, d17r1 SOR Dirichlet-Poisson'
    call fe_read('nr1986_pascal_ch19.nml', ok, gnx, gny, ns, si, sj, sv, iname, ival, tol, &
                 src_i, src_j, src_val, m_resp)
    if (.not. ok) return
    if (gnx /= gny .or. gnx < 3 .or. src_i < 2 .or. src_j < 2) then
      call fe_check('grid', .false., 'the sine-transform solve needs a square lattice and an interior source')
      return
    end if
    n = gnx - 2
    allocate(sm(n,n), fm(n,n), tm(n,n), g(gnx,gny))
    ! orthonormal DST-I: sm(k,i) = sqrt(2/(n+1)) sin(pi k i/(n+1)); sm.sm = I
    c = sqrt(2.0_c_double / real(n + 1, c_double))
    do k = 1, n
      do i = 1, n
        sm(k,i) = c * sin(KH_PI * real(k*i, c_double) / real(n + 1, c_double))
      end do
    end do
    fm = 0.0_c_double
    fm(src_i - 1, src_j - 1) = src_val          ! interior index = lattice index - 1
    tm = matmul(sm, matmul(fm, sm))
    do i = 1, n
      do j = 1, n
        tm(i,j) = tm(i,j) / (2.0_c_double * cos(real(i, c_double) * KH_PI / real(n + 1, c_double)) + &
                             2.0_c_double * cos(real(j, c_double) * KH_PI / real(n + 1, c_double)) - 4.0_c_double)
      end do
    end do
    g = 0.0_c_double
    g(2:n+1, 2:n+1) = matmul(sm, matmul(tm, sm))  ! lattice, 1-based (book indices)

    max_res = 0.0_c_double
    do i = 2, gnx - 1
      do j = 2, gny - 1
        f = 0.0_c_double
        if (i == src_i .and. j == src_j) f = src_val
        resid = g(i+1,j) + g(i-1,j) + g(i,j+1) + g(i,j-1) - 4.0_c_double * g(i,j)
        max_res = max(max_res, abs(resid - f))
      end do
    end do
    write(line,'(A,ES9.2)') 'max |residual - f| = ', max_res
    call fe_check('5-point equation on all interior nodes', max_res <= 1.0e-12_c_double, trim(line))

    worst = 0.0_c_double
    do s = 1, ns
      book = sv(s); got = g(si(s), sj(s))
      if (abs(book) > 1.0e-30_c_double) then
        rel = abs(got - book) / abs(book)
      else
        rel = abs(got - book)
      end if
      worst = max(worst, rel)
      write(line,'(A,F6.2,A,F10.6,A,ES9.3,A,ES8.1,A)') 'book ', book, '  got ', got, '  rel ', rel, &
           '  (gate ', tol, ')'
      write(iname,'(A,I0,A,I0,A)') 'u[', si(s), ',', sj(s), ']'
      call fe_check(trim(iname), rel <= tol, trim(line))
    end do
    resid = g(src_i+1,src_j) + g(src_i-1,src_j) + g(src_i,src_j+1) + g(src_i,src_j-1) - 4.0_c_double * g(src_i,src_j)
    write(line,'(F15.12,A,F4.1)') resid, ' vs ', ival
    call fe_check('laplacian_residual_at_center', abs(resid - ival) <= 1.0e-9_c_double, trim(line))
    write(*,'(A,I0,A,ES10.3)') '  E2: ', ns, ' printed points, worst rel-err ', worst
  end subroutine cp7fe_e2_poisson

  subroutine cp7fe_e1_convlv()
    integer, parameter :: MAXS = 64
    logical :: ok
    integer :: gnx, gny, ns, si(MAXS), sj(MAXS), src_i, src_j, m_resp, n, i, r, lag, s
    real(c_double) :: sv(MAXS), ival, tol, src_val, total, worst, book, got
    character(len=64)  :: iname
    character(len=160) :: line
    complex(8), allocatable :: a(:,:), h(:,:)   ! the kind kh_fft takes

    write(*,'(A)') 'CP-7F-E/E1: NR in Pascal s12.4, d12r6 convlv'
    call fe_read('nr1986_convlv.nml', ok, gnx, gny, ns, si, sj, sv, iname, ival, tol, &
                 src_i, src_j, src_val, m_resp)
    if (.not. ok) return
    n = gnx
    tol = max(tol, 1.0e-9_c_double)
    allocate(a(0:n-1, 0:0), h(0:n-1, 0:0))
    a = (0.0_c_double, 0.0_c_double); h = (0.0_c_double, 0.0_c_double)
    ! data[i] = 1 for n/2-n/8 <= i <= n/2+n/8 (1-based; 6..10 at n = 16)
    do i = n/2 - n/8, n/2 + n/8
      a(i-1, 0) = (1.0_c_double, 0.0_c_double)
    end do
    ! respns[3..6] = 1 over m; NR stores lag r-1 for r-1 <= m/2, else r-1-m
    do r = 3, 6
      if (r - 1 <= m_resp / 2) then
        lag = r - 1
      else
        lag = r - 1 - m_resp
      end if
      h(modulo(lag, n), 0) = (1.0_c_double, 0.0_c_double)
    end do
    call kh_fft_forward_2d(a, n, 1)
    call kh_fft_forward_2d(h, n, 1)
    a = a * h
    call kh_fft_inverse_2d(a, n, 1)

    total = sum(real(a(:,0), c_double))
    worst = 0.0_c_double
    do s = 1, ns
      book = sv(s); got = real(a(si(s) - 1, 0), c_double)
      worst = max(worst, abs(got - book))
      write(line,'(A,F4.1,A,F18.15)') 'book ', book, '  got ', got
      write(iname,'(A,I0,A)') 'ans[', si(s), ']'
      call fe_check(trim(iname), abs(got - book) <= tol, trim(line))
    end do
    write(line,'(F15.12,A,F4.1)') total, ' vs ', ival
    call fe_check('sum_ans', abs(total - ival) <= tol, trim(line))
    write(*,'(A,I0,A,ES9.2,A,ES8.1,A)') '  E1: ', ns, ' channels, worst abs-err ', worst, ' (gate ', tol, ')'
  end subroutine cp7fe_e1_convlv

  ! ===========================================================================
  ! CP-7F-A: MOVED 2026-10-12 to tests/integration/test_canonical_physics_cp7fa.py
  ! (MASTER-BRIEF s3.1, E9), which runs one specification for the oracle and every
  ! backend. The per-mode runs that stood here (12 seeded modes on the double layer at
  ! delta = Ly/20, where the two layers couple - assessment 2026-09-15, s7 - with the
  ! gate [0.180, 0.200]) are in the git history of this file.
  ! ===========================================================================
  subroutine cp7fa_moved()
    write(*,'(A)') 'CP-7F-A: SKIP -- moved to tests/integration/test_canonical_physics_cp7fa.py (MASTER s3.1 E9)'
  end subroutine cp7fa_moved

  ! ===========================================================================
  ! CP-7F-C: MOVED 2026-10-07 to tests/integration/test_canonical_physics_cp7fc.py
  ! (MASTER-BRIEF s3.3, E6), which runs one specification against every backend.
  ! The two-maxima tracker that stood here never fired on a merged vortex, and its
  ! nu = a^2/Re made Re_Gamma ~ 1840 instead of 5000 (a 3.5a pair then merges);
  ! see _audit/FINDING-CP7FC-MERGE-DETECTION-2026-10-06.md.
  ! ===========================================================================
  subroutine cp7fc_moved()
    write(*,'(A)') 'CP-7F-C: SKIP -- moved to tests/integration/test_canonical_physics_cp7fc.py (MASTER s3.3 E6)'
  end subroutine cp7fc_moved

end program test_canonical_physics
