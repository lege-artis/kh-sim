unit kh_forcing;
{ kh-sim phase-7f forcing-field construction (Pascal).
  Port of backends/fortran/src/kh_forcing.f90 (the verified reference). The
  deterministic phase math (sqrt(2),sqrt(3) seeds; sqrt(11),sqrt(13) rotation
  rates) is RNG-free and reproduces bit-near-identically in every language --
  that portability is the whole point, so this is translated 1:1 from Fortran.

  Refs:
    PHASE-7F-FORCING-API-DESIGN-v0.1.md sec.3 (mode/band) + sec.7 (v0.2 rotating)
    Lesieur (2008) Turbulence in Fluids 4th ed ch.9 sec.9.3   -- CP-7F-B context
    Brachet et al. (1983) J. Fluid Mech. 130:411-452          -- CP-7F-D anchor

  Array layout: flat row-major, field[i*ny + j], (i,j)=(x,y) -- SAME layout as
  omega in kh_physics.pas. Grid: x_i = i*lx/nx (i=0..nx-1), y_j = j*ly/ny.

  SCOPE (2026-06-16; authored before the native FPC build):
    This unit provides the pure forcing-FIELD builders only. The forced RK4
    integrator (Fortran's rk4_step_forced) needs the solver's FFT/Poisson/
    velocity internals, which kh_physics.pas does NOT export. The brief says add
    forcing "like scala/rust did" (thread it through the Solver). So on the FPC
    pass, choose ONE and do it under an FPC compile + the MANDATORY phase-6c/7d
    bit-identity regression before any green tag:
      (A) extend kh_physics.pas: add forcing_* fields to TSimRequest (default
          'none'/0 => byte-identical pre-7f path) + apply (build via this unit)
          + linear drag in the RK4 RHS, guarded; OR
      (B) export the FFT2D/SolvePoisson/VelFromPsi/VortRHS internals from
          kh_physics.pas and add rk4_step_forced here (mirrors Fortran exactly).
    (A) matches scala/rust + the design doc schema; (B) matches the Fortran
    "don't touch the core solver" rule. Recommend (A) for cross-backend parity.

  ASCII-only source (KB-039). }
{$mode objfpc}{$H+}

interface

uses
  Math;

type
  TForceReal = array of Double;
  TForceInt  = array of Integer;

{ CP-7F-D steady Taylor-Green-mode forcing:
  F(x,y) = A*cos(2*pi*kx_mode*x/lx)*cos(2*pi*ky_mode*y/ly).
  Returns a flat (nx*ny) field. Ref: design sec.3 "mode"; Brachet (1983). }
function BuildModeForcing(nx, ny: Integer; lx, ly, amp: Double;
  kx_mode, ky_mode: Integer): TForceReal;

{ Deterministic band-mode phase seeds:
  phi0(ix,iy) = 2*pi*frac(ix*sqrt(2) + iy*sqrt(3)) over the shell
  kmin <= sqrt(ix^2+iy^2) <= kmax, half-plane ix>=0 skipping (ix=0, iy<=0).
  Fills ix_arr/iy_arr/phi0_arr (resized); returns the mode count M.
  Ref: design sec.3 "band". }
function BuildBandPhi0(kmin, kmax: Double;
  out ix_arr, iy_arr: TForceInt; out phi0_arr: TForceReal): Integer;

{ Per-mode rotation rates for v0.2 rotating-phase forcing:
  Omega(ix,iy) = 2*pi*frac(ix*sqrt(11) + iy*sqrt(13))*forcing_omega0.
  forcing_omega0 = 0 reproduces v0.1 (time-independent). Ref: design sec.7. }
function BuildBandOmegas(const ix_arr, iy_arr: TForceInt;
  forcing_omega0: Double): TForceReal;

{ Evaluate the band forcing field at time t_now:
  F = (A/sqrt(M)) * sum_m cos(2*pi*(ix*x/lx + iy*y/ly) + phi0_m + Omega_m*t).
  omg_arr all-zero reproduces v0.1 (time-independent). Returns flat (nx*ny).
  Ref: design sec.3 + sec.7. }
function EvalBandForcing(nx, ny: Integer; lx, ly, amp: Double;
  const ix_arr, iy_arr: TForceInt; const phi0_arr, omg_arr: TForceReal;
  t_now: Double): TForceReal;

implementation

const
  TWO_PI = 6.28318530717958647692; { 2*pi, matches Fortran 8*atan(1) }

function BuildModeForcing(nx, ny: Integer; lx, ly, amp: Double;
  kx_mode, ky_mode: Integer): TForceReal;
var
  i, j: Integer;
  xi, yj: Double;
begin
  SetLength(Result, nx * ny);
  for i := 0 to nx - 1 do
  begin
    xi := i * lx / nx;
    for j := 0 to ny - 1 do
    begin
      yj := j * ly / ny;
      Result[i * ny + j] := amp * Cos(TWO_PI * kx_mode * xi / lx)
                                * Cos(TWO_PI * ky_mode * yj / ly);
    end;
  end;
end;

function BuildBandPhi0(kmin, kmax: Double;
  out ix_arr, iy_arr: TForceInt; out phi0_arr: TForceReal): Integer;
var
  ix_max, ix, iy, mm, cap: Integer;
  sqrt2, sqrt3, kr, fv: Double;
begin
  sqrt2  := Sqrt(2.0);
  sqrt3  := Sqrt(3.0);
  ix_max := Trunc(kmax) + 1;
  cap    := (2 * ix_max + 2) * (ix_max + 2);  { generous upper bound }
  SetLength(ix_arr, cap);
  SetLength(iy_arr, cap);
  SetLength(phi0_arr, cap);
  mm := 0;
  for ix := 0 to ix_max do
    for iy := -ix_max to ix_max do
    begin
      if (ix = 0) and (iy <= 0) then Continue;
      kr := Sqrt(Sqr(Double(ix)) + Sqr(Double(iy)));
      if (kr < kmin) or (kr > kmax) then Continue;
      ix_arr[mm] := ix;
      iy_arr[mm] := iy;
      fv := ix * sqrt2 + iy * sqrt3;
      fv := fv - Floor(fv);
      phi0_arr[mm] := TWO_PI * fv;
      Inc(mm);
    end;
  SetLength(ix_arr, mm);
  SetLength(iy_arr, mm);
  SetLength(phi0_arr, mm);
  Result := mm;
end;

function BuildBandOmegas(const ix_arr, iy_arr: TForceInt;
  forcing_omega0: Double): TForceReal;
var
  m, nmode: Integer;
  s11, s13, fv: Double;
begin
  nmode := Length(ix_arr);
  SetLength(Result, nmode);
  s11 := Sqrt(11.0);
  s13 := Sqrt(13.0);
  for m := 0 to nmode - 1 do
  begin
    fv := ix_arr[m] * s11 + iy_arr[m] * s13;
    fv := fv - Floor(fv);
    Result[m] := TWO_PI * fv * forcing_omega0;
  end;
end;

function EvalBandForcing(nx, ny: Integer; lx, ly, amp: Double;
  const ix_arr, iy_arr: TForceInt; const phi0_arr, omg_arr: TForceReal;
  t_now: Double): TForceReal;
var
  i, j, m, nmode: Integer;
  scale, xi, yj, phase: Double;
begin
  nmode := Length(ix_arr);
  SetLength(Result, nx * ny);
  FillChar(Result[0], nx * ny * SizeOf(Double), 0);
  if nmode < 1 then Exit;
  scale := amp / Sqrt(Double(nmode));
  for m := 0 to nmode - 1 do
    for i := 0 to nx - 1 do
    begin
      xi := i * lx / nx;
      for j := 0 to ny - 1 do
      begin
        yj := j * ly / ny;
        phase := TWO_PI * (ix_arr[m] * xi / lx + iy_arr[m] * yj / ly)
               + phi0_arr[m] + omg_arr[m] * t_now;
        Result[i * ny + j] := Result[i * ny + j] + scale * Cos(phase);
      end;
    end;
end;

end.
