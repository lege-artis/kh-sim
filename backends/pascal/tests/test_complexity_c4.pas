program test_complexity_c4;
{ kh-sim Pascal backend - C4 complexity-claim audit (LEGE-ARTIS language doctrine
  v0.1 s4.4 C4) for the hand-rolled FFT.

  Claim under audit: FFT1D (src/kh_physics.pas) is O(N log N) - its C3 header,
  after Press et al. (1986), Numerical Recipes in Pascal S12.2.

  Procedure (doctrine C4, literally): N = 2^6 .. 2^16; 5 runs per size, the median
  is used; execution time t(N) of one transform; least-squares slope of ln t over
  ln(N log2 N). PASS when 0.9 <= slope <= 1.1. An O(N^2) defect (for example a
  naive bit reversal) gives a slope near 2 over this range, an O(N^2 log N) one
  higher still.

  One unit of work is a forward transform followed by the inverse one (which
  normalises), so the data stay bounded however often the pair is repeated; its
  cost is 2 transforms + N divisions, O(N log N). Each run repeats the pair until
  at least MIN_RUN_S seconds have elapsed, so the clock resolution is irrelevant.
  Execution time is process CPU time (clock_gettime CLOCK_PROCESS_CPUTIME_ID) on
  Linux, so preemption by other processes on a shared runner does not enter the
  fit; elsewhere it falls back to wall time (SysUtils.Now).

  Build (release, no checks - C4 measures run time; doctrine: release-no-LTO):
    fpc -O2 -Fusrc -FUbuild/c4 -FEbuild/c4 tests/test_complexity_c4.pas
  Run:   ./build/c4/test_complexity_c4      exit 0 = PASS
  Output: one line per size, then "C4 FFT1D slope = s" and RESULT: PASS/FAIL.
  ASCII-only source. }
{$mode objfpc}{$H+}

uses SysUtils, DateUtils, Math, {$IFDEF LINUX} Linux, UnixType, {$ENDIF} kh_physics;

const
  NMIN_EXP  = 6;
  NMAX_EXP  = 16;
  RUNS      = 5;
  MIN_RUN_S = 0.10;
  SLOPE_LO  = 0.9;
  SLOPE_HI  = 1.1;

type
  TDoubles = array of Double;

function Median(v: TDoubles): Double;
var i, j: Integer; t: Double;
begin
  for i := 1 to High(v) do begin
    t := v[i]; j := i - 1;
    while (j >= 0) and (v[j] > t) do begin v[j+1] := v[j]; Dec(j); end;
    v[j+1] := t;
  end;
  Result := v[Length(v) div 2];
end;

{ Execution-time clock in seconds. }
function ClockS: Double;
{$IFDEF LINUX}
var ts: TTimeSpec;
begin
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, @ts);
  Result := ts.tv_sec + ts.tv_nsec * 1e-9;
end;
{$ELSE}
begin
  Result := MilliSecondsBetween(Now, 0) / 1000.0;
end;
{$ENDIF}

{ Seconds per forward+inverse pair, one run. }
function TimePair(var a: TComplexArr; n: Integer): Double;
var reps: Int64; t0, el: Double;
begin
  reps := 0;
  t0 := ClockS;
  repeat
    FFT1D(a, n, False);
    FFT1D(a, n, True);
    Inc(reps);
    el := ClockS - t0;
  until el >= MIN_RUN_S;
  Result := el / reps;
end;

var
  e, i, r, n, cnt: Integer;
  a: TComplexArr;
  samples: TDoubles;
  xs, ys: array of Double;
  sx, sy, sxx, sxy, slope, t: Double;
  ok: Boolean;
begin
  WriteLn('C4 complexity audit - FFT1D (kh_physics.pas), claim O(N log N)');
  WriteLn(Format('  N = 2^%d..2^%d, %d runs per size (median), >= %.2f s CPU per run',
                 [NMIN_EXP, NMAX_EXP, RUNS, MIN_RUN_S]));
  cnt := NMAX_EXP - NMIN_EXP + 1;
  SetLength(xs, cnt); SetLength(ys, cnt); SetLength(samples, RUNS);
  for e := NMIN_EXP to NMAX_EXP do begin
    n := 1 shl e;
    SetLength(a, n);
    for i := 0 to n-1 do a[i] := CmplxOf(Sin(0.37*i), Cos(0.11*i));
    FFT1D(a, n, False); FFT1D(a, n, True);          { warm-up }
    for r := 0 to RUNS-1 do samples[r] := TimePair(a, n);
    t := Median(Copy(samples, 0, RUNS));
    xs[e-NMIN_EXP] := Ln(Double(n) * e);           { N log2 N }
    ys[e-NMIN_EXP] := Ln(t);
    WriteLn(Format('  N = %6d   t = %10.3f us   t/(N log2 N) = %8.3f ns',
                   [n, t*1e6, t*1e9/(Double(n)*e)]));
  end;
  sx := 0; sy := 0; sxx := 0; sxy := 0;
  for i := 0 to cnt-1 do begin
    sx := sx + xs[i]; sy := sy + ys[i];
    sxx := sxx + xs[i]*xs[i]; sxy := sxy + xs[i]*ys[i];
  end;
  slope := (cnt*sxy - sx*sy) / (cnt*sxx - sx*sx);
  ok := (slope >= SLOPE_LO) and (slope <= SLOPE_HI);
  WriteLn(Format('C4 FFT1D slope = %.3f  (gate [%.1f, %.1f])', [slope, SLOPE_LO, SLOPE_HI]));
  if ok then begin WriteLn('RESULT: PASS'); Halt(0); end
  else begin WriteLn('RESULT: FAIL'); Halt(1); end;
end.
