unit kh_cp7fe;
{ kh-sim phase 7f -- CP-7F-E, the Pascal-binding strict-C1 anchor (MASTER s3.5).

  Two sub-anchors from Numerical Recipes in Pascal, read from the golden vectors
  the oracle reads (tests/integration/test_canonical_physics_python.py is the
  answer key; framing record CP7F-E-TRANSCRIPTION-DECISION-2026-06-17.md s8):

  E2  Press, Flannery, Teukolsky, Vetterling, Numerical Recipes in Pascal (CUP
      1989) s17.5 + Example Book (Pascal) driver d17r1, pp.217-219: the SOR
      Dirichlet-Poisson model problem, 11x11 lattice, u = 0 on the boundary,
      source 2.0 at the centre node, printed solution grid (2 decimals).
      Solved here by the Fourier method for Dirichlet boundaries (NR s17.4):
      an orthonormal sine transform (DST-I) on both axes, division by the
      5-point eigenvalues, the inverse transform. Gate: the JSON's
      tolerance_rel (2e-2, printed-precision exception OQ-7F-PAS-01), the
      residual at the source node = 2 (1e-9), and the 5-point equation on
      every interior node (1e-12, the discrete solution itself).
  E1  Same edition s12.4 + Example Book (Pascal) driver d12r6, pp.154-155:
      convlv, cyclic convolution of data[6..10] = 1 (n = 16) with respns[3..6]
      = 1 (m = 9, NR wrap-around lag convention). Computed with the backend's
      own FFT1D (kh_physics, the production radix-2 Cooley-Tukey): ans =
      IFFT(FFT(data) * FFT(h)). Gate: every channel within the JSON's
      tolerance (1e-6, exact integers), sum = 20.

  The golden vectors are found through KH_GOLDEN_DIR, else by walking up from
  the working directory and from the executable's directory to
  shared/golden-vectors (or kh-sim/shared/golden-vectors). A missing file is a
  FAIL: the transcription exists since 2026-06-17, so a SKIP would hide it.

  ASCII-only source (KB-039). }
{$mode objfpc}{$H+}

interface

{ Runs E2 and E1; prints one line per check; True when every check passed. }
function RunCp7fE(const goldenDirArg: string): Boolean;

implementation

uses
  SysUtils, Math, Classes, fpjson, jsonparser, kh_physics;

const
  CH19_FILE   = 'nr1986-pascal-ch19.json';
  CONVLV_FILE = 'nr1986-pascal-ch12-convlv.json';

var
  ok: Boolean;

procedure Check(const what: string; pass: Boolean; const detail: string);
begin
  if pass then WriteLn('  PASS  ', what, '  ', detail)
  else begin WriteLn('  FAIL  ', what, '  ', detail); ok := False; end;
end;

function FindGoldenDir(const arg: string): string;
var
  starts: array[0..1] of string;
  s, d, cand: string;
  up, k: Integer;
begin
  if arg <> '' then exit(IncludeTrailingPathDelimiter(arg));
  d := GetEnvironmentVariable('KH_GOLDEN_DIR');
  if d <> '' then exit(IncludeTrailingPathDelimiter(d));
  starts[0] := GetCurrentDir;
  starts[1] := ExtractFileDir(ExpandFileName(ParamStr(0)));
  for k := 0 to 1 do
  begin
    s := ExcludeTrailingPathDelimiter(starts[k]);
    for up := 0 to 6 do
    begin
      cand := IncludeTrailingPathDelimiter(s) + 'shared' + PathDelim + 'golden-vectors' + PathDelim;
      if FileExists(cand + CH19_FILE) then exit(cand);
      cand := IncludeTrailingPathDelimiter(s) + 'kh-sim' + PathDelim + 'shared' + PathDelim
              + 'golden-vectors' + PathDelim;
      if FileExists(cand + CH19_FILE) then exit(cand);
      d := ExtractFileDir(s);
      if (d = '') or (d = s) then break;
      s := d;
    end;
  end;
  Result := '';
end;

function LoadJson(const path: string): TJSONObject;
var
  sl: TStringList;
begin
  sl := TStringList.Create;
  try
    sl.LoadFromFile(path);
    Result := GetJSON(sl.Text) as TJSONObject;
  finally
    sl.Free;
  end;
end;

{ ---- E2: Dirichlet-Poisson by the sine transform (NR s17.4) ---------------- }
procedure RunE2(const path: string);
type
  TMat = array of array of Double;
var
  root, prob, params, exp, sp, it: TJSONObject;
  samples, integrals: TJSONArray;
  nx, ny, n, i, j, k, si, sj, s: Integer;
  srcVal, tol, c, got, book, rel, worst, resid, maxRes, f: Double;
  Sm, F2, T, W: TMat;
  g: TMat;   { full lattice, 1-based (Pascal / book indices), 0 on the boundary }

  function Mul(const A, B: TMat; m: Integer): TMat;
  var r, q, p: Integer; acc: Double;
  begin
    SetLength(Result, m + 1, m + 1);
    for r := 1 to m do
      for q := 1 to m do
      begin
        acc := 0.0;
        for p := 1 to m do acc := acc + A[r][p] * B[p][q];
        Result[r][q] := acc;
      end;
  end;

begin
  WriteLn('CP-7F-E/E2: NR in Pascal s17.5, d17r1 SOR Dirichlet-Poisson (', path, ')');
  root := LoadJson(path);
  try
    prob := root.Objects['problem'];
    params := prob.Objects['params'];
    exp := root.Objects['expected'];
    nx := prob.Integers['grid_nx']; ny := prob.Integers['grid_ny'];
    if (nx <> ny) or (nx < 3) then
    begin
      Check('grid', False, Format('%dx%d: the sine-transform solve needs a square lattice', [nx, ny]));
      exit;
    end;
    si := params.Integers['source_i_1based']; sj := params.Integers['source_j_1based'];
    srcVal := params.Floats['source_value'];
    tol := exp.Floats['tolerance_rel'];
    n := nx - 2;

    { orthonormal DST-I: Sm[k][i] = sqrt(2/(n+1)) sin(pi k i/(n+1)); Sm.Sm = I }
    SetLength(Sm, n + 1, n + 1);
    { Double(...): FPC types the literal 2.0 as Single, which cost 7 digits here }
    c := Sqrt(Double(2.0) / (n + 1));
    for k := 1 to n do
      for i := 1 to n do
        Sm[k][i] := c * Sin(Pi * k * i / (n + 1));
    { source on the interior (interior index = lattice index - 1) }
    SetLength(F2, n + 1, n + 1);
    for i := 1 to n do for j := 1 to n do F2[i][j] := 0.0;
    F2[si - 1][sj - 1] := srcVal;
    T := Mul(Mul(Sm, F2, n), Sm, n);
    for i := 1 to n do
      for j := 1 to n do
        T[i][j] := T[i][j] / (2.0 * Cos(i * Pi / (n + 1)) + 2.0 * Cos(j * Pi / (n + 1)) - 4.0);
    W := Mul(Mul(Sm, T, n), Sm, n);
    SetLength(g, nx + 2, ny + 2);
    for i := 0 to nx + 1 do for j := 0 to ny + 1 do g[i][j] := 0.0;
    for i := 1 to n do for j := 1 to n do g[i + 1][j + 1] := W[i][j];

    { the discrete equation on every interior node }
    maxRes := 0.0;
    for i := 2 to nx - 1 do
      for j := 2 to ny - 1 do
      begin
        if (i = si) and (j = sj) then f := srcVal else f := 0.0;
        resid := g[i + 1][j] + g[i - 1][j] + g[i][j + 1] + g[i][j - 1] - 4.0 * g[i][j];
        maxRes := Max(maxRes, Abs(resid - f));
      end;
    Check('5-point equation on all interior nodes', maxRes <= 1e-12, Format('max |residual - f| = %.2e', [maxRes]));

    samples := exp.Arrays['sample_points'];
    worst := 0.0;
    for s := 0 to samples.Count - 1 do
    begin
      sp := samples.Objects[s];
      i := sp.Integers['i']; j := sp.Integers['j']; book := sp.Floats['value'];
      got := g[i][j];
      if Abs(book) > 1e-30 then rel := Abs(got - book) / Abs(book) else rel := Abs(got - book);
      worst := Max(worst, rel);
      Check(Format('u[%d,%d]', [i, j]), rel <= tol,
            Format('book %.2f  got %.6f  rel %.3e  (gate %.1e)', [book, got, rel, tol]));
    end;
    integrals := exp.Arrays['integrals'];
    for s := 0 to integrals.Count - 1 do
    begin
      it := integrals.Objects[s];
      if it.Strings['name'] = 'laplacian_residual_at_center' then
      begin
        resid := g[si + 1][sj] + g[si - 1][sj] + g[si][sj + 1] + g[si][sj - 1] - 4.0 * g[si][sj];
        Check('laplacian_residual_at_center', Abs(resid - it.Floats['value']) <= 1e-9,
              Format('%.12f vs %.1f', [resid, it.Floats['value']]));
      end;
    end;
    WriteLn(Format('  E2: %d printed points, worst rel-err %.3e', [samples.Count, worst]));
  finally
    root.Free;
  end;
end;

{ ---- E1: convlv through the production FFT --------------------------------- }
procedure RunE1(const path: string);
var
  root, prob, params, exp, sp, it: TJSONObject;
  samples, integrals: TJSONArray;
  n, m, i, r, lag, s, lo, hi: Integer;
  tol, got, book, sum, worst: Double;
  a, h: TComplexArr;
begin
  WriteLn('CP-7F-E/E1: NR in Pascal s12.4, d12r6 convlv (', path, ')');
  root := LoadJson(path);
  try
    prob := root.Objects['problem'];
    params := prob.Objects['params'];
    exp := root.Objects['expected'];
    n := prob.Integers['grid_nx'];
    m := params.Integers['m_response'];
    tol := Max(exp.Floats['tolerance_rel'], 1e-9);
    SetLength(a, n); SetLength(h, n);
    for i := 0 to n - 1 do begin a[i] := CmplxOf(0.0, 0.0); h[i] := CmplxOf(0.0, 0.0); end;
    { data[i] = 1 for n/2-n/8 <= i <= n/2+n/8 (1-based; 6..10 at n = 16) }
    lo := n div 2 - n div 8; hi := n div 2 + n div 8;
    for i := lo to hi do a[i - 1] := CmplxOf(1.0, 0.0);
    { respns[3..6] = 1 over m; NR stores lag r-1 for r-1 <= m/2, else r-1-m }
    for r := 3 to 6 do
    begin
      if (r - 1) <= m div 2 then lag := r - 1 else lag := r - 1 - m;
      h[((lag mod n) + n) mod n] := CmplxOf(1.0, 0.0);
    end;
    FFT1D(a, n, False);
    FFT1D(h, n, False);
    for i := 0 to n - 1 do a[i] := a[i] * h[i];
    FFT1D(a, n, True);

    samples := exp.Arrays['sample_points'];
    worst := 0.0; sum := 0.0;
    for i := 0 to n - 1 do sum := sum + a[i].re;
    for s := 0 to samples.Count - 1 do
    begin
      sp := samples.Objects[s];
      i := sp.Integers['i']; book := sp.Floats['value'];
      got := a[i - 1].re;
      worst := Max(worst, Abs(got - book));
      Check(Format('ans[%d]', [i]), Abs(got - book) <= tol,
            Format('book %.1f  got %.15f', [book, got]));
    end;
    integrals := exp.Arrays['integrals'];
    for s := 0 to integrals.Count - 1 do
    begin
      it := integrals.Objects[s];
      if it.Strings['name'] = 'sum_ans' then
        Check('sum_ans', Abs(sum - it.Floats['value']) <= tol, Format('%.12f vs %.1f', [sum, it.Floats['value']]));
    end;
    WriteLn(Format('  E1: %d channels, worst abs-err %.2e (gate %.1e)', [samples.Count, worst, tol]));
  finally
    root.Free;
  end;
end;

function RunCp7fE(const goldenDirArg: string): Boolean;
var
  dir: string;
begin
  ok := True;
  WriteLn('CP-7F-E: Numerical Recipes in Pascal, Pascal-binding strict-C1 anchor (MASTER s3.5)');
  dir := FindGoldenDir(goldenDirArg);
  if (dir = '') or not FileExists(dir + CH19_FILE) or not FileExists(dir + CONVLV_FILE) then
  begin
    WriteLn('  FAIL  golden vectors not found (set KH_GOLDEN_DIR to kh-sim/shared/golden-vectors);',
            ' looked for ', CH19_FILE, ' and ', CONVLV_FILE, ' from "', GetCurrentDir, '"');
    exit(False);
  end;
  try
    RunE2(dir + CH19_FILE);
    RunE1(dir + CONVLV_FILE);
  except
    on e: Exception do
    begin
      WriteLn('  FAIL  ', e.ClassName, ': ', e.Message);
      ok := False;
    end;
  end;
  if ok then WriteLn('CP-7F-E: PASS (E2 + E1)') else WriteLn('CP-7F-E: FAIL');
  Result := ok;
end;

end.
