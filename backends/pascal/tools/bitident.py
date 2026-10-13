"""Bit-identity regression: run every configuration on two kh_server builds and
compare every field array and diagnostic as raw IEEE-754 bit patterns; for the
nine golden vectors also report the distance from the oracle's value.

Usage:  python bitident.py <old kh_server binary> <new kh_server binary>
Needs port 8005 free (kh_server binds it). Workspace harness used for finding
v0.4, 2026-09-26: 13 configurations, 127041 values, 0 differing.
"""
import glob, json, struct, subprocess, sys, time, urllib.request
from pathlib import Path

GOLDEN_DIR = Path(__file__).resolve().parents[3] / "shared" / "golden-vectors"

ORIG, FIXED = sys.argv[1], sys.argv[2]
configs = {}
GOLDEN = {}
for p in sorted(glob.glob(str(GOLDEN_DIR / "v*.json"))):
    gv = json.load(open(p)); prm = dict(gv["parameters"])
    prm["steps"] = prm.pop("nsteps")          # golden vectors say nsteps; the REST API says steps
    configs[Path(p).stem] = prm
    GOLDEN[Path(p).stem] = gv["diagnostics"]["kinetic_energy"]
n = 64 * 32
configs["forcing-mode"] = {"grid_nx": 64, "grid_ny": 32, "steps": 30, "forcing_type": "mode",
                           "forcing_amplitude": 0.5, "forcing_kx": 2, "forcing_ky": 1}
configs["forcing-band-rot-drag"] = {"grid_nx": 64, "grid_ny": 64, "domain_lx": 1.0, "domain_ly": 1.0,
                           "steps": 30, "forcing_type": "band", "forcing_amplitude": 0.3,
                           "forcing_kmin": 2, "forcing_kmax": 5, "forcing_omega0": 0.7, "linear_drag": 0.05}
configs["init_omega-legacy-key"] = {"grid_nx": 64, "grid_ny": 32, "steps": 20,
                           "init_omega": [0.01 * ((i * 7919) % 101 - 50) for i in range(n)]}
configs["nondefault-mix"] = {"grid_nx": 32, "grid_ny": 16, "steps": 17, "dt": 0.0008,
                           "reynolds_number": 800.0, "velocity_shear": 1.2, "perturbation_amplitude": 0.02,
                           "perturbation_mode": 3, "domain_lx": 1.5, "domain_ly": 0.75}

def run_all(binary):
    p = subprocess.Popen([binary], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); time.sleep(1)
    out = {}
    for name, body in configs.items():
        r = urllib.request.Request("http://127.0.0.1:8005/simulate", data=json.dumps(body).encode(), method="POST")
        out[name] = json.load(urllib.request.urlopen(r, timeout=300))
    p.kill(); p.wait()
    return out

def bits(x): return struct.pack("<d", float(x))

a, b = run_all(ORIG), run_all(FIXED)
total_vals = 0; bad = 0
for name in configs:
    ra, rb = a[name], b[name]
    diffs = []
    for k in ("kinetic_energy", "enstrophy", "max_vorticity", "divergence_rms", "cfl_peak"):
        total_vals += 1
        if bits(ra["diagnostics"][k]) != bits(rb["diagnostics"][k]): diffs.append(k)
    for f in ("u_velocity", "v_velocity", "vorticity", "pressure"):
        total_vals += len(ra[f])
        if len(ra[f]) != len(rb[f]) or any(bits(x) != bits(y) for x, y in zip(ra[f], rb[f])): diffs.append(f)
    bad += bool(diffs)
    grid = "%sx%s" % (rb.get("grid_nx"), rb.get("grid_ny"))
    ke = rb['diagnostics']['kinetic_energy']
    gold = f"vs golden {abs(ke-GOLDEN[name])/GOLDEN[name]:.1e}" if name in GOLDEN else ""
    print(f"{'IDENTICAL' if not diffs else 'DIFFERS  '} {name:26} {grid:8} steps {rb.get('steps_completed'):>4}  KE {ke:.15E}  {gold} {diffs or ''}")
print(f"\n{len(configs)} configurations, {total_vals} values compared bit-for-bit, {bad} differing")
sys.exit(1 if bad else 0)
