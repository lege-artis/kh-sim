"""All hostile cases on ONE long-lived process, in sequence. After each case:
/health must answer, and the reference simulation must return a kinetic energy
bit-identical to the one taken before any hostile input (no poisoning)."""
import json, subprocess, sys, time, urllib.request, urllib.error

BIN = sys.argv[1]
REF = {"grid_nx": 64, "grid_ny": 32, "steps": 20}
BASE = {"grid_nx": 64, "grid_ny": 32, "steps": 5}

def post_raw(data, timeout=20):
    r = urllib.request.Request("http://127.0.0.1:8005/simulate", data=data, method="POST")
    try:
        with urllib.request.urlopen(r, timeout=timeout) as resp:
            return resp.status, resp.read().decode(errors="replace")
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode(errors="replace")
    except Exception as e:
        return 0, type(e).__name__

def health():
    try:
        with urllib.request.urlopen("http://127.0.0.1:8005/health", timeout=3) as r: return r.status
    except Exception: return 0

def ref_ke():
    c, t = post_raw(json.dumps(REF).encode())
    return json.loads(t)["diagnostics"]["kinetic_energy"] if c == 200 else None

cases = []
for k, vals in [("grid_nx", [0, -64, 3, 48, 96, 100, 65537, 2**30]), ("grid_ny", [0, -32, 24, 48]),
                ("domain_lx", [0.0]), ("domain_ly", [0.0]), ("reynolds_number", [0.0]),
                ("dt", [1.0]), ("velocity_shear", [0.0]), ("perturbation_amplitude", [0.0]),
                ("forcing_type", ["bogus"])]:
    for v in vals:
        b = dict(BASE); b[k] = v; cases.append((f"{k}={v}", json.dumps(b).encode()))
cases += [
    ("grid 4096x4096 (cell cap)", json.dumps({"grid_nx": 4096, "grid_ny": 4096, "steps": 1}).encode()),
    ("band forcing lx!=ly", json.dumps(dict(BASE, forcing_type="band", forcing_amplitude=1.0, forcing_kmin=1, forcing_kmax=3)).encode()),
    ("init_omega wrong length", json.dumps(dict(BASE, init_omega=[0.0] * 10)).encode()),
    ("init_omega right length", json.dumps(dict(BASE, init_omega=[0.0] * 2048)).encode()),
    ("body: not JSON", b"this is not json"),
    ("body: JSON array", b"[1,2,3]"),
    ("body: empty", b""),
    ("grid_nx as string", json.dumps(dict(BASE, grid_nx="sixty-four")).encode()),
]

p = subprocess.Popen([BIN], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
for _ in range(50):
    if health() == 200: break
    time.sleep(0.05)
ke0 = ref_ke()
print(f"reference KE before hostile inputs: {ke0!r}\n")
bad = 0
for name, data in cases:
    code, txt = post_raw(data)
    try:
        j = json.loads(txt); info = j.get("error") and (f"{j.get('field') or j.get('backend')}: {j['error']}") or "200 ok"
    except Exception:
        info = txt[:50]
    h = health(); ke = ref_ke()
    ok = (h == 200) and (ke == ke0) and code in (200, 400, 500)
    bad += not ok
    flag = "ok  " if ok else "FAIL"
    print(f"{flag} {name:28} HTTP {code:3}  health {h:3}  ref-identical {str(ke == ke0):5}  {info[:70]}")
p.kill()
print(f"\n{len(cases)} cases, {bad} failures")
sys.exit(1 if bad else 0)
