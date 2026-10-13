"""KB-039: transliterate non-ASCII characters in Pascal sources to ASCII.
Comments only - the script refuses to touch a file if any non-ASCII character
occurs outside a { } or (* *) comment or a // line comment."""
import sys, re

MAP = {"─": "-", "—": "-", "·": "*", "̂": "_hat", "ψ": "psi",
       "ω": "omega", "→": "->", "∂": "d", "²": "^2", "π": "pi",
       "ν": "nu", "∇": "nabla", "û": "u_hat", "§": "S"}

def comment_spans(s):
    spans, i, n = [], 0, len(s)
    while i < n:
        if s[i] == "'":                                   # string literal
            j = s.index("'", i + 1)
            while j + 1 < n and s[j + 1] == "'": j = s.index("'", j + 2)
            i = j + 1
        elif s[i] == "{":
            j = s.index("}", i); spans.append((i, j)); i = j + 1
        elif s.startswith("(*", i):
            j = s.index("*)", i); spans.append((i, j + 1)); i = j + 2
        elif s.startswith("//", i):
            j = s.find("\n", i); j = n if j < 0 else j; spans.append((i, j)); i = j
        else:
            i += 1
    return spans

for path in sys.argv[1:]:
    s = open(path, encoding="utf-8").read()
    spans = comment_spans(s)
    inside = lambda k: any(a <= k <= b for a, b in spans)
    outside = [k for k, ch in enumerate(s) if ord(ch) > 127 and not inside(k)]
    if outside:
        sys.exit(f"{path}: non-ASCII outside comments at offsets {outside[:5]} - refusing")
    unknown = sorted({ch for ch in s if ord(ch) > 127 and ch not in MAP})
    if unknown:
        sys.exit(f"{path}: no transliteration for {unknown} - refusing")
    # a combining circumflex follows its base letter: "ψ̂" -> "psi_hat"
    out = "".join(MAP.get(ch, ch) for ch in s)
    open(path, "w", encoding="ascii", newline="").write(out)
    print(f"{path}: {sum(1 for ch in s if ord(ch) > 127)} characters transliterated")
