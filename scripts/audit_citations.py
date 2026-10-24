#!/usr/bin/env python3
"""audit_citations.py - doctrine C3 audit: textbook citations in the reference sources.

LEGE-ARTIS language doctrine s4.4 C3: every non-trivial routine of the hand-rolled
reference backends carries a comment block naming the textbook (author, year,
section, equation), its sign or convention, and its complexity. The authoritative
list is _audit/textbook-citations-allowlist.yaml; this script checks the sources
against it and fails on any mismatch.

Checks
  1. Every hand-rolled entry under `routines:` is defined in its file, and the
     comment block directly above the definition (blank lines between them allowed)
     contains: the first author's surname and the year of each citation; the
     section (written S12.2, sec. 12.2 or with the section sign); the equation and
     page when the allowlist gives them; the sign or convention text when given;
     the complexity claim. Library-backed entries are listed, not checked.
  2. Inventory: every public routine of the Fortran and Pascal sources (Fortran
     `public ::` lists; Pascal unit interface sections) appears in the allowlist -
     cited under `routines:`, or under `pending:` or `exempt:` with a reason. A new
     public routine without an entry fails the audit.
  3. Every `pending:` and `exempt:` symbol exists in its file (no stale entries).
  `pending:` entries (an anchor still to be chosen) and citations not yet verified
  against the printed book (`verified: false`) are reported but do not fail, except
  with --strict, the gate for a strict doctrine C3.

Usage (from the repository root, or with --root)
  python3 scripts/audit_citations.py               exit 0 = PASS
  python3 scripts/audit_citations.py --strict      pending and unverified fail too
  python3 scripts/audit_citations.py --self-test   proves the audit fails on defects
Needs PyYAML.
"""
import argparse
import os
import re
import shutil
import sys
import tempfile

import yaml

ALLOWLIST = '_audit/textbook-citations-allowlist.yaml'
LAST_ERRORS = []


def norm(s):
    return re.sub(r'\s+', ' ', s).strip()


def read_lines(path):
    with open(path, encoding='utf-8') as f:
        return f.read().splitlines()


# --- locating a definition ----------------------------------------------------------

def find_definition(lines, symbol, language):
    name = re.escape(symbol)
    if language == 'fortran':
        pat = re.compile(r'^\s*(?!end\b)(?:[\w(),=*\s]*?\s)?(?:subroutine|function)\s+' + name + r'\s*(\(|$)', re.I)
        start = 0
    elif language == 'pascal':
        pat = re.compile(r'^\s*(?:procedure|function)\s+' + name + r'\b', re.I)
        start = 0
        for i, l in enumerate(lines):
            if re.match(r'^\s*implementation\b', l, re.I):
                start = i + 1
                break
    elif language == 'cpp':
        pat = re.compile(r'^\s*(?:static\s+|inline\s+|constexpr\s+)*[\w:<>,&*\s]*?\b' + name + r'\s*\(')
        start = 0
    else:
        raise ValueError('unknown language ' + language)
    for i in range(start, len(lines)):
        if not pat.search(lines[i]):
            continue
        if language == 'cpp':
            tail = ' '.join(lines[i:i + 6])
            brace, semi = tail.find('{'), tail.find(';')
            if brace < 0 or (0 <= semi < brace):
                continue
        return i
    return None


def comment_block(lines, i, language):
    """Comment lines directly above line i (blank lines between block and line i allowed)."""
    j = i - 1
    while j >= 0 and not lines[j].strip():
        j -= 1
    out = []
    in_brace = False
    while j >= 0:
        s = lines[j].strip()
        if language == 'fortran':
            if not s.startswith('!'):
                break
            out.append(s)
        elif language == 'cpp':
            if not s.startswith('//'):
                break
            out.append(s)
        else:  # pascal: { } blocks and // lines
            if in_brace:
                out.append(s)
                k = s.find('{')
                if k >= 0 and not s.startswith('{$', k):
                    in_brace = False
            elif s.startswith('//'):
                out.append(s)
            elif s.endswith('}') and not s.startswith('{$'):
                out.append(s)
                k = s.find('{')
                in_brace = not (k >= 0 and not s.startswith('{$', k))
            else:
                break
        j -= 1
    return '\n'.join(reversed(out))


# --- inventory of public routines -----------------------------------------------------

def fortran_public(lines):
    text = '\n'.join(lines)
    text = re.sub(r'&\s*\n\s*&?', ' ', text)            # join continuation lines
    names = set()
    has_private = False
    for l in text.splitlines():
        code = l.split('!', 1)[0]
        if re.match(r'^\s*private\s*$', code, re.I):
            has_private = True
        m = re.match(r'^\s*public\s*::\s*(.+)$', code, re.I)
        if m:
            names.update(n.strip().lower() for n in m.group(1).split(',') if n.strip())
    defined = set(m.group(1).lower() for m in re.finditer(
        r'^\s*(?!end\b)(?:[\w(),=*\s]*?\s)?(?:subroutine|function)\s+(\w+)', text, re.I | re.M))
    if not has_private:
        return defined
    return names & defined


def pascal_public(lines):
    if not any(re.match(r'^\s*unit\b', l, re.I) for l in lines):
        return set()
    names, inside = set(), False
    for l in lines:
        if re.match(r'^\s*interface\b', l, re.I):
            inside = True
            continue
        if re.match(r'^\s*implementation\b', l, re.I):
            break
        if inside:
            m = re.match(r'^\s*(?:procedure|function)\s+(\w+)', l, re.I)
            if m:
                names.add(m.group(1).lower())
    return names


# --- the audit ------------------------------------------------------------------------

def surname(author):
    return author.split(',')[0].strip()


def check_entry(root, e, refs, errors):
    tag = '%s :: %s' % (e['file'], e['symbol'])
    path = os.path.join(root, e['file'])
    if not os.path.isfile(path):
        errors.append('%s: file not found' % tag)
        return False
    lines = read_lines(path)
    if e.get('status') != 'hand-rolled':          # library-backed: listed, not checked
        if not re.search(r'\b' + re.escape(e['symbol']) + r'\b', '\n'.join(lines)):
            errors.append('%s: symbol not found in the file (stale allowlist entry?)' % tag)
            return False
        return True
    i = find_definition(lines, e['symbol'], e['language'])
    if i is None:
        errors.append('%s: routine not defined in the file (stale allowlist entry?)' % tag)
        return False
    block = comment_block(lines, i, e['language'])
    if not block:
        errors.append('%s (line %d): no comment block above the definition - missing reference' % (tag, i + 1))
        return False
    flat = norm(block)
    missing = []
    for c in e.get('citations') or []:
        key = c.get('reference_key')
        ref = refs.get(key)
        if ref is None:
            missing.append('reference_key %s not in references:' % key)
            continue
        author = surname(ref['authors'][0])
        if author not in flat:
            missing.append('author "%s" (%s)' % (author, key))
        if str(ref['year']) not in flat:
            missing.append('year %s (%s)' % (ref['year'], key))
        sec = c.get('section')
        if sec:
            num = re.escape(str(sec).lstrip('\u00a7S').strip())
            if not re.search(r'(\u00a7\s*|\bS|\bsec\.?\s*|\bsection\s+|\bch\.?\s*|\bchapter\s+)' + num + r'(?![.\d])', flat, re.I):
                missing.append('section %s (%s)' % (sec, key))
        eq = c.get('equation')
        if eq and str(eq) not in flat:
            missing.append('equation %s (%s)' % (eq, key))
        page = c.get('page')
        if page and not re.search(r'\bp(p|age)?\.?\s*' + re.escape(str(page)) + r'\b', flat):
            missing.append('page %s (%s)' % (page, key))
    conv = e.get('sign_convention')
    if conv and norm(conv) not in flat:
        missing.append('convention "%s"' % conv)
    cx = e.get('complexity_claim')
    if cx and norm(cx) not in flat:
        missing.append('complexity "%s"' % cx)
    if missing:
        errors.append('%s (line %d): citation header lacks %s' % (tag, i + 1, '; '.join(missing)))
        return False
    return True


def run(root, strict=False, quiet=False):
    say = (lambda *a: None) if quiet else print
    with open(os.path.join(root, ALLOWLIST), encoding='utf-8') as f:
        al = yaml.safe_load(f)
    refs = al.get('references') or {}
    routines = al.get('routines') or []
    pending = al.get('pending') or []
    exempt = al.get('exempt') or []
    errors = []

    say('C3 citation audit - %s' % ALLOWLIST)
    n_ok = 0
    for e in routines:
        ok = check_entry(root, e, refs, errors)
        if ok:
            n_ok += 1
            say('  %-4s %-9s %s :: %s' % ('ok', e.get('status', ''), e['file'], e['symbol']))
        else:
            say('  FAIL %-9s %s :: %s' % (e.get('status', ''), e['file'], e['symbol']))

    listed = set()
    for e in routines:
        listed.add((e['file'], e['symbol'].lower()))
    for group, label in ((pending, 'pending'), (exempt, 'exempt')):
        for g in group:
            if not g.get('reason'):
                errors.append('%s %s: no reason given' % (label, g.get('file')))
            path = os.path.join(root, g['file'])
            if not os.path.isfile(path):
                errors.append('%s %s: file not found' % (label, g['file']))
                continue
            lines = read_lines(path)
            lang = g.get('language') or ('pascal' if path.endswith('.pas') else 'fortran' if path.endswith('.f90') else 'cpp')
            for s in g.get('symbols') or []:
                if find_definition(lines, s, lang) is None:
                    errors.append('%s %s :: %s: routine not defined in the file (stale entry)' % (label, g['file'], s))
                listed.add((g['file'], s.lower()))

    n_pub = 0
    for lang, sub, ext, fn in (('fortran', 'backends/fortran/src', '.f90', fortran_public),
                               ('pascal', 'backends/pascal/src', '.pas', pascal_public)):
        d = os.path.join(root, sub)
        if not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            if not name.endswith(ext):
                continue
            rel = sub + '/' + name
            for sym in sorted(fn(read_lines(os.path.join(d, name)))):
                n_pub += 1
                if (rel, sym) not in listed:
                    errors.append('%s :: %s: public routine not in the allowlist (cite it under routines:, '
                                  'or list it under pending: or exempt: with a reason)' % (rel, sym))

    n_pending = sum(len(g.get('symbols') or []) for g in pending)
    n_exempt = sum(len(g.get('symbols') or []) for g in exempt)
    n_hand = sum(1 for e in routines if e.get('status') == 'hand-rolled')
    say('  inventory: %d public Fortran/Pascal routines; allowlist: %d hand-rolled (checked), '
        '%d library-backed, %d pending, %d exempt'
        % (n_pub, n_hand, len(routines) - n_hand, n_pending, n_exempt))
    for g in pending:
        say('  pending  %s :: %s  (%s)' % (g['file'], ', '.join(g.get('symbols') or []), g['reason']))
    cits = [(e, c) for e in routines if e.get('status') == 'hand-rolled' for c in (e.get('citations') or [])]
    unverified = [(e, c) for e, c in cits if c.get('verified') is not True]
    say('  citations verified against the book: %d of %d' % (len(cits) - len(unverified), len(cits)))
    if strict and n_pending:
        errors.append('%d pending routine(s) - --strict requires every one cited' % n_pending)
    if strict and unverified:
        errors.append('%d citation(s) not yet verified against the book - --strict requires all' % len(unverified))
    LAST_ERRORS[:] = errors
    for err in errors:
        say('ERROR ' + err)
    if errors:
        say('RESULT: FAIL (%d error(s))' % len(errors))
        return 1
    say('RESULT: PASS (%d hand-rolled cited, %d pending, %d exempt)' % (n_hand, n_pending, n_exempt))
    return 0


# --- self-test: the audit must fail on each kind of defect ----------------------------

def self_test(root):
    with open(os.path.join(root, ALLOWLIST), encoding='utf-8') as f:
        al = yaml.safe_load(f)
    cited = [e for e in al['routines'] if e.get('status') == 'hand-rolled']
    first = cited[0]

    def mutate_drop_line(tmp, needle):
        p = os.path.join(tmp, first['file'])
        lines = read_lines(p)
        lines = [l for l in lines if needle not in l]
        open(p, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

    def mutate_rename(tmp):
        p = os.path.join(tmp, ALLOWLIST)
        d = yaml.safe_load(open(p, encoding='utf-8'))
        d['routines'][0]['symbol'] = d['routines'][0]['symbol'] + '_renamed'
        yaml.safe_dump(d, open(p, 'w', encoding='utf-8'), sort_keys=False, allow_unicode=True)

    def mutate_new_public(tmp):
        p = os.path.join(tmp, 'backends/fortran/src/kh_poisson.f90')
        s = open(p, encoding='utf-8').read()
        s = s.replace('contains', 'public :: kh_audit_probe\n\ncontains\n\n  subroutine kh_audit_probe()\n  end subroutine kh_audit_probe\n', 1)
        open(p, 'w', encoding='utf-8').write(s)

    cases = [
        ('complexity line removed from %s' % first['symbol'], lambda t: mutate_drop_line(t, 'Complexity:')),
        ('reference line removed from %s' % first['symbol'], lambda t: mutate_drop_line(t, 'Reference:')),
        ('allowlist symbol renamed (stale entry)', mutate_rename),
        ('new public routine without an allowlist entry', mutate_new_public),
    ]
    print('C3 audit self-test: the unmodified tree must pass, each mutation must fail')
    ok = run(root, quiet=True) == 0
    print('  %-4s unmodified tree passes' % ('ok' if ok else 'FAIL'))
    for label, fn in cases:
        tmp = tempfile.mkdtemp(prefix='c3-selftest-')
        try:
            for sub in ('_audit', 'backends'):
                shutil.copytree(os.path.join(root, sub), os.path.join(tmp, sub),
                                ignore=shutil.ignore_patterns('build', 'target', '*.o', '*.ppu', '*.mod'))
            fn(tmp)
            caught = run(tmp, quiet=True) != 0
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        print('  %-4s %s -> audit fails' % ('ok' if caught else 'FAIL', label))
        for err in LAST_ERRORS[:2]:
            print('         ' + err)
        ok = ok and caught
    print('RESULT: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--root', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    ap.add_argument('--strict', action='store_true')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    sys.exit(self_test(root) if a.self_test else run(root, strict=a.strict))


if __name__ == '__main__':
    main()
