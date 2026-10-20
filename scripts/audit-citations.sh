#!/bin/sh
# audit-citations.sh - doctrine C3 citation audit (see scripts/audit_citations.py).
# Usage: sh scripts/audit-citations.sh [--strict | --self-test]      exit 0 = PASS
exec python3 "$(dirname "$0")/audit_citations.py" "$@"
