#!/usr/bin/env bash
# Download, checksum-verify and compile plantri (third-party; NOT redistributed here).
#   plantri by Gunnar Brinkmann and Brendan McKay (with contributions by others),
#   https://users.cecs.anu.edu.au/~bdm/plantri/ , Apache License 2.0.
# Usage: scripts/get_plantri.sh [55|58]      (default 58)
#   58 = plantri 5.8 (March 4, 2026), used by the independent re-verification
#   55 = plantri 5.5 (May 17, 2024),  used by the primary computation
# Result: scripts/_build/plantri<ver>/plantri   (the _build/ directory is not part of the repo)
set -euo pipefail
VER="${1:-58}"
case "$VER" in
  58) SHA=e78a944116fec9f2c9f5e484206276cc2b0043bae803e9815f4b2683614629b8 ;;
  55) SHA=911cdf5bcca7294eb80f8f79fefc148183f7ba81da15b3aa4d6d2401a3bc7ded ;;
  *) echo "unknown version $VER" >&2; exit 2 ;;
esac
HERE="$(cd "$(dirname "$0")" && pwd)"
B="$HERE/_build"; mkdir -p "$B"; cd "$B"
TGZ="plantri${VER}.tar.gz"
[ -f "$TGZ" ] || curl -sSL -o "$TGZ" "https://users.cecs.anu.edu.au/~bdm/plantri/$TGZ"
echo "$SHA  $TGZ" | sha256sum -c -
rm -rf "plantri${VER}"; tar xzf "$TGZ"
cd "plantri${VER}" && ${CC:-cc} -o plantri -O3 plantri.c
./plantri 2>&1 | head -0 || true
echo "built $B/plantri${VER}/plantri"
