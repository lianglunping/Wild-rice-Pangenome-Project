#!/usr/bin/env bash
set -euo pipefail

PREFIX="${1:-$HOME/.local/opt/DupGen_finder}"
REPOSITORY="https://github.com/qiao-xin/DupGen_finder.git"

if [[ -e "$PREFIX" ]]; then
  echo "Refusing to overwrite existing path: $PREFIX" >&2
  exit 2
fi

git clone "$REPOSITORY" "$PREFIX"
make -C "$PREFIX"

cat <<MSG
Installed DupGen_finder under:
  $PREFIX

Add the directory to PATH:
  export PATH="$PREFIX:\$PATH"

Record the source commit:
  git -C "$PREFIX" rev-parse HEAD
MSG
