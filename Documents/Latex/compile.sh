#!/usr/bin/env sh
set -eu

cd "$(dirname "$0")"

if command -v latexmk >/dev/null 2>&1 && command -v lualatex >/dev/null 2>&1; then
  exec latexmk -lualatex -interaction=nonstopmode main.tex
fi

if command -v tectonic >/dev/null 2>&1; then
  exec tectonic main.tex
fi

if [ -x "$HOME/.local/bin/tectonic" ]; then
  exec "$HOME/.local/bin/tectonic" main.tex
fi

echo "No LaTeX compiler found. Install latexmk+lualatex or Tectonic." >&2
exit 127