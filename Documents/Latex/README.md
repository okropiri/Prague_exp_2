# NCAL LaTeX Build

Build from the repository root with:

```bash
Documents/Latex/compile.sh
```

The helper prefers a full TeX Live setup with `latexmk` and `lualatex` when both are available. If they are not installed, it falls back to `tectonic`, which is sufficient for the current working document on this Mac.

For this workspace, `tectonic` is installed in:

```bash
~/.local/bin/tectonic
```

The output PDF is written to:

```bash
Documents/Latex/main.pdf
```

The theme keeps LuaLaTeX-specific Urbanist variable-font weight settings when compiled with LuaLaTeX. The Tectonic/XeTeX fallback uses the same bundled Urbanist font files with portable filenames.