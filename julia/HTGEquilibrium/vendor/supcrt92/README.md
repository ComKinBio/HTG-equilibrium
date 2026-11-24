# SUPCRT92 water code (third-party)

`H2O92D.f` is the water routine of SUPCRT92 (J. W. Johnson, Lawrence Livermore
National Laboratory; file timestamp 15 November 1991): the Haar–Gallagher–Kell /
Levelt Sengers equation of state for water and the Johnson–Norton dielectric
constant. HTGEquilibrium compiles it (`deps/build.jl`) and calls it from Julia.

- **Copied from:** the Python package pychnosz 1.5.4, file
  `pychnosz/data/extdata/src/H2O92D.f.orig`, as installed in this repository's
  `.venv`. pychnosz describes it as the original SUPCRT92 file without
  modification.
- **Copied on:** 24 November 2026.
- **Changes:** none. The file is byte for byte identical to the source
  (sha256 `885cf390d235024db0e96c249c3e6d533cd08c31be672b26c50b90786cbc3e46`,
  checked by the package tests).
- **Licence:** the file has no licence header. The owner accepted using it as is
  (PR #1).
