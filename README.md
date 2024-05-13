This file and the files listed in <MANIFEST.md> comprise a software bundle
designed to support the numerical experiments performed in [1]. Some (all) files
may appear in support bundles for other publications.

Main contributors:

- Cédric M. Campos
- Jose Torrente Teruel

Instructions:

1. Install [Julia](https://julialang.org/downloads/).
2. Install Plots.jl. In Julia's repl, run the following command:
   ```julia
   import Pkg; Pkg.add("Plots")
   ```
3. Run the test script. Assuming that the bundle's files are at the current
   path, in Julia's repl, run the following command:
   ```julia
   include("test_momdescso3.jl")
   ```
   The script will reproduce the last set experiment or, by default, the one
   corresponding to Fig. 5.2 in [1]. To set a particular experiment, set the
   variable `example` to a value among 11, 12, 21, 22, 31, or 32, like so:
   ```julia
   example = 21; include("test_momdescso3.jl")
   ```

References:

[1] C.M. Campos, D. Martín de Diego, J. Torrente  
    Momentum-based gradient descent methods for Lie groups

Copyright 2021-2024 Cédric M. Campos
