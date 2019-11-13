# BeurlingLasso.jl

[![Build Status](https://travis-ci.com/rpetit/BeurlingLasso.jl.svg?branch=master)](https://travis-ci.com/rpetit/BeurlingLasso.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/rpetit/BeurlingLasso.jl?svg=true)](https://ci.appveyor.com/project/rpetit/BeurlingLasso-jl)
[![Codecov](https://codecov.io/gh/rpetit/BeurlingLasso.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rpetit/BeurlingLasso.jl)
[![Coveralls](https://coveralls.io/repos/github/rpetit/BeurlingLasso.jl/badge.svg?branch=master)](https://coveralls.io/github/rpetit/BeurlingLasso.jl?branch=master)

Numerical solvers for the Beurling LASSO
<br/>

## Available algorithms
* Sliding Frank-Wolfe
* Lasserre hierarchies

## Disclaimer

This package is under active development. Tests and documentation still
have to be written.
<br/>

## References

[1] Q. Denoyelle, V. Duval, G. Peyré and E. Soubies, the sliding Frank-Wolfe
algorithm and its application to super-resolution microscopy, Inverse Problems,
2019 <br/>

[2] Y. De Castro, F. Gamboa, D. Henrion, J.-B. Lasserre, exact solutions to super-resolution on semi-algebraic domains in higher dimensions, IEEE Transactions on Information Theory, 2017 <br/>

See also these repositories :

* [qdenoyelle/sfw4blasso](https://github.com/qdenoyelle/sfw4blasso)
* [nboyd/SparseInverseProblems.jl](https://github.com/nboyd/SparseInverseProblems.jl)
* [JuliaOpt/SumOfSquares.jl](https://github.com/JuliaOpt/SumOfSquares.jl)
