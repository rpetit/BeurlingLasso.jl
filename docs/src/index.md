# Introduction

BeurlingLasso.jl is a Julia package for solving the following optimization problem over measures, which is known as the Beurling LASSO:

```math
\begin{equation*}
  \underset{ \mu\in\mathcal{M}(\mathcal{X})}{\text{inf}} ~~
  \frac{1}{2}||\Phi\mu-y||_{\mathcal{H}}^2 + \lambda |\mu|(\mathcal{X})
  \label{eqn:noisy_blasso}
\end{equation*}
```

Two solvers are available: one relies on semi-definite relaxations (Lasserre hierarchies), while the other implements the sliding Frank-Wolfe algorithm.
