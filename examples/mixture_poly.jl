using DynamicPolynomials
using MultivariateMoments
using SemialgebraicSets

const MM = MultivariateMoments

include("../src/BeurlingLasso.jl")
using .BeurlingLasso

@polyvar x[1:2]

# only valid for odd numbers
function doublefactorial(n::Integer)
    if n <= 0
        1
    else
        prod(2*k+1 for k=0:div(n-1, 2))
    end
end

nmoments = 15

φ = [sum(binomial(n, 2*k) * x[1]^(n-2*k) * x[2]^(2*k) * doublefactorial(2*k-1)
     for k=0:div(n, 2)) for n=0:nmoments-1]

params = [[-0.5, 0.5], [0.6, 0.3]]
weights = [0.7, 0.3]

natoms = length(params)

μ0 = AtomicMeasure(x, [WeightedDiracMeasure(params[i], weights[i]) for i=1:natoms])
y = [MM.expectation(μ0, φ[n]) for n=1:nmoments]

# λ = 0
#
# domain = @set x[1] >= -1 && x[1] <= 1 && x[2] >= -1 && x[2] <= 1
#
# prob1 = blasso(y, φ, λ, domain)
# solver1 = momsossolver(1e-3)
# solve!(prob1, solver1)

λ = 1e-5

φ_func(u, n) = φ[n](x=>u)
bounds = [[-1 1]; [-1 1]]

prob2 = blasso(y, φ_func, λ, bounds)
solver2 = sfwsolver(2, 1000)
solve!(prob2, solver2)
