using DynamicPolynomials
using MultivariateMoments
using Distributions
using PDMats

const MM=MultivariateMoments

include("../src/BeurlingLasso.jl")
using .BeurlingLasso

params = [[-0.8], [0.7]]
weights = [0.3, 0.7]

domain = reshape([-1, 1], (1, 2))

natoms = length(params)
ndims = length(params[1])

@polyvar x[1:ndims]
μ0 = AtomicMeasure(x, [WeightedDiracMeasure(params[i], weights[i]) for i=1:natoms])

gmm = MixtureModel(Normal[Normal(param[1], 0.5) for param in params],
                   weights)

λ = 1e-3
σ = 0.5

freq_distrib = Normal(0, 2.0)
nfreqs = 20
freqs = rand(freq_distrib, nfreqs)

function φ(x, n)
   μ = x[1]
   ω = freqs[n]
   exp(im*ω*μ - σ^2*ω^2/2)
end

y = [BeurlingLasso.expectation(μ0, x -> φ(x, n)) for n=1:nfreqs]

# function θ(t, n)
#    ω = freqs[n]
#    exp(im*ω*t)
# end

# nsamples = 1000
# samples = rand(gmm, nsamples)
# y = [sum(θ(t, n) for t in samples) / nsamples for n=1:nfreqs]

gridsize = 1000
niter = 10

prob = blasso(y, φ, λ, domain)
solver = sfwsolver(niter, gridsize)

solve!(prob, solver)
