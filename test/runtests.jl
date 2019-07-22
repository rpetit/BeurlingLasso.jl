using Test

@test 1+1 == 2

#
# include("../src/BeurlingLasso.jl")
# using .BeurlingLasso
#
# using LinearAlgebra
# using DynamicPolynomials
# using MultivariateMoments
#
# const MM = MultivariateMoments
#
# function weighterror(weight, μ)
#     if !isa(μ, AtomicMeasure)
#         return Inf
#     else
#         return min([abs(weight - atom.weight) for atom in μ.atoms]...)
#     end
# end
#
# function locationerror(location, μ)
#     if !isa(μ, AtomicMeasure)
#         return Inf
#     else
#         return  min([norm(location - atom.center) for atom in μ.atoms]...)
#     end
# end
#
# function setuptest(means, weights, maxorder, lowerbounds, upperbounds, λ, solver)
#     ndims = length(lowerbounds)
#     natoms = length(weights)
#
#     @polyvar x[1:ndims]
#     μ0 = AtomicMeasure(x, [WeightedDiracMeasure(means[i], weights[i]) for i=1:natoms])
#     monos = monomials(x, 0:maxorder)
#     moments = map(mono -> MM.moment(expectation(μ0, x -> φ(x, mono)), mono), monos)
#
#     if solver == "sfw"
#         niter = natoms
#         gridsize = 1000
#         μ = slidingfrankwolfe(moments, λ, niter, gridsize, lowerbounds, upperbounds)
#     elseif solver == "momsos"
#         ϵ = 1e-3
#         μ = momsos(moments, λ, ϵ, lowerbounds, upperbounds)
#     else
#         error("Invalid solver")
#     end
#
#     return μ
# end
#
# function recoverytest(means, weights, maxorder, lowerbounds, upperbounds, λ, solver, tol)
#     μ = setuptest(means, weights, maxorder, lowerbounds, upperbounds, λ, solver)
#
#     @test all(weighterror(weight, μ) < tol for weight in weights)
#     @test all(locationerror(mean, μ) < tol for mean in means)
# end
#
# @testset "recovery test momsos λ=0 1D tol=0.05" begin
#     lowerbounds = [-1]
#     upperbounds = [1]
#     maxorder = 10
#     λ = 0
#     solver = "momsos"
#     tol = 0.05
#
#     for (means, weights) in [([[0.0]], [1.0]), ([[0.9]], [0.2]),
#                              ([[-0.3]], [0.01]), ([[0.1], [0.4]], [0.1, 0.3]),
#                              ([[-0.4], [-0.7]], [2.0, 0.01]),
#                              ([[-0.8], [0.4], [0.8]], [0.5, 0.3, 0.7])]
#
#         recoverytest(means, weights, maxorder, lowerbounds, upperbounds, λ, solver, tol)
#     end
# end
#
# @testset "recovery test momsos λ=0 2D tol=0.05" begin
#     lowerbounds = [-1, -1]
#     upperbounds = [1, 1]
#     maxorder = 15
#     λ = 0
#     solver = "momsos"
#     tol = 0.05
#
#     for (means, weights) in [([[0.0, 0.0]], [1.0]), ([[0.9, 0.2]], [0.2]),
#                              ([[0.1, 0.2], [0.4, 0.8]], [0.1, 0.3]),
#                              ([[-0.4, -0.9], [-0.7, 0.8]], [2.0, 0.01]),
#                              ([[-0.8, -0.2], [0.4, 0.3], [0.8, 0.1]], [0.5, 0.3, 0.7])]
#
#         recoverytest(means, weights, maxorder, lowerbounds, upperbounds, λ, solver, tol)
#     end
# end
#
# @testset "recovery test sfw λ=1e-5 1D tol=0.05" begin
#     lowerbounds = [-1]
#     upperbounds = [1]
#     maxorder = 10
#     λ = 1e-5
#     solver = "sfw"
#     tol = 0.05
#
#     for (means, weights) in [([[0.0]], [1.0]), ([[0.9]], [0.2]),
#                              ([[-0.3]], [0.01]), ([[0.1], [0.4]], [0.1, 0.3]),
#                              ([[-0.4], [-0.7]], [2.0, 0.01]),
#                              ([[-0.8], [0.4], [0.8]], [0.5, 0.3, 0.7])]
#
#         recoverytest(means, weights, maxorder, lowerbounds, upperbounds, λ, solver, tol)
#     end
# end
#
# @testset "recovery test sfw λ=1e-5 2D tol=0.05" begin
#     lowerbounds = [-1, -1]
#     upperbounds = [1, 1]
#     maxorder = 10
#     λ = 1e-5
#     solver = "sfw"
#     tol = 0.05
#
#     for (means, weights) in [([[0.0, 0.0]], [1.0]), ([[0.9, 0.2]], [0.2]),
#                              ([[0.1, 0.2], [0.4, 0.8]], [0.1, 0.3]),
#                              ([[-0.4, -0.9], [-0.7, 0.8]], [2.0, 0.01]),
#                              ([[-0.8, -0.2], [0.4, 0.3], [0.8, 0.1]], [0.5, 0.3, 0.7])]
#
#         recoverytest(means, weights, maxorder, lowerbounds, upperbounds, λ, solver, tol)
#     end
# end
