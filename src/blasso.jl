export blasso, sfwsolver, momsossolver, solve!

mutable struct BLASSO{T, FT, DT}
    y::Array{T, 1}
    φ::FT
    λ::Real
    domain::DT

    μ::Union{AbstractMeasureLike, Nothing}
    objvalue::Union{Real, Nothing}
    dualsol::Union{Vector{Real}, Nothing}
    mommat::Union{MomentMatrix, Nothing}
end

function blasso(y::Array{T, 1}, φ::Array{PT, 1}, λ::Real,
                domain::AbstractSemialgebraicSet) where {T, PT<:AbstractPolynomialLike}
    BLASSO(y, φ, λ, domain, nothing, nothing, nothing, nothing)
end

function blasso(y::Array{T, 1}, φ::Function, λ::Real,
                domain::Array{BT, 2}) where {T, BT}
    BLASSO(y, φ, λ, domain, nothing, nothing, nothing, nothing)
end

struct Solver
    name::String
    options::Dict
end

function sfwsolver(niter::Integer, gridsize::Integer)
    Solver("sfw", Dict("niter" => niter, "gridsize" => gridsize))
end

function momsossolver(ϵ::Real, maxdeg::Union{Int, String}="defaut")
    Solver("momsos", Dict("ϵ" => ϵ, "maxdeg" => maxdeg))
end

function solve!(blasso::BLASSO, solver::Solver)
    if solver.name == "momsos"
        momsos!(blasso, solver.options["ϵ"], solver.options["maxdeg"])
    elseif solver.name == "sfw"
        sfw!(blasso, solver.options["niter"], solver.options["gridsize"])
    else
        throw(ArgumentError(solver.name,
                            "solver name must be either 'momsos' or 'sfw'"))
    end
end
