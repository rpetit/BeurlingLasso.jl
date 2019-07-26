function momsos!(blasso::BLASSO{T, FT, DT},
                 ϵ::Real,
                 maxdeg::Union{Int, String}) where {T,
                                                    PT<:AbstractPolynomialLike,
                                                    FT<:Array{PT, 1},
                                                    DT<:AbstractSemialgebraicSet}

    model = SOSModel(with_optimizer(Mosek.Optimizer, QUIET=true))

    y = blasso.y
    φ = blasso.φ
    m = length(y)

    @variable(model, p[1:m])  # dual variable

    if maxdeg == "defaut"
        maxdeg = max([maxdegree(monomials(φ[i])) for i=1:m]...)
    end

    @constraint(model, cminus,
                1 + sum(φ[i] * p[i] for i=1:m) in SOSCone(),
                domain=blasso.domain, maxdegree=maxdeg)

    @constraint(model, cplus,
                1 - sum(φ[i] * p[i] for i=1:m) in SOSCone(),
                domain=blasso.domain, maxdegree=maxdeg)

    if blasso.λ == 0  # noiseless formulation
        @objective(model, Max, sum(y[i] * p[i] for i=1:m))
    else  # regularized formulation
        # workaround to avoid MosekTools error with quadratic objective
        @variable(model, t)
        @constraint(model, sum(p[i]^2 for i=1:m) <= t)

        @objective(model, Max, sum(y[i] * p[i] for i=1:m) - blasso.λ / 2 * t)
    end

    optimize!(model)
    blasso.objvalue = objective_value(model)
    blasso.dualsol = [value(p[i]) for i=1:m]

    # finds a primal solution given the dual solution
    νplus = moment_matrix(cplus)
    μplus = extractatoms(νplus, ϵ)

    νminus = moment_matrix(cminus)
    μminus = extractatoms(νminus, ϵ)

    # this will raise an error if extraction failed
    # (i.e. μplus or μminus is of type Nothing)
    atoms = μplus.atoms
    append!(atoms, [WeightedDiracMeasure(atom.center, -atom.weight)
                    for atom in μminus.atoms])

    blasso.μ = AtomicMeasure(μplus.variables, atoms)
end
