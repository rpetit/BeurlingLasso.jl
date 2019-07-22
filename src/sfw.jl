export sfw!

# construction of a regular grid with given bounds and number of points
function makegrid(lowerbounds, upperbounds, gridsize::Int64)
    ndims = length(lowerbounds)
    grid = collect(Iterators.product(
        [[lowerbounds[n] + (upperbounds[n] - lowerbounds[n]) * (i-1) / (gridsize-1)
          for i=1:gridsize] for n=1:ndims]...))
end

# finds the maximum of a function on a grid
function gridsearch(η, lowerbounds, upperbounds, gridsize)
    grid = makegrid(lowerbounds, upperbounds, gridsize)
    argmax_temp = grid[1]

    for x in grid
        if abs(η(x)) > abs(η(argmax_temp))
            argmax_temp = x
        end
    end

    collect(argmax_temp)  # tuple to array conversion
end

# attempts to find the maximum of a function using Newton algorithm, initialized
# with the output of a grid search (no guarantee to reach a global optimum)
function findargmax(η, lowerbounds, upperbounds, gridsize)
    x0 = gridsearch(η, lowerbounds, upperbounds, gridsize)
    obj(x) = -abs(η(x))

    # TODO: check autodiff vs analytic
    Optim.minimizer(optimize(obj, x0, Newton(); autodiff=:forward))
end

# adds a spike with 0 amplitude at a given location
function addspike(μ, new_x)
    new_atom = WeightedDiracMeasure(new_x, 0.0)
    AtomicMeasure(μ.variables, append!(μ.atoms, [new_atom]))
end

# finds optimal amplitudes with fixed locations by solving a LASSO problem
function updateamplitudes(μ, y, φ, λ)
    natoms = length(μ.atoms)
    nmoments = length(y)

    x = [atom.center for atom in μ.atoms]
    Φx = [φ(x[i], n) for n=1:nmoments, i=1:natoms]

    # constraints = hcat([[0 Inf]' for i=1:natoms]...)  # positive spikes
    # lasso = glmnet(Φx, y, lambda=[λ], constraints=constraints)
    lasso = glmnet(vcat(real(Φx), imag(Φx)), vcat(real(y), imag(y)), lambda=[λ])
    new_a = lasso.betas[:]

    new_atoms = [WeightedDiracMeasure(μ.atoms[i].center, new_a[i]) for i=1:natoms]
    AtomicMeasure(μ.variables, new_atoms)
end

# attempts to find optimal amplitudes and locations by solving a non convex
# problem with BFGS algorithm (no guarantee to reach a global optimum)
function performsliding(μ, y, φ, λ)
    natoms = length(μ.atoms)
    ndims = length(μ.variables)
    nmoments = length(y)

    function obj(arg)
        # amplitudes and locations are parsed from the arg vector
        a = arg[1:natoms]
        # x = reshape(arg[natoms+1:length(arg)], (ndims, natoms))'
        x = [[arg[natoms + (i-1) * ndims + j] for j=1:ndims] for i=1:natoms]

        nspikes = length(μ.atoms)
        Φx = [sum(a[i] * φ(x[i], n) for i=1:natoms)
              for n=1:nmoments]

        sum(abs(Φx[n] - y[n])^2 for n=1:nmoments) / 2 + λ * sum(abs(a[i]) for i=1:nspikes)
    end

    # amplitudes and locations are concatenated in a single vector to cope
    # with Optim.jl API
    cat_centers = vcat([atom.center for atom in μ.atoms]...)
    x0 = vcat([atom.weight for atom in μ.atoms]..., cat_centers)

    # TODO: check autodiff vs analytic
    res = Optim.minimizer(optimize(obj, x0, BFGS(); autodiff=:forward))

    new_a = res[1:natoms]
    new_x = [[res[natoms + (i-1) * ndims + j] for j=1:ndims] for i=1:natoms]
    new_atoms = [WeightedDiracMeasure(new_x[i], new_a[i]) for i=1:natoms]

    AtomicMeasure(μ.variables, new_atoms)
end

# Sliding Frank-Wolfe algorithm
function sfw!(blasso::BLASSO{T, FT, DT},
              niter::Integer,
              gridsize::Integer) where {T, FT<:Function, BT, DT<:Array{BT, 2}}

    lowerbounds = blasso.domain[:, 1]
    upperbounds = blasso.domain[:, 2]

    φ = blasso.φ
    y = blasso.y

    ndims = length(lowerbounds)
    nmoments = length(y)

    @polyvar var[1:ndims]

    for k=1:niter
        if k == 1
            Φμ = zeros(nmoments)  # workaround to simulate a null atomic measure
        else
            Φμ = [expectation(μ, x -> φ(x, n)) for n=1:nmoments]
        end

        η(x) = sum((y[n] - Φμ[n]) * φ(x, n) for n=1:length(y)) / blasso.λ
        new_x = findargmax(η, lowerbounds, upperbounds, gridsize)

        if k == 1
            global μ
            μ = AtomicMeasure(var, [WeightedDiracMeasure(new_x, 0.0)])
        else
            μ = addspike(μ, new_x)
        end

        if abs(η(new_x)) <= 1  # check convergence
            @goto convergence_label
        end

        μ = updateamplitudes(μ, y, φ, blasso.λ)
        μ = performsliding(μ, y, φ, blasso.λ)
    end

    @label convergence_label
    blasso.μ = μ
end
