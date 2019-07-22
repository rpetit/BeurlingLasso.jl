export expectation, computehermite, unravel_index, rieszmap, extractmoments

# expectation of a function with respect to an atomic measure
function expectation(μ::AbstractMeasureLike, f::Function)
    sum(δ -> δ.weight * f(δ.center), μ.atoms)
end

# computation of Hermite polynomials' coefficients
function computehermite(maxdegree::Int64)
    a = zeros(Int64, maxdegree+1, maxdegree+1)
    a[1, 1] = 1
    a[2, 2] = 2

    for n=2:maxdegree
        a[n+1, 1] = - a[n, 2]
        a[n+1, n] = 2 * a[n, n-1]
        a[n+1, n+1] = 2 * a[n, n]
        for k=2:maxdegree
            a[n+1, k] = 2 * a[n, k-1] - k * a[n, k+1]
        end
    end

    return a
end

function unravel_index(x, n)
    return (div(x - 1, n) + 1, mod(x - 1, n) + 1)
end

# it is implicitly assumed that a given monomial cannot be found twice in y
function rieszmap(y::Array{MT, 1}, φ::AbstractPolynomialLike) where {MT<:AbstractMoment}
    sum(moment_value(y[i]) * coefficient(φ, monomial(y[i])) for i=1:length(y))
end

function extractmoments(M::MomentMatrix)
    n = M.Q.n
    σ, monos = sortmonovec(vcat([[M.x[i] * M.x[j] for j=1:n] for i=1:n]...))
    [moment(M.Q[unravel_index(σ[i], n)...], monos[i]) for i=1:length(σ)]
end
