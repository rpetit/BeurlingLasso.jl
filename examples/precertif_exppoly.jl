using BlockArrays
using LinearAlgebra
using JuMP
using Ipopt
using Plots
using PyPlot
using Plots.PlotMeasures
using Formatting
using Printf

pyplot()
PyPlot.matplotlib.rc("text", usetex=true)
PyPlot.matplotlib.rc("font", family="serif")
Plots.scalefontsizes(1.5)

φ(x, n, σ) = exp(-x^2 / (2*σ^2)) * x^n

function dφ(x, n, σ)
    if n == 0
        exp(-x^2 / (2*σ^2)) * (-x/σ^2)
    else
        exp(-x^2 / (2*σ^2)) * (n * x^(n-1) - x^(n+1)/σ^2)
    end
end

K(x, y, σ, N) = sum(φ(x, n, σ) * φ(y, n, σ) for n=0:N)
∂K(x, y, σ, N) = sum(φ(x, n, σ) * dφ(y, n, σ) for n=0:N)
∂∂K(x, y, σ, N) = sum(dφ(x, n, σ) * dφ(y, n, σ) for n=0:N)

function D(x, σ, N)
    s = length(x)
    D0 = [K(x[i], x[j], σ, N) for i=1:s, j=1:s]
    D1 = [∂K(x[i], x[j], σ, N) for i=1:s, j=1:s]
    D2 = [∂∂K(x[i], x[j], σ, N) for i=1:s, j=1:s]
    [D0 D1; D1 D2]
end

function compute_ηV(x0, σ, N)
    s = length(x0)

    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(model, p[1:N+1])

    for i=1:s
        @constraint(model, sum(p[n+1] * φ(x0[i], n, σ) for n=0:N) == 1)
        @constraint(model, sum(p[n+1] * dφ(x0[i], n, σ) for n=0:N) == 0)
    end

    @objective(model, Min, sum(p[n]^2 for n=1:N+1))
    optimize!(model)

    pV = [value(p[n]) for n=1:N+1]
    ηV(x) = sum(φ(x, n, σ) * pV[n+1] for n=0:N)
end

fetitle = FormatExpr(L"$\sigma={}$")

function plot_ηV(x0, ηV, σ)
    s = length(x0)
    numpoints = 500

    x_tab = range(-1, stop=1, length=numpoints)
    y_tab = [ηV(x) for x in x_tab]

    Plots.plot(x_tab, y_tab, gridcolor="lightgrey", title=format(fetitle, σ))
    hline!([1], linestyle=:dash, linecolor=:black)

    if min(y_tab...) < -0.05
        hline!([-1], linestyle=:dash, linecolor=:black)
    end

    plot!(x0, ones(s), line=:stem, linestyle=:dash, marker=:circle, markersize=5,
          color="red")
    plot!(xlabel=L"$x$", legend=false, yformatter = yi -> @sprintf("%0.2f", yi))
    plot!(bottom_margin=5mm, right_margin=10mm, left_margin=5mm, top_margin=5mm)
end

N = 9

l = @layout[a; [b c d]; [e f g]]

p0 = Plots.plot((1:1)', label=L"$\eta_V$")
plot!(p0, (1:1)', label=L"$x_0$", color="red", linestyle=:dash,
      legend=:bottom, framestyle=:none)

x0 = [-0.8, 0.1, 0.3, 0.9]

σ = 0.285
ηV = compute_ηV(x0, σ, N)
p1 = plot_ηV(x0, ηV, σ)

σ = 0.315
ηV = compute_ηV(x0, σ, N)
p2 = plot_ηV(x0, ηV, σ)

σ = 0.375
ηV = compute_ηV(x0, σ, N)
p3 = plot_ηV(x0, ηV, σ)

x0 = [-0.1, 0.0, 0.05, 0.1]

σ = 0.05
ηV = compute_ηV(x0, σ, N)
p4 = plot_ηV(x0, ηV, σ)

σ = 0.3
ηV = compute_ηV(x0, σ, N)
p5 = plot_ηV(x0, ηV, σ)

σ = 0.5
ηV = compute_ηV(x0, σ, N)
p6 = plot_ηV(x0, ηV, σ)

Plots.plot(p0, p1, p2, p3, p4, p5, p6,
           layout=l, link=:x, size=(1000, 750),
           left_margin=5mm, right_margin=10mm)
