using Distributions
import Statistics as stats
using DataFrames
using CSV
using Expectations
using Cumulants
using Plots 
gr()
using Polynomials
using LinearAlgebra, ForwardDiff
import PyPlot
import Contour: contours, levels, level, lines, coordinates


xs_ys(vs) = Tuple(eltype(vs[1])[vs[i][j] for i in 1:length(vs)] for j in eachindex(first(vs)))
xs_ys(v,vs...) = xs_ys([v, vs...])
xs_ys(r::Function, a, b, n=100) = xs_ys(r.(range(a, stop=b, length=n)))

function arrow!(p, v; kwargs...)
  if length(p) == 2
     quiver!(xs_ys([p])..., quiver=Tuple(xs_ys([v])); kwargs...)
  elseif length(p) == 3
    # 3d quiver needs support
    # https://github.com/JuliaPlots/Plots.jl/issues/319#issue-159652535
    # headless arrow instead
    plot!(xs_ys(p, p+v)...; kwargs...)
	end
end

using ForwardDiff
function D(f, n::Int=1)
    n < 0 && throw(ArgumentError("n is a non-negative integer"))
    n == 0 && return f
    n == 1 && return t -> ForwardDiff.derivative(f, float(t))
    D(D(f), n-1)
end
Base.adjoint(r::Function) = D(r)

r(t) = [sin(t), cos(t), t]
plot(xs_ys(r, 0, 4pi)...)
