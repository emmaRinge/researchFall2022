using Distributions
import Statistics as stats
using DataFrames
using CSV
using LinearAlgebra
using Expectations
using Cumulants
using Plots
using Polynomials

#generate data
mu1 = rand(0.1:1)
sigma1 = rand(0.1:5)
mu2 = rand(0.1:2)
sigma2 = rand(0.1:2)
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
x = rand(d, 10)
y = rand(e, 10)
z = x .* y

xlin = range(findmin(z)[1], findmax(z)[1], length = 10)
ys = convert(AbstractVector, @.sin(z))
f = Polynomials.fit(z, ys, 3)
f2 = Polynomials.fit(z, ys, 6)

scatter(z, ys, markersize=3, label="data")
plot!(xlin, f.(xlin), label = "f(x) = x3", dpi = 200)
plot!(xlin, f2.(xlin), label = "f(x) = x3", dpi = 200)
