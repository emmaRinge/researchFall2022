using Distributions
import Statistics as stats
using DataFrames
using CSV
using LinearAlgebra
using Expectations
using Cumulants
using Plots 
gr()
using Polynomials

#generate data
mu1 = rand(0:1)
sigma1 = rand(0.1:5)
mu2 = rand(0:1)
sigma2 = rand(0.1:5)
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
x = rand(d, 10)
y = rand(e, 10)
z = x .* y

xlin = range(findmin(z)[1], findmax(z)[1], length = 10)
ys = convert(AbstractVector, @.sin(z))
f = Polynomials.fit(z, ys, 3)
f2 = Polynomials.fit(z, ys, 6)

x = range(0, 10, length = 10)
y = convert(AbstractVector, @.sin(x))
f12 = Polynomials.fit(x, y, 3)
f21 = Polynomials.fit(x, y, 6)

scatter(z, ys, markersize=3, legend =:bottomright, label="data")
#plot!(xlin, f.(xlin), label = "f(x) = x3", dpi = 200)
plot!(x, f12.(x), label = "f(x) = x3 : 1", dpi = 200)
plot!(sin, label = "sin")
#plot!(xlin, f2.(xlin), label = "f(x) = x3", dpi = 200)
plot!(x, f21.(x), label = "f(x) = x3 : 2", dpi = 200)
