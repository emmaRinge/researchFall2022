using Distributions
using StatsBase
import Statistics as stats
using DataFrames
using CSV
using LinearAlgebra
using Expectations
using Cumulants
using Plots
using Polynomials
using CurveFit

mu1 = rand(0:5)
sigma1 = rand(1:5)
mu2 = rand(0:5)
sigma2 = rand(1:5)
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
x = rand(d, 100)
y = rand(e, 100)
z = x .* y

xs = range(0, 10, length = 10)
ys = @.exp(-xs)
f2 = poly_fit(xs, ys, 9) # degree = 2
  
scatter(xs, ys, markerstrokewidth = 0, label = "Data")
plot!(f2)