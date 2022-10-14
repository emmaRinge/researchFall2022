using Distributions
import Statistics as stats
using DataFrames
using CSV
using LinearAlgebra
using Expectations
using Cumulants
using Plots
using Polynomials
using Random

mu1 = rand(0.1:1)
sigma1 = rand(0.1:1)
mu2 = rand(0.1:1)
sigma2 = rand(0.1:1)
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
d2 = Distributions.product_distribution(d)
e2 = Distributions.product_distribution(e)
data1 = Matrix(undef, 5, 5)
data1 = Matrix(undef, 5, 5)
data2 = reshape([0,0,0,0,0,0,0,0,0],(3,3))
x = rand(d2, 25)
y = rand(e2, 25)
x2 = reshape(x, 5, 5)
y2 = reshape(x, 5, 5)
z = (x2 .* y2)
