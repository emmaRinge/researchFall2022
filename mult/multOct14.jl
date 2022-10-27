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

#Method 0: 2 Gaussian distributions reshaped into matrixes and multiplied
mu1 = rand(0.1:1)
sigma1 = rand(0.1:1)
mu2 = rand(0.1:1)
sigma2 = rand(0.1:1)
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
x = rand(d, 25)
y = rand(e, 25)
x2 = reshape(x, 5, 5)
y2 = reshape(x, 5, 5)
z = (x2 .* y2)

#Method 1: Gaussian distribution weighted by a Gaussian in a matrix
m1 = rand(0.1:1)
s1 = rand(0.1:1)
m2 = rand(0.1:1)
s2 = rand(0.1:1)
g = Normal(m1, s1)
h = Normal(m2, s2)

data1 = Matrix(undef, 10, 10)
g2 = rand(g, 10)
h2 = rand(h, 10)

for i in 1:size(data1, 1)
    for j in 1:size(data1, 2)
        data1[i, j] = (g2[j] .* h2[i])
    end
end


#Method 2: Non Gaussian distribution weighted by a Gaussian in a matrix
mA = rand(0.1:1)
sA = rand(0.1:3)
mB = rand(0.1:1)
sB = rand(0.1:3)
k = Normal(mA, sA)
l = Normal(mB, sB)
v = rand(k, 5)
w = rand(l, 5)
z = v .* w

mC = rand(0.1:1)
sC = rand(0.1:1)

n = Normal(mC, sC)
n2 = rand(n, 5)
data2 = Array{Float64}(undef, 5, 5)

for p in 1:size(data2, 1)
    for q in 1:size(data2, 1)
        data2[p, q] = (z[q] .* n2[p])
    end
end

moment(Matrix(data2), 3)






