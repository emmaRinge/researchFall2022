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


mu1 = 0
sigma1 = 1
mu2 = 0
sigma2 = 1
muAdd = mu1 + mu2
sigmaCom = sigma1^2 + sigma2^2
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
g = Normal(muAdd, sigmaCom)
x = rand(d, 10000)

#Store Moments to condense data
function storeMoments(data)
    return [stats.mean(data), stats.std(data)]
end

#histogram(x)

f2 = @.sin(x)
plot!(f2)

arr = storeMoments(x)
h = Normal(arr[1], arr[2])
y = rand(h, 10000)
f3 = @.sin(y)
plot!(f3)

#xs = range(0, 10, length = 10)
#ys = @.exp(-xs)
#f2 = poly_fit(xs, ys, 3) # degree = 2
#scatter(xs, ys, markerstrokewidth = 0, label = "Data")
#plot!(f2)
