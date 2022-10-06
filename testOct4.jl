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

#generate data
mu1 = 0
sigma1 = 1
mu2 = 0
sigma2 = 1
muAdd = mu1 + mu2
sigmaCom = sigma1^2 + sigma2^2
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
g = Normal(muAdd, sigmaCom)
x = rand(d, 10000000)

#Store Moments to condense data
function storeMoments(data)
    sum2 = 0
    sum3 = 0
    for x in data 
        sum2 += x^2
        sum3 += x^3
    end
    sum2 = sum2 / length(data)
    sum3 = sum3 / length(data)
    return [stats.mean(data), sum2, sum3]
end

#true expectation with original data
function trueExp(data) 
    sum = 0
    for x in data 
        sum += @.sin(x)
    end
    sum = sum / length(data)
    return sum
end

println(trueExp(x))

#first method with gaussian approximation
arr = storeMoments(x)
h = Normal(arr[1], arr[2]-(arr[1]^2))
E = expectation(h)
println(E(y -> @.sin(y)))

#second method of approximation

#third method of approximation

#compute cumulants
function calcCumulants(arr)
    cm1 = arr[1]
    cm2 = (arr[2]-(arr[1]^2))
    cm3 = arr[3] - (3 * arr[2] * arr[1]) + 2(arr[1])^3
    return [cm1, cm2, cm3]
end

#compute pseudomoments
function calcPseudomoments(arr)
    m1 = arr[1]
    m2 = arr[2] + (arr[1]^2)
    m3 = arr[3] + (3 * arr[2] * arr[1]) + 2(arr[1])^3
    m4 = (4 * arr[3] * arr[1]) + 3(arr[2])^2 + (6 * arr[2] * (arr[1])^2) + (arr[1])^4
    m5 = (10 * arr[3] * arr[2]) + (10 * arr[3] * (arr[1])^2) + (15 * ((arr[2])^2) * arr[1])
            + (10 * arr[2] * (arr[1])^3) + (arr[1])^5
    m6 = (10(arr[3]^2)) + (60 * arr[3] * arr[2] * arr[1]) + (20 * arr[3] * (arr[1])^3)
            + (15 * (arr[2])^3) + (45 * (arr[2])^2 * (arr[1])^2) + (15 * arr[2] * (arr[1])^4)
            + (arr[1])^6
    return [m1, m2, m3, m4, m5, m6]
end



#xs = range(0, 10, length = 10)
#ys = @.exp(-xs)
#f2 = poly_fit(xs, ys, 3) # degree = 2
#scatter(xs, ys, markerstrokewidth = 0, label = "Data")
#plot!(f2)
#f2 = @.sin(x)
#plot(f2)

