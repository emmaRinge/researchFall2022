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
mu1 = rand(0:5)
sigma1 = rand(1:5)
mu2 = rand(0:5)
sigma2 = rand(1:5)
muAdd = mu1 + mu2
sigmaCom = sigma1^2 + sigma2^2
d = Normal(mu1, sigma1)
e = Normal(mu2, sigma2)
x = rand(d, 10000)
y = rand(e, 10000)
z = x .* y

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

println(trueExp(z))

#first method with gaussian approximation
arr = storeMoments(z)
h = Normal(arr[1], arr[2]-(arr[1]^2))
E = expectation(h)
println(E(y -> @.sin(y)))

#second method of approximation
println(E(x -> (@.sin(x) + (@.cos(x) * x) - ((@.sin(x)/2) * x^2) - ((@.cos(x)/6) * x^3))))




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
    m5 = ((10 * arr[3] * arr[2]) + (10 * arr[3] * (arr[1])^2) + (15 * ((arr[2])^2) * arr[1])
            + (10 * arr[2] * (arr[1])^3) + (arr[1])^5)
    m6 = ((10(arr[3]^2)) + (60 * arr[3] * arr[2] * arr[1]) + (20 * arr[3] * (arr[1])^3)
            + (15 * (arr[2])^3) + (45 * (arr[2])^2 * (arr[1])^2) + (15 * arr[2] * (arr[1])^4)
            + (arr[1])^6)
    return [m1, m2, m3, m4, m5, m6]
end

#calc pseudomoments2
function calcPseudomoments2(arr)
    m1 = arr[1]
    m2 = arr[2] + (arr[1]^2)
    m3 = (3 * arr[2] * arr[1]) + 2(arr[1])^3
    m4 = 3(arr[2])^2 + (6 * arr[2] * (arr[1])^2) + (arr[1])^4
    m5 = (15 * ((arr[2])^2) * arr[1]) + (10 * arr[2] * (arr[1])^3) + (arr[1])^5
    m6 = (15 * (arr[2])^3) + (45 * (arr[2])^2 * (arr[1])^2) + (15 * arr[2] * (arr[1])^4) + (arr[1])^6
    return [m1, m2, m3, m4, m5, m6]
end

#approx with pseudomoments
function approxWPseudomoments(arr)
    approx = (sin(arr[1]) - ((sin(arr[1])/2)*arr[2]) - ((cos(arr[1])/6)*arr[3]) + ((sin(arr[1])/24)*arr[4])
        + ((cos(arr[1])/120)*arr[5]) - ((sin(arr[1])/720)*arr[6]))
    return approx
end

function approxWPseudomoments2(arr)
    approx = sin(arr[1]) - ((sin(arr[1])/2) * arr[2])
    return approx
end

arr2 = calcCumulants(arr)
arr3 = calcPseudomoments2(arr2)
println(approxWPseudomoments2(arr3))

