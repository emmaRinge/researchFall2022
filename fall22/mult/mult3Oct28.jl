using Distributions
using StatsBase
import Statistics as stats
using DataFrames
using CSV
using LinearAlgebra
using Expectations
using Cumulants
using Plots
using Combinatorics
using SymmetricTensors
using Polynomials

function storeMoments(data)
    return [moment(data, 1), moment(data,2), moment(data, 3)]
end

function calcCumulants(arr)
    moms2cums!(m)
end

function calcPseudomoms(data)
    m4 = moment(zeros(5,5), 4)
    m = [moment(data, 1), moment(data,2), moment(data, 3)]
    moms2cums!(m)
    m2 = [m[1], m[2], m[3], m4]
    cums2moms(m2)
end

function trueExp(data, yvals) 
    sum = 0
    y = 1
    while y < length(yvals)
        for x in data[y] 
            sum += @.sin(x) .* @.cos(yvals[y])
        end
        y = y + 1
    end
    sum = sum / length(data[1]) / length(data)
    return sum
end

#= function trueExp2(data) 
    sum = 0
    for y in data
        for x in data[y] 
            sum += @.sin(x) .* @.cos(yvals[y])
        end
    end
    sum = sum / length(data[1])^2
    return sum
end =#

#=function polyExp(fit, arr)
    sum = 0
    sum += fit[1] * arr[1]
    sum += fit[2] * arr[2]
    sum += fit[3] * arr[3]
    return sum
end

function polyExpEx(fit, arr)
    sum = 0
    sum += fit[1] * arr[1]
    sum += fit[2] * arr[2]
    sum += fit[3] * arr[3]
    sum += fit[4] * arr[4]
    sum += fit[5] * arr[5]
    sum += fit[6] * arr[6]
    return sum
end

function findZ(std)
    mu1 = rand(0.1:1)
    sigma1 = rand(0.1:std)
    mu2 = rand(0.1:1)
    sigma2 = rand(0.1:std)
    d = Normal(mu1, sigma1)
    e = Normal(mu2, sigma2)
    x = rand(d, 1000)
    y = rand(e, 1000)
    z = x .* y
    return z
end
=#
#=
function approxZ2(z)
    xlin = range(findmin(z)[1], findmax(z)[1], length = 10)
    ys = convert(AbstractVector, @.sin(z))
    x = range(0, 10, length = 10)
    y = convert(AbstractVector, @.sin(x))
    return f = Polynomials.fit(x, y, 3)
end


function approxZ3(z)
    xlin = range(findmin(z)[1], findmax(z)[1], length = 10)
    ys = convert(AbstractVector, @.sin(z))
    x = range(0, 10, length = 10)
    y = convert(AbstractVector, @.sin(x))
    return f2 = Polynomials.fit(x, y, 6)
end
=#

function trialZ(num, std)
    mA = rand(0.1:1)
    sA = rand(0.1:std)
    mB = rand(0.1:1)
    sB = rand(0.1:std)
    k = Normal(mA, sA)
    l = Normal(mB, sB)
    v = rand(k, num)
    w = rand(l, num)
    z = v .* w

    mC = rand(0.1:std)
    sC = rand(0.1:std)

    n = Normal(mC, sC)
    n2 = rand(n, num)
    data2 = Array{Float64}(undef, num, num)

    for p in 1:size(data2, 1)
        for q in 1:size(data2, 1)
            data2[p, q] = (z[q] .* n2[p])
        end
    end
    return [n2, data2]
end
