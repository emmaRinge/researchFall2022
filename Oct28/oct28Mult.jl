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

function findSupport(std)
    mu1 = rand(0.1:1)
    sigma1 = rand(0.1:std)
    mu2 = rand(0.1:1)
    sigma2 = rand(0.1:std)
    d = Normal(mu1, sigma1)
    e = Normal(mu2, sigma2)
    x = rand(d, 5)
    y = rand(e, 5)
    z = x .* y
    return z
end

function calcZ(x, y) 
    z = zeros(length(x), 1)
    for q in 1:length(x)
        z[q] = @.sin(x[q]) .* @.cos(y[q])
    end
    return z
end

function createMonomial(a, b)
    @polyvar x y 
    X = monomials([x, y], 0:8)
    X1 = monomials([x, y], 0:4)
    X2 = [X1[15], X1[14], X1[13], X1[12], X1[11], X1[10], X1[9], X1[8], X1[7], X1[6], X1[5], X1[4], X1[3], X1[2], X1[1], X[29], X[28], X[27], X[26], X[22], X[21], X[20], X[13], X[14], X[5]]
    p = subs(X2, x => a, y=> b)
    return p
end

function createA(xVect, yVect)
    A = Array{Float64}(undef, length(xVect), 25)
    index = 1
    for q in 1:length(xVect)
        h = createMonomial(xVect[q], yVect[q])
        for r in 1:length(h)
            temp = convert(Float64, h[r])
            A[index] = temp
            index += 1
        end
    end
    return A
end

function leastSquaresReg(A, z)
    b = A .\ z
    return b
end

function main(std) 
    x = findSupport(std)
    y = findSupport(std)
    z = calcZ(x, y)
    A = createA(x, y)
    return leastSquaresReg(A, z)
end

main(3)
#=
S = rand(2,2)
println(S)
dist = Wishart(5, S)
mat = rand(dist,n)
=#

#= function trueExp(data, yvals) 
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
=#

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
