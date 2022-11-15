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
using DynamicPolynomials
using CumulantsUpdates

#calculate moments, cumulants, and pseudoemoments
function storeMoments(data)
    #println([moment(data, 1), moment(data,2), moment(data, 3)])
    return [moment(data, 1), moment(data,2), moment(data, 3)]
end


function calcCumulants(arr)
    return moms2cums!(arr)
end

#FIX
function calcPseudomoms(data, samples)
    m4 = moment(zeros(samples,samples), 4)
    m = [moment(data, 1), moment(data,2), moment(data, 3), moment(data, 4)]
    m4 = [[0 0; 0 0;;; 0 0; 0 0;;;; 0 0; 0 0;;; 0 0; 0 0];;;;]
    moms2cums!(m)
    m3 = [m[1], m[2], m[3], m4]
    cums2moms(m)
end

#FIX
function calcPseudomoms(data, samples)
    m4 = moment(zeros(samples,samples), 4)
    m = [moment(data, 1), moment(data,2), moment(data, 3), moment(data, 4)]
    m4 = [[0 0; 0 0;;; 0 0; 0 0;;;; 0 0; 0 0;;; 0 0; 0 0];;;;]
    moms2cums!(m)
    m3 = [m[1], m[2], m[3], m4]
    cums2moms(m)
end

#find the random variables for the distribution and calculate the true expectation
function findSupport(std,  samples)
    mu1 = rand(1:3)
    sigma1 = rand(0.1:std)
    mu2 = rand(1:3)
    sigma2 = rand(0.1:std)
    d = Normal(mu1, sigma1)
    e = Normal(mu2, sigma2)
    x = rand(d, samples)
    y = rand(e, samples)
    z = x .* y
    return z
end

#see which one is the right way of calculating
function trueExp(x, y)
    z = zeros(length(x), 1)
    sum = 0
    for q in (1:length(x))
        z[q] = @.sin(x[q]) .* @.cos(y[q])
        sum += @.sin(x[q]) .* @.cos(y[q])
    end
    println(sum/length(x))
    return z
end

#=function trueExpect(xs, ys)
    sum = 0
    for a in xs, b in ys 
        sum += @.sin(a) .* @.cos(b)
    end
    sum /= length(xs)
end =#

function trueExpect(xs, ys)
    sum = 0
    for a in (1:length(xs))
        sum += @.sin(xs[a]) .* @.cos(ys[a])
    end
    sum /= length(xs)
end

#x and y from min to max with l data point for sin(x) * cos(y) ex. 0, 5, 5
function calcGen(min, max, l, c) 
    x = range(min, max, length = l)
    y = range(min, max, length = l)
    z = calcGenSupport(x, y, l)
    A = calcGenA(x, y, c)
    return leastSquaresReg(A, z)
end

function calcGenSupport(x, y, l) 
    z = zeros(l, 1)
    for q in (1:l)
        z[q] = @.sin(x[q]) .* @.cos(y[q])
    end
    return z
end

function calcGenA(xVect, yVect, c)
    A = Array{Float64}(undef, length(xVect), binomial(c+2, c))
    index = 1
    for q in length(xVect)
        h = createMonomialOrderC(xVect[q], yVect[q], c)
        for r in length(h)
            temp = convert(Float64, h[r])
            A[index] = temp
            index += 1
        end
    end
    return A
end

function leastSquaresReg(A, z)
    b = A\z
    return b
end

#create a monomial of a designated order and two variables
function createMonomialOrderC(a, b, c)
    @polyvar x y 
    X2 = monomials([x, y], 0:c)
    p = subs(X2, x => a, y=> b)
    return p
end

#puts x and y next to eachother in a matrix
function calcZ2(x, y) 
    z = zeros(length(x), 2)
    for q in (1:length(x))
        r=2*q
        z[r-1] = x[q] 
        z[r] = y[q]
        q = q + 1
    end
    return z
end

function getMomentsVector(data)
    data1 = data[1]
    data2 = data[2]
    data3 = data[3]
    data4 = data[4]
    return [data4[1,1,1,1], data4[1,1,1,2], data4[1,1,2,2], 
        data4[1,2,2,2], data4[2,2,2,2], 
        data3[1,1,1], data3[1,1,2], data3[1,2,2], data3[2,2,2], 
        data2[1,1], data2[1,2], data2[2,2],
        data1[1], data1[2], 0]
end

#main function that throws everything together
function main(std, min, max, samples, order) 
    #random data generated
    x = findSupport(std, samples)
    y = findSupport(std, samples)
    z2 = calcZ2(x,y)

    #calc moments, cumulants, pseudoemoments, and expectation
    m = storeMoments(z2)
    c = calcCumulants(m)
    p = calcPseudomoms(z2, samples)
    z = trueExp(x, y)
    z1 = trueExpect(x, y)

    #calculate approximations
    b = calcGen(min, max, samples, order)

    arr = getMomentsVector(p)
    
    #println(z2)
    println(z)
    println(z1)
    return dot(b, arr)
    
end

#actually calling the function
main(3, 0, 5, 5, 4)
