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

#=Calculates specified number of pseudoemoments from number of stored moments
    parameter: data - the data to calculate the stored moments on
    parameter: order - will calculate pseudoemoments 1 to order
    parameter: stored - will calculate "stored" moments 1 to stored
    return: SymmetricTensor of pseudoemoment tensors=#
function calcAnyPseudomom(data, order, stored)
    m = [moment(data,1), moment(data,2)]
    for i in (1:(stored-2))
        push!(m, moment(data, i+2))
    end
        
    moms2cums!(m)

    for j in (1:(order - stored))
        dim = ones(Int8, 1, stored + (j))
        for i in (1:(stored + (j)))
            dim[i] *= 2
        end
        Tuple(dim)
        c = SymmetricTensor(zeros(Float64, Tuple(dim)))
        push!(m, c)
    end
    return cums2moms(m)
end

#=Calculates a random variable that does not follow a Gaussian distribution
    parameter: std - the maximum value of the standard deviation that can be randomly generated
    parameter: samples - the number of samples generated 
    return: the number of samples of the random variable =#
function findSupport(std,  samples)
    mu1 = rand(0:1)
    sigma1 = rand(0.1:std)
    mu2 = rand(0:1)
    sigma2 = rand(0.1:std)
    d = Normal(mu1, sigma1)
    e = Normal(mu2, sigma2)
    x = rand(d, samples)
    y = rand(e, samples)
    z = x .* y
    return z
end

#=Calculates the true expectation of the function applied to two random variables
    parameter: xs - the first random variable
    parameter: samples -  the second random variable
    return: the expectation of sin(xs)*cos(ys) =#
function trueExpect(xs, ys)
    sum = 0
    for a in (1:length(xs))
        sum += @.sin(xs[a]) .* @.cos(ys[a])
    end
    return sum /= length(xs)
end

#=Calculates regression of the function over a specified range  
    parameter: min - the minimum value of the range
    parameter: max - the maximum value of the range
    parameter: l - the number of samples to take on the range
    parameter: c - the highest order of pseudoemoments calculated
    parameter: samples -  the second random variable
    return: matrix resulting from the linear regression =#    
function calcGen(min, max, l, c) 
    x = LinRange(min, max, l)
    y = LinRange(min, max, l)
    z = calcGenSupport(x, y, l)
    A = calcGenA(x, y, c)
    return leastSquaresReg(A, z)
end

#=Calculates the vector to be used in an approximation of the function over a specified range 
    parameter: x - the first LinRange
    parameter: y - the second LinRange
    parameter: l - the number of samples to take on the range
    return: a vector approximation of the function over a specified range=#
function calcGenSupport(x, y, l) 
    z = zeros(l, 1)
    for q in (1:length(x))
        z[q] = @.sin(x[q]) * @.cos(y[q])
    end
    return z
end

#=Calculates the matrix to be used in approximation of the function over a specified range 
    parameter: xVect - the first LinRange
    parameter: yVect - the second LinRange
    parameter: c - the highest order of pseudoemoments calculated
    return: coefficient matrix to use in regression=#
function calcGenA(xVect, yVect, c)
    A = Array{Float64}(undef, length(xVect), binomial(c+2, c))
    for q in (1:length(xVect))
        h = createMonomialOrderC(xVect[q], yVect[q], c)
        for r in (1:length(h))
            A[q,r] = h[r]
        end
    end
    println(typeof(A))
    return A
end

#=Calculates the least squares regression  
    parameter: A - the coefficient matrix
    parameter: z - the vector
    return: vector result of least squares regression=#
function leastSquaresReg(A, z)
    b = A\z
    return b
end

#=creates a monomial of a designated order and two variables
    parameter: a - the value of the first variable
    parameter: b - the value of the first variable
    parameter: c - the order of the monomial
    return: monomial=#
function createMonomialOrderC(a, b, c)
    @polyvar x y 
    X2 = monomials([x, y], 0:c)
    p = subs(X2, x => a, y=> b)
    return p
end

#=Puts two random variables x and y next to eachother in a matrix
    parameter: x - the first random variable
    parameter: y - the second random variable
    return: matrix=#
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

#=Selects the correct moments to use in the estimation using the result of the least squares approximation
    parameter: data - SymmetricTensor of pseudoemoment tensors
    parameter: c - the maximum order of the pseudoemoments calculated
    return: vector of selected moments=#
function getMomentsVector(data, c)
    A = Vector{Float64}(undef, binomial(c+2, c))
    B = Vector{Int}(undef, c)
    count = 1;
    for i in (0:c-1)
        if (c-i-1 >=0)
            B[i+1] =(binomial((c-i)+2, (c-i))) - (binomial(((c-1)-i)+2, ((c-1)-i)))
        else
            B[i+1] =(binomial((c-i)+2, (c-i)))
        end

    end
  
    for i in (1:c)
        for j in (1:B[i])
            A[count]= data[(c-i+1)][getMomentsIndex(c-i+1,j)...]
            count = count + 1
        end
    end
    A[binomial(c+2, c)] = 0
    
    return A

end

#=Helper method that calculates the index of the moment to select for getMomentsVector
    parameter: i - the order of the moment tensor that is being selected from
    parameter: j - the jth moment to select from the order of the moment tensor
    return: index to select=#
function getMomentsIndex(i,j)
    A =  Vector{Int}(undef, i)
    for l in (1:i)
        A[l] = 2
    end
    for l in (1:i - (j - 1))
        A[l] = 1
    end
    return A
end

#=Main function that throws everything together
    parameter: std - the maximum value of the standard deviation that can be randomly generated
    parameter: min - the minimum value of the range
    parameter: max - the maximum value of the range
    parameter: samples - the number of samples generated for each random variable
    parameter: order - will calculate pseudoemoments 1 to order
    parameter: stored - will calculate "stored" moments 1 to stored
    return: comparison of least squares approximation to true expectation=#
function main(std, min, max, samples, order, stored) 
    #random data generated
    x = findSupport(std, samples)
    y = findSupport(std, samples)
    z2 = calcZ2(x,y)

    #calc pseudoemoments and expectation
    p = calcAnyPseudomom(z2, order, stored)
    z1 = trueExpect(x, y)

    #calculate approximations
    b = calcGen(min, max, samples, order)
    arr2 = getMomentsVector(p, order)

    println(z1)
    return dot(b, arr2)
end

#actually calling the function
main(1, 0, 3, 5, 5, 3)
