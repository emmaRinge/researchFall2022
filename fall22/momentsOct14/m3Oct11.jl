using Distributions
import Statistics as stats
using DataFrames
using CSV
using LinearAlgebra
using Expectations
using Cumulants
using Plots
using Polynomials

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

function calcCumulants(arr)
    cm1 = arr[1]
    cm2 = (arr[2]-(arr[1]^2))
    cm3 = arr[3] - (3 * arr[2] * arr[1]) + 2(arr[1])^3
    return [cm1, cm2, cm3]
end

function calcPseudomoments(arr)
    m1 = arr[1]
    m2 = arr[2] + (arr[1]^2)
    m3 = arr[3] + (3 * arr[2] * arr[1]) + (arr[1])^3
    m4 = ((4 * arr[3] * arr[1]) 
        + 3(arr[2])^2 
        + (6 * arr[2] * (arr[1])^2) 
        + (arr[1])^4)

    m5 = ((10 * arr[3] * arr[2]) 
        + (10 * arr[3] * (arr[1])^2) 
        + (15 * ((arr[2])^2) * arr[1])
        + (10 * arr[2] * (arr[1])^3) 
        + (arr[1])^5)

    m6 = ((10(arr[3])^2) 
        + (60 * arr[3] * arr[2] * arr[1]) 
        + (20 * arr[3] * (arr[1])^3)
        + (15 * (arr[2])^3) 
        + (45 * (arr[2])^2 * (arr[1])^2) 
        + (15 * arr[2] * (arr[1])^4)
        + (arr[1])^6)
    return [m1, m2, m3, m4, m5, m6]
end

function trueExp(data) 
    sum = 0
    for x in data 
        sum += @.sin(x)
    end
    sum = sum / length(data)
    return sum
end

function polyExp(fit, arr)
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

function approxZ2(z)
    #xlin = range(findmin(z)[1], findmax(z)[1], length = 10)
    #ys = convert(AbstractVector, @.sin(z))
    x = range(0, 10, length = 10)
    y = convert(AbstractVector, @.sin(x))
    return f = Polynomials.fit(x, y, 3)
end
function approxZ3(z)
    #xlin = range(findmin(z)[1], findmax(z)[1], length = 10)
    #ys = convert(AbstractVector, @.sin(z))
    x = range(0, 10, length = 10)
    y = convert(AbstractVector, @.sin(x))
    return f2 = Polynomials.fit(x, y, 6)
end

function trialZ(num, std)
    a = 0
    b = 0
    c = 0
    x = 1
    data1 = Matrix(undef, num, 4) 
    while x <= num
            z = findZ(std)
            f = approxZ2(z)
            f2 = approxZ3(z)
            arr = storeMoments(z)
            cm = calcCumulants(arr)
            mom = calcPseudomoments(cm)
            a = trueExp(z)
            b = polyExp(f, arr)
            c = polyExpEx(f2, mom)
            data1[x, 1] = stats.std(z)
            data1[x, 2] = a
            data1[x, 3] = b
            data1[x, 4] = c
            a = 0
            b = 0
            c = 0
            x+=1
        end
    return data1
end

df = DataFrame(trialZ(100,1), :auto)
CSV.write("/Users/emmaringe/m3October11.csv", df)
