using Distributions
using StatsBase
import Statistics as stats
using Cumulants
using Combinatorics
using SymmetricTensors

data = reshape(collect(1.:15.),(5,3))
m = moment(data, 3)
Array(m)
println(m[1,2,3])
c = cumulants(data, 3)
Array(c[2])
d = c[2]


function calcMoment(data, m)
    part = Vector(collect(1:m))
    e = collect(partitions(Array(part)))
    long = length(e)
    sum = 0
    i = 1

    while i <= long
        long2 = length(e[i])
        mult = 1
        j = 1
        while j <= long2
            l3 = length(e[i][j]) 
            c = Array(cumulants(data, l3))
            mult *= c[l3][e[i][j]...]
            j = j + 1
        end
        sum += mult
        i = i + 1
    end
    println(sum)
end

calcMoment(data, 3)
