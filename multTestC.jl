using Distributions
using StatsBase
import Statistics as stats
using Cumulants
using Combinatorics
using SymmetricTensors

data = reshape(collect(1.:15.),(5,3))
m = moment(data, 3)
Array(m)
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
            long3 = length(e[i][j]) 
            c = cumulants(data, long3)
            d = Array(c[long3])
            if long3 == 4
                mult *= d[1, 2, 3, 4]
            elseif long3 == 3
                mult *= d[e[i][j][1], e[i][j][2], e[i][j][3]]
            elseif long3 == 2
                mult *= d[e[i][j][1], e[i][j][2]]
            else
                mult *= d[e[i][j][1]]
            end
            j = j + 1
        end
        sum += mult
        i = i + 1
    end
    println(sum)
end

calcMoment(data, 3)
