using LinearAlgebra

# === Modified Algorithm 11.2 from K&W ===
function traverseSimplex(A, b, B, q)
    n = size(A, 2)
    B = sort(B)
    V = sort!(setdiff(1:n, B))

    AB = A[:,B]
    d, xB = AB\A[:,V[q]], AB\b

    p, xq_ = 0, Inf
    for i in 1:Int(length(d))
        if d[i] > 0
            v = xB[i]/d[i]
            if v < xq_
                p, xq_ = i, v
            end
        end
    end
    return p, xq_
end

function stepSimplex(A, b, c, B)
    n = size(A, 2)
    B = sort!(B)
    V = sort!(setdiff(1:n, B))

    AB, AV = A[:,B], A[:,V]
    xB = AB\b
    cB = c[B]
    lambda = AB'\cB
    cV = c[V]
    muV = cV - AV'*lambda

    q, p, xq_, delta = 0, 0, Inf, Inf
    for i in 1:Int(length(muV))
        if muV[i] < 0
            p_i, x_i = traverseSimplex(A, b, B, i)
            if muV[i]*x_i < delta
                q, p, xq_, delta = i, p_i, x_i, muV[i]*x_i
            end
        end
    end
    if q == 0
        return (B, true)
    end

    j = findfirst(isequal(B[p]), B)
    B[j] = V[q]
    return (B, false)
end