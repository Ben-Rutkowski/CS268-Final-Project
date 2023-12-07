using LinearAlgebra

# --- Calculate x from B ---
function xFrom_xB(A, b, c, B)
    B = [1, 3]

    d = size(A, 2)
    AB = A[:,B]
    xB = AB\b

    # return xB
    # return c[B]'*xB
end

# --- Solves a LP in equality form ---
function equalityLP(A, b, c)
    m, n = size(A)
    zz = [(b[i]>=0 ? 1 : -1) for i in 1:m]
    Z = diagm(zz)
    A_init = hcat(A, Z)
    c_init = vcat(zeros(n), ones(m))
    B, isopt = [n+i for i in 1:m], false

    # Find initial feasible value
    while !isopt
        B, isopt = stepSimplex(A_init, b, c_init, B)
    end

    A_prob_top = hcat(A, I(m))
    A_prob_bot = hcat(zeros(m,n), I(m))
    A_prob = vcat(A_prob_top, A_prob_bot)
    b_prob = vcat(b, zeros(m))
    c_prob = vcat(c, zeros(m))

    isopt = false
    while !isopt
        B, isopt = stepSimplex(A_prob, b_prob, c_prob, B)
    end

    return xFrom_xB(A_prob, b_prob, c_prob, B)
end

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
    println("Checked q=", q, " to exchange.")
    println("Exchanged p=", p, " with xq=", xq_)
    return p, xq_
end

# === Modified Algorithm 11.3 from K&W ===
function stepSimplex(A, b, c, B)
    n = size(A, 2)
    B = sort!(B)
    V = sort!(setdiff(1:n, B))

    println("B: ", B, " V: ", V)

    AB, AV = A[:,B], A[:,V]
    cB = c[B]
    lambda = AB'\cB
    cV = c[V]
    muV = cV - AV'*lambda

    println("lambda: ", lambda)
    println("muV: ", muV)

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

    if isinf(xq_)
        return (B, false)
    end

    j = findfirst(isequal(B[p]), B)
    B[j] = V[q]
    return (sort(B), false)
end