using LinearAlgebra

# --- Solves a general position LP ---
function simplex(A, b, c)
    n = size(A,2)
    A_eq, b_eq, c_eq = convertToEqaulityLP(A, b, c)
    x = equalityLP_Minimize(A_eq, b_eq, c_eq)
    return x[1:n]
end

# --- Takes an inequality LP and returns the equality form and reduces ---
#     it to be linearly independent.
function convertToEqaulityLP(A, b, c)
    m = size(A, 1)

    A_eq = hcat(A, I(m))
    b_eq = b
    c_eq = vcat(c, zeros(m))

    return A_eq, b_eq, c_eq
end

# --- Converts xB to x ---
function xBTo_x(A, b, B_idx)
    n = size(A, 2)
    AB = A[:,B_idx]
    xB = AB\b

    x = []
    for i = 1:n
        if in(i, B_idx)
            term = popfirst!(xB)
            push!(x, term)
        else
            push!(x, 0.0)
        end
    end
    return x
end

# --- Optimizes an equality LP ---
function equalityLP_Minimize(A, b, c)
    m, n = size(A)

    # --- Set up initial LP ---
    z = [(b[i] >= 0 ? 1 : -1) for i in 1:m]
    Z = diagm(z)
    
    A_init = hcat(A, Z)
    b_init = b
    c_init = vcat(zeros(n), ones(m))
    B_idx  = [n+i for i = 1:m]
    
    # --- Solve initial LP to find a feasible partition ---
    B_idx = equalityLP_MinimizePartition(A_init, b_init, c_init, B_idx)

    # --- Set up final LP ---
    A_final_top = hcat(A, I(m))
    A_final_bot = hcat(zeros(m,n), I(m))
    A_final = vcat(A_final_top, A_final_bot)
    b_final = vcat(b, zeros(m))
    c_final = vcat(c, zeros(m))

    # --- Solve the final LP ---
    B_idx = equalityLP_MinimizePartition(A_final, b_final, c_final, B_idx)

    x = xBTo_x(A_final, b_final, B_idx)

    return x[1:n]
end

# --- Minimizes an LP in equality form given a partition ---
function equalityLP_MinimizePartition(A, b, c, B_idx)
    is_optimal = false
    while !is_optimal
        B_idx, is_optimal = stepVertex(A, b, c, B_idx)
    end
    return B_idx
end

# --- Finds the value of xq of the vertex taken from exchanging q with the best p ---
#     If p gives a finite value for xq, then returns 0 and Inf
function testVertex(A, b, B_idx, q)
    n = size(A, 2)
    B_idx = sort(B_idx)
    V_idx = sort!(setdiff(1:n, B_idx))

    # --- Calculate the numerator and denominator of xq ---
    AB = A[:,B_idx]
    num, denom = AB\b, AB\A[:,V_idx[q]]

    # --- Finds the value of p to be exchanged with q ---
    p, xq = 0, Inf
    for i = 1:Int(length(denom))
        if denom[i] > 0.0
            value = num[i] / denom[i]
            if value < xq
                p, xq = i, value
            end
        end
    end

    return p, xq
end

# --- Takes a single step across the best vertex. Returns true if ---
#     already at the best vertex.
function stepVertex(A, b, c, B_idx)
    # println("Starting with partition: ", B_idx)

    n = size(A, 2)
    B_idx = sort(B_idx)
    V_idx = sort!(setdiff(1:n, B_idx))

    # --- Calculate Constants ---
    AB, AV = A[:,B_idx], A[:,V_idx]
    cB = c[B_idx]
    cV = c[V_idx]
    lambda = AB'\cB
    muV = cV - AV'*lambda

    # println("Current muV = ", muV)

    # --- Test each adjacent vertex and take the one with greatest descent ---
    q, p, xq, delta = 0, 0, Inf, Inf
    for i = 1:Int(length(muV))

        if muV[i] < 0
            p_cur, x_cur = testVertex(A, b, B_idx, i)
            if muV[i]*x_cur < delta
                q, p, xq, delta = i, p_cur, x_cur, muV[i]*x_cur
            end
        end
    end

    # --- If no better vertex is found retrun partition with true ---
    if q == 0
        return B_idx, true
    end

    # --- If found a q but the xq is still Inf, the problem must be unbounded ---
    if isinf(xq)
        error("Problem is unbounded")
    end

    # --- return B with the swapped value ---
    j = findfirst(isequal(B_idx[p]), B_idx)
    B_idx[j] = V_idx[q]
    return B_idx, false
end