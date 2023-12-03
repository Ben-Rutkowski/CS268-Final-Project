# ==== Model ====
mutable struct Model
    G   # Adjacency Matrix of the graph, g_ij = 0 if no edge between nodes i and j,
        # g_ij = d if there is an edge between i and j and the distance between is d
            
    n   # Number of nodes in N
    f   # Array of charging capacity for each node i
    F   # Array of demand requirement for each node i
    D   # The distance an EV can drive on a full battery
    c   # Costs at each site
    alpha # The discount factor ( See Paper )

    # === Private attributes ===
    G_hat  # Edges of G in which the distances are less than D
    N_dist # A vector of N_i^(alpha,D) ( See Paper )
end

# ==== Cost Function ====
function calcN_i(M::Model, i)
    alphaD = M.alpha*M.D
    N_i = []
    n = size(M.G)[1]
    for j = 1:n
        if G[i,j] <= alphaD
            push!(N_i, j)
        end
    end
    return N_i
end

function init!(M::Model)
    # Setting n
    row, col = size(G_hat)
    M.n = row

    # Construct G_hat
    G_hat, D = M.G, M.D
    for i = 1:row
        for j = 1:col
            if G_hat[i, j] > D
                G_hat[i, j] = 0
            end
        end
    end
    M.G_hat = G_hat

    # Construct the N_i^(alpha,D)
    N_dist = []
    N_i = []
    n = row
    for i = 1:n
        N_i = calcN_i(M, i)
        push!(N_dist, N_i)
    end

    M.N_dist = N_dist
end

function costFunction(M::Model, x)
    n = length(x)
    sum = 0
    for i = 1:n
        sum += c[i]*x[i]
    end
    return sum
end

# ==== Design Variable ====
# x = [x1 ... xn y1_11, y1_12]