mutable struct Model
    G   # Adjacency Matrix of the graph, g_ij = 0 if no edge between nodes i and j,
        # g_ij = d if there is an edge between i and j and the distance between is d

    G_hat # Edges of G in which the distances are less than D

    f   # Array of charging capacity for each node i
    F   # Array of demand requirement for each node i
    D   # The distance an EV can drive on a full battery
end

function calcGHat!(M::Model)
    G_hat, D = M.G, M.D
    row, col = size(G_hat)
    for i = 1:row
        for j = 1:col
            if G_hat[i, j] > D
                G_hat[i, j] = 0
            end
        end
    end
    M.G_hat = G_hat
end