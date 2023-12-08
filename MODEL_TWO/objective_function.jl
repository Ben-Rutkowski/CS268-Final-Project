include("/Users/benjaminrutkowski/School Code Projects/CS268-Final-Project/Cross_Entropy/simplex_algorithm.jl")

# ================ Demand System ================
mutable struct DemandSystem
    Nodes # 2N length vector of xy-coords per node
    F     # N length vector of demand per node

    c # Initial cost to build a station
    f # Scaling cost to scale a station
    r # Rolloff rate for the reach of a station

    x_min
    x_max
    y_min
    y_max
    demand_max

    function DemandSystem(Nodes_in, F_in, c_in, f_in, r_in=0.5)
        x_values = Nodes_in[1:2:end]
        y_values = Nodes_in[2:2:end]
        x_order = sortperm(x_values)
        y_order = sortperm(y_values)

        x_min, y_min = Nodes_in[2*x_order[begin]-1], Nodes[2*y_order[begin]]
        x_max, y_max = Nodes_in[2*x_order[end]-1], Nodes[2*y_order[end]]

        demand_max = maximum([float(F_in[i]) for i = 1:Int(length(Nodes_in)/2)] )

        new(
            Nodes_in, F_in, c_in, f_in, r_in,
            x_min, x_max, y_min, y_max, demand_max
        )
    end
end;

# --- Gets the coordinates of the ith node ---
function nCrd(M::DemandSystem, i)
    return M.Nodes[2*i-1], M.Nodes[2*i]
end;

# --- Gets the demand of the ith node ---
function nDmd(M::DemandSystem, i)
    return M.F[i]
end;


# ================ Design Variable ================
# --- Gets the coordinates of the kth station ---
function sCrd(x, k)
    return x[3*k-2], x[3*k-1]
end;

# --- Gets the capacity of the kth station ---
function sCap(x, k)
    return x[3*k]
end;

# --- Quantizes the b elements to be in 1,...,M ---
#     and orders stations by x-value
#     pushes the z component to be nonnegative
function correct(x, B_FIX=-1)
    M, b = (length(x)-1)/3, round(x[end])
    x_correct = Float64[]

    order = sortperm( [x[3*i-2] for i in 1:Int(M)] )
    for k = 1:Int(M)
        append!(x_correct, x[3*order[k]-2:3*order[k]])
    end

    b = B_FIX != -1 ? B_FIX : b
    b = b < 1.0    ? 1.0 :
        b > M      ? float(M) : 
        float(b)
    push!(x_correct, b)

    return x_correct
end

# --- Takes a matrix of uncorrected design variables and ---
#     returns a list of corrected variables 
function correctMatrix(M, B_FIX=-1)
    col = size(M)[2]
    matrix = correct(M[:,1], B_FIX)
    for j = 2:col
        x = correct(M[:,j], B_FIX)
        matrix = hcat(matrix, x)
    end

    return matrix
end

# --- Returns first m legal vectors in a matrix ---
function first_mLegal(D::DemandSystem, mat, m_elite)
    m, n = size(mat)
    holder = zeros(m, 1)
    count = 0
    for j = 1:n
        x = mat[:,j]
        if isLegal(D, x)
            holder = hcat(holder, x)
            count += 1
            if count == m_elite
                break
            end
        end
    end

    return holder[:,2:end]
end

# --- Returns whether or not variable is legal ---
function isLegal(D::DemandSystem, x)
    x_min, x_max, y_min, y_max = D.x_min, D.x_max, D.y_min, D.y_max
    demand_max = D.demand_max
    b = Int(x[end])

    for k = 1:b
        sx, sy = sCrd(x, k)
        sz = sCap(x, k)

        legal = x_min <= sx && sx <= x_max &&
                y_min <= sy && sy <= y_max &&
                0.0 <= sz && sz <= demand_max
        
        if !legal
            return false
        end
    end
    return true
end

# --- Takes a capacity vector and returns a new design variable with capacities ---
function zToX(z, x)
    M = Int((length(x)-1)/3)
    for k = 1:M
        x[3*k] = z[k]
    end
    return x
end


# ================ Objective Function ================
# --- OBJECTIVE FUNCTION: Total cost of all built stations ---
function costFunc(M::DemandSystem, x)
    b, sum = Int(x[end]), 0
    c, f = M.c, M.f

    for k = 1:b
        sum += c + f*sCap(x, k)
    end
    return sum
end;


# ================ Constraint Functions ================
# --- Squared distance between a node and a station ---
function i_kDist(nu, nv, sx, sy)
    return (nu-sx)*(nu-sx) + (nv-sy)*(nv-sy)
end;

# ================ Distance Coefficient ================
function i_kDistCoef(M::DemandSystem, x, i, k)
    r = M.r
    nu, nv = nCrd(M, i)
    sx, sy = sCrd(x, k)
    dist = i_kDist(nu, nv, sx, sy)

    return exp(-r*dist)
end

# --- INEQUALITY CONSTRAINT i: The ith node's demand is met by the stations ---
#       (  There are N constraints of this form. One for each i = 1,...,N  )
#       (  In the form of g_i(x) <= 0  )
function nodeConstraint_i(M::DemandSystem, x, i)
    r, b, sum = M.r, Int(x[end]), 0
    
    for k = 1:b
        sz = sCap(x, k)
        sum += sz*i_kDistCoef(sys, x, i, k)
    end
    
    return nDmd(M, i) - sum
end;

# --- INEQUALITY CONSTRAINT k: Bound constraint all positions and cost must be nonnegative ---
function variableConstraint_k(x, k)
    return -x[k]
end


# ================ Generating Variables ================
function pushToFeasible(D::DemandSystem, x)
    N = Int(length(D.Nodes)/2)
    b = Int(x[end])
    A_lp = zeros(N,b)
    for i = 1:N
        for k = 1:b
            A_lp[i,k] = -i_kDistCoef(sys, x, i, k)
        end
    end

    b_lp = [-nDmd(D, i) for i = 1:N]
    c_lp = D.f*ones(b)

    z = simplex(A_lp, b_lp, c_lp)

    return zToX(z, x)
end

function pushToFeasibleMatrix(D::DemandSystem, X)
    n = size(X, 2)
    for j = 1:n
        x = pushToFeasible(D, X[:,j])
        X = hcat(X, x)
    end
    return X[:,n:end]
end

# --- Finds the maximum capacity a single station needs to have ---
function findMaxCap(x_diff, y_diff, r)
    max_demand = maximum([nDmd(D,i) for i = 1:Int(length(D.Nodes)/2)] )
    max_dist = sqrt(x_diff*x_diff + y_diff*y_diff)
    return max_demand*exp(Dr*max_dist)
end


# ================ Penalties ================
function quadPenalty(M::DemandSystem, x, r=0.5)
    N = Int(length(M.Nodes)/2)
    penalty = 0
    for i = 1:N
        gi = max(nodeConstraint_i(M, x, i), 0)
        # println("Node ", i, " penalty: ", nodeConstraint_i(M, x, i))
        penalty += gi*gi
    end

    for k = 1:Int(length(x)-1)
        gk = max(variableConstraint_k(x, k), 0)
        # println("Variable ", k, " penalty: ", gk)
        penalty += gk*gk
    end

    return penalty
end


# ================ Reading Data ================
# --- Takes a design variable and interprates it into meaningful data ---
#       (  Returns list of coordinates and capacity for built stations only  )
#       (  Returns [ (x,y,z), ..., (x,y,z) ] for built stations only  )
#       (  Must be in valid for, i.e. b = 0,...,M  )
function readVariable(x)
    result = []
    b = Int(x[end])
    for k = 1:b
        sx, sy = sCrd(x, k)
        sz = sCap(x, k)

        push!(result, (sx, sy, sz))
    end

    return result
end;


# --- EQUALITY CONSTRAINT: The number of stations built is an integer in 1,...,M ---
#       (  In the form of h(x) = 0  )
# function stationConstraint_k(x)
#     M = (length(x)-1)/3
#     b = x[end]
#     for k = 1:M
#         b *= b-k
#     end
#     return (1/M)*b
# end;

# --- Takes a design variable and returns just the capacity ---
# function xToZ(x)
#     M = Int((length(x)-1)/3)
#     z = []
#     for k = 1:M
#         push!(z, x[3*k])
#     end
#     return z
# end

# ================ Data Processing ================
# --- Returns the width and height of the Nodes ---
# function findBounds(Nodes)
#     x_values = Nodes[1:2:end]
#     y_values = Nodes[2:2:end]
#     x_order = sortperm(x_values)
#     y_order = sortperm(y_values)
    
#     x_min = nCrd(D, x_order[begin])[1]
#     y_min = nCrd(D, y_order[begin])[2]
#     x_max = nCrd(D, x_order[end])[1]
#     y_max = nCrd(D, y_order[end])[2]

#     return x_min, x_max, y_min, y_max
# end

# function pushToLegal(D::DemandSystem, x)
#     x_min, x_max = D.x_min, D.x_max
#     y_min, y_max = D.y_min, D.y_max
#     demand_max = D.demand_max

#     M = Int((length(x)-1)/3)
#     for k = 1:M
#         x[3*k-2] = max(x_min, min(x[3*k-2], x_max))
#         x[3*k-1] = max(y_min, min(x[3*k-1], y_max))
#         x[3*k]   = max(demand_max, min(x[3*k], demand_max))
#     end

#     return x
# end