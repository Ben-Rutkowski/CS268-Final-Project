# ================ Demand System ================
mutable struct DemandSystem
    Nodes # 2N length vector of xy-coords per node
    F     # N length vector of demand per node

    c # Initial cost to build a station
    f # Scaling cost to scale a station

    function DemandSystem(Nodes_in, F_in, c_in, f_in)
        new(Nodes_in, F_in, c_in, f_in)
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

# --- Gets the number of stations built ---
# function sB(x)
#     return x[end]
# end;

# --- Quantizes the b elements to be in 0,...,M ---
function correct(x)
    M, b = (length(x)-1)/3, round(x[end])

    b = b < 0 ? 0.0 :
        b > M ? float(M) :
        b
end

# --- Corrects an array of design variables ---
# function correctArray(X)
#     X_out = X
#     len = length(X)
#     for j = 1:len
#         X_out[j] = correct(X[j])
#     end
#     return X_out
# end

# --- Takes a matrix of uncorrected design variables and ---
#     returns a list of corrected variables 
function correctMatrix(M)
    x_out = []
    row, col = size(M)

    for j = 1:col
        x = correct(M[:,j])
        push!(x_out , x)
    end

    return x_out
end


# ================ Objective Function ================
# --- OBJECTIVE FUNCTION: Total cost of all built stations ---
function costFunc(M::DemandSystem, x)
    b, sum = Int(x[end]), 0

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

# --- INEQUALITY CONSTRAINT i: The ith node's demand is met by the stations ---
#       (  There are N constraints of this form. One for each i = 1,...,N  )
#       (  In the form of g_i(x) <= 0  )
function nodeConstraint_i(M::DemandSystem, x, i)
    nu, nv = nCrd(M, i)
    b, sum = Int(x[end]), 0
    
    for k = 1:b
        sx, sy = sCrd(x, k)
        sz     = sCap(x, k)

        dist = i_kDist(nu, nv, sx, sy)
        sum += sz/dist
    end
    
    return nDmd(M, i) - sum
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