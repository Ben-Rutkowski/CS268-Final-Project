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
    return x[4*k-3], x[4*k-2]
end;

# --- Gets the capacity of the kth station ---
function sCap(x, k)
    return x[4*k-1]
end;

# --- Gets the building boolean of the kth station ---
function sB(x, k)
    return x[4*k]
end;

# --- Quantizes the b_k elements to be 0 or 1 ---
function correct(x)
    x_out = x
    M = length(x)/4
    for k = 1:M
        x_out[4*k] = sB(x, k)>=5 ? 1 : 0
    end
    return x_out
end


# ================ Objective Function ================
# --- OBJECTIVE FUNCTION: Total cost of all stations ---
function costFunc(M::DemandSystem, x)
    c, f = M.c, M.f
    sz, sb = NaN, NaN
    sum = 0

    M = length(x)/4
    for k = 1:M
        sz = sCap(x, k)
        sb = sB(x, k)

        sum += (c + f*sz)*sb
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
    nx, ny = nCrd(M, i)
    F = nDmd(M, i)
    sum = 0

    M = length(x)/4
    for k = 1:M
        sx, sy = sCrd(x, k)
        sz     = sCap(x, k)
        sb     = sB(x, k)
        dist = i_kDist(nx, ny, sx, sy)

        sum += (sz/dist)*sb
    end

    return F - sum
end;

# --- EQUALITY CONSTRAINT k: The kth station's b_k is a boolean --
#       (  There are M constraints of this form. One for each k = 1,...,M  )
#       (  In the form of h_k(x) = 0  )
function stationConstraint_k(x, k)
    sb = sB(x, k)
    return sb*sb - sb
end;


# ================ Reading Data ================
# --- Takes a design variable and interprates it into meaningful data ---
#       (  Returns list of coordinates and capacity for built stations only  )
#       (  Returns [ (x,y,z), ..., (x,y,z) ] for built stations only  )
function readVariable(x)
    TOL = 1e-3
    
    result = []
    M = length(x)/4
    for k = 1:M
        sx, sy = sCrd(x, k)
        sz = sCap(x, k)
        sb = sB(x, k)

        if (abs(sb-1) <= TOL)
            push!(result, (sx, sy, sz))
        end
    end

    return result
end;