# include("/Users/benjaminrutkowski/School Code Projects/CS268-Final-Project/MODEL_TWO/objective_function.jl.jl")
using Distributions
using LinearAlgebra

# ================ Gather Stats ================
# --- Calculate the mean of a list of vectors of length n ---
function calcMean(vv)
    len  = length(vv)
    n    = length(vv[1]) 
    mean = zeros(n)
    for i = 1:len
        mean += vv[i]
    end
    return mean/len
end

# --- Calculates the covariance matrix given a list of n-vecs and a mean vec ---
function calcCov(vv, mean)
    len = length(vv)
    n   = length(mean)
    cov = zeros(n,n)
    for i = 1:len
        v = vv[i] - mean
        cov += v*v'
    end
    return cov/len
end


# ================ Calculating Distributions and Sampling ================
# --- Split Distribution ---
mutable struct SplitDist
    MvNorms
    Bi

    function SplitDist(norms, bi)
        return new(norms, bi)
    end
end

function sample(D::SplitDist, n)
    M = length(D.MvNorms)
    x = rand(D.MvNorms[1], n)
    for k = 2:M
        m = rand(D.MvNorms[k], n)
        x = vcat(x, m)
    end
    b = transpose(rand(D.Bi,n)+ones(n))
    x = vcat(x, b)
    return x
end

function recalc!(D::SplitDist, v)
    M = (length(v[1])-1)/3
end

# --- Multivariate Normal Distribution ---
mutable struct MultiNormal
    Dist

    function MultiNormal(dist_in)
        new(dist_in)
    end
end

function sample(D::MultiNormal, n)
    dist = D.Dist
    return rand(dist, n)
end

# ================ Cross Entropy with Independant Distributions ================
function crossEntropyIndp(func, P, k_max, m=100, m_elite=10, B_FIX=-1)
    # for k = 1:k_max
    #     x = correctMatrix(rand(P, m), B_FIX)
    #     return x[1]
    # end
end

function subCrossEntropy(func, calcDist, P, k_max, m, m_elite)
    for k = 1:k_max
        # --- Gather m samples ---
        x = rand(P, m)
    
        # --- Find best m_elite samples
        order = sortperm( [func(x[:,i] for i in 1:m)] )
        println("Best values: ", x[order[1:m_elite]])

        # --- Recalculate distribution ---
        P = calcDist(x[order[1:m_elite]])
    end

    return P
end


# ================ Penalties ================
function quadPenalty(M::DemandSystem, x)
    N = Int(length(M.Nodes)/2)
    penalty = 0
    for i = 1:N
        gi = max(nodeConstraint_i(M, x, i), 0)
        println("Node ", i, " penalty: ", gi)
        penalty += gi*gi
    end
    return penalty
end


function inversePenalty(M::DemandSystem, x)
    N = Int(length(M.Nodes)/2)
    penalty = 0
    for i = 1:N
        println("Node ", i, " penalty: ", -1.0/nodeConstraint_i(M, x, i))
        penalty -= 1.0/nodeConstraint_i(M, x, i)
    end
    return penalty
end


# ================ Cost Function with barrier ================
# function costWithQuadPenalty(M::DemandSystem, x, p=1.0)
#     N = Int(length(M.Nodes)/2)
#     penalty = 0
#     for i = 1:N
#         gi = max(nodeConstraint_i(M, x, i), 0)
#         println("Node ", i, " penalty: ", gi)
#         penalty += gi*gi
#     end
#     cost = costFunc(M, x)
#     println("Cost: ", cost, " | Total penalty: ", p*penalty)
#     return cost + p*penalty
# end

# function costWithInverseBarrier(M::DemandSystem, x, p)
#     N = Int(length(M.Nodes)/2)
#     barrier = 0

#     for i = 1:N
#         println("Node ", i, " constraint: ", -1.0/nodeConstraint_i(M,x,i))
#         barrier -= 1/nodeConstraint_i(M,x,i)
#     end
#     barrier *= 1/p

#     println("Cost: ", costFunc(M,x), " | Penatly: ", barrier)

#     return costFunc(M,x) + barrier
# end;

# function getDist(x_elite)
#     m = length(x_elite)
#     l = length(x_elite[1])
#     mean = zeros(l)
#     for i = 1:m
#         mean += x_elite[i]
#     end
#     mean *= 1.0/float(m)
#     cov = zeros(l,l)
#     for i = 1:m
#         v = x_elite[i] - mean
#         cov += v*v'
#     end
#     cov *= 1.0/float(m)
#     # return MultivariateNormal(mean, cov)
#     return MvNormal(mean, cov)
# end