# include("/Users/benjaminrutkowski/School Code Projects/CS268-Final-Project/MODEL_TWO/objective_function.jl.jl")
using Distributions
using LinearAlgebra

# ================ Gather Stats ================
# --- Calculate the mean of a list of vectors of length n ---
function calcMean(v)
    n = length(v) 
    mean = 0
    for i = 1:n
        mean += v[i]
    end
    return mean/n
end

# --- Calculates the covariance matrix given a list of n-vecs and a mean vec ---
function calcVar(v, mean)
    n   = length(v)
    var = 0
    for i = 1:n
        var += (v[i]-mean)*(v[i]-mean)
    end
    return var/n
end

# --- Calculate the mean of a matrix of n rows ---
function calcMeanMat(mat)
    row, col = size(mat)
    mean = zeros(row)
    for i = 1:col
        mean += mat[:,i]
    end
    return mean/col
end

# --- Calculates the covariance matrix given a matrix with n rows and a mean vec ---
function calcCovMat(mat, mean)
    row, col = size(mat)
    cov = zeros(row,row)
    for i = 1:col
        v = mat[:,i] - mean
        cov += v*v'
    end
    out = cov/col
    try 
        cholesky(out)
    catch
        # println(out)
        # println("fail")
        return I(row)
    end
    return out
end


# ================ Calculating Distributions and Sampling ================
# --- Split Distribution ---
mutable struct SplitDist
    MvNorms
    DiscNorm

    function SplitDist(norms, disc)
        return new(norms, bi)
    end
end

function sample(D::SplitDist, n, B_FIX=-1)
    M = length(D.MvNorms)
    x = rand(D.MvNorms[1], n)
    for k = 2:M
        m = rand(D.MvNorms[k], n)
        x = vcat(x, m)
    end
    b = transpose(rand(D.DiscNorm,n))
    x = vcat(x, b)
    x = correctMatrix(x, B_FIX)
    return x
end

function recalc!(D::SplitDist, mat)
    M = Int((size(mat)[1]-1)/3)
    mvnorms = []
    for i = 1:M
        mean = calcMeanMat(mat[3*i-2:3*i,:])
        cov  = calcCovMat(mat[3*i-2:3*i,:], mean)
        dist = MvNormal(mean, cov)
        push!(mvnorms, dist)
    end

    mean = calcMean(mat[end,:])
    std  = sqrt(calcVar(mat[end,:], mean))
    discnorm   = Normal(mean, std)

    D.MvNorms  = mvnorms
    D.DiscNorm = discnorm
end;

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
function crossEntropy(func, P, k_max, m=100, m_elite=10, B_FIX=-1)
    # TODO: add sub routine to nudge the elite samples into the feasible range
    #       by calling subCrossEntropy again on a different anonymous function
end

function subCrossEntropy(func, P, k_max, m=100, m_elite=10, B_FIX=-1)
    for k = 1:k_max
        # --- Gather m samples ---
        x = sample(P, m, B_FIX)
    
        # --- Find best m_elite samples
        order = sortperm( [func(x[:,i]) for i in 1:m] )

        # --- Recalculate distribution ---
        recalc!(P, x[:,order[1:m_elite]])
    end

    return P
end


# ================ Penalties ================
function quadPenalty(M::DemandSystem, x)
    N = Int(length(M.Nodes)/2)
    penalty = 0
    for i = 1:N
        gi = max(nodeConstraint_i(M, x, i), 0)
        # println("Node ", i, " penalty: ", gi)
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