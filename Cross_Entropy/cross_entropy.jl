include("/Users/benjaminrutkowski/School Code Projects/CS268-Final-Project/Cross_Entropy/new_objective_function.jl")
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
    mean = mean/col
    return mean
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
        return out + (1e-2)*I(row)
    end
    return out
end

function recalc(x_elite)
    mean = calcMeanMat(x_elite)
    cov  = calcCovMat(x_elite, mean)
    return MultivariateNormal(mean, cov)
end

function initialSamples(D::DemandSystem, M, m, B_FIX=-1)
    x_min, x_max = D.x_min, D.x_max
    y_min, y_max = D.y_min, D.y_max
    demand_max = D.demand_max

    x_values = rand(Uniform(x_min, x_max), 1, m)
    y_values = rand(Uniform(y_min, y_max), 1, m)
    dem_values = rand(Uniform(0.0, demand_max), 1, m)

    x_sample = vcat(x_values, y_values)
    x_sample = vcat(x_sample, dem_values)

    for k = 2:M
        x_values = rand(Uniform(x_min, x_max), 1, m)
        y_values = rand(Uniform(y_min, y_max), 1, m)
        dem_values = rand(Uniform(0.0, demand_max), 1, m)

        x_sample = vcat(x_sample, x_values)
        x_sample = vcat(x_sample, y_values)
        x_sample = vcat(x_sample, dem_values)   
    end

    if B_FIX == -1
        b_samples = float(rand(1:M, 1, m))
    else
        b_samples = B_FIX*ones(1, m)
    end

    # x_sample = vcat(x_sample, b*ones(1, m))
    x_sample = vcat(x_sample, b_samples)

    return x_sample
end


# ================ Cross Entropy Methods ================
function crossEntropyLegalCorrection(D::DemandSystem, func, M, k_max, B_FIX=-1, m=100, m_elite=10)
    x_sample = initialSamples(D, M, m, B_FIX)
    order    = sortperm( [func(x_sample[:,i]) for i in 1:m] )
    x_legal  = first_mLegal(D, x_sample, order, m_elite)
    P = recalc(x_legal)

    for k = 1:k_max
        x_sample = correctMatrix(rand(P, m), B_FIX)
        order = sortperm( [func(x_sample[:,i]) for i in 1:m] )

        x_legal = first_mLegal(D, x_sample, order, m_elite)
        x_legal = pushToFeasibleMatrix(D, x_legal)
        P = recalc(x_legal)
    end

    x_sample = correctMatrix(rand(P, m), B_FIX)
    order    = sortperm( [func(x_sample[:,i]) for i in 1:m] )
    x_legal  = first_mLegal(D, x_sample, order, m_elite)
    
    # return x_legal[:,1], pushToFeasible(D, x_legal[:,1])
    return x_legal[:,1]
end


function crossEntropy(D::DemandSystem, func, M, k_max, B_FIX=-1, m=100, m_elite=10)
    x_sample = initialSamples(D, M, m, B_FIX)
    order    = sortperm( [func(x_sample[:,i]) for i in 1:m] )
    x_elite  = first_mLegal(D, x_sample, order, m_elite)
    P = recalc(x_elite)

    for k = 1:k_max
        x_sample = correctMatrix(rand(P, m), B_FIX)
        order    = sortperm( [func(x_sample[:,i]) for i in 1:m] )

        x_elite = x_sample[:,order[1:m_elite]]        
        x_elite = pushToFeasibleMatrix(D, x_elite)
        P = recalc(x_elite)
    end

    x_sample = correctMatrix(rand(P, m), B_FIX)
    order    = sortperm( [func(x_sample[:,i]) for i in 1:m] )
    x_elite  = x_sample[:,order[1:m_elite]]        
    
    # return x_elite[:,1], pushToFeasible(D, x_elite[:,1])
    return x_elite[:,1]
end


# ================ Cost Function with barrier ================
# function inversePenalty(M::DemandSystem, x)
#     N = Int(length(M.Nodes)/2)
#     penalty = 0
#     for i = 1:N
#         println("Node ", i, " penalty: ", -1.0/nodeConstraint_i(M, x, i))
#         penalty -= 1.0/nodeConstraint_i(M, x, i)
#     end
#     return penalty
# end

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

# --- Creates an initial mulitvariate distribution from a ---
#     starting point. The 90% quantile is still the width of
#     the data set.
# function initialDistribution(D::DemandSystem, x_start)
#     M = Int((length(x_start)-1)/3)
#     cov_vec = Float64[]
#     x_width = D.x_max - D.x_min
#     y_width = D.y_max - D.y_min
#     d_width = D.demand_max
    
#     for k = 1:M
#         push!(cov_vec, x_width)
#         push!(cov_vec, y_width)
#         push!(cov_vec, d_width)
#     end
#     push!(cov_vec, x_start[end])

#     return MvNormal(x_start, diagm(cov_vec))
# end

# function findBestStart(D::DemandSystem, M)
#     order = sortperm(D.F)
#     x = Float64[]
#     for k = 1:M
#         top_i = pop!(order)
#         sx, sy = nCrd(D, top_i)
#         x = vcat(x, [sx, sy, 0.0])
#     end
#     push!(x, float(M))

#     return pushToFeasible(D, x)
# end

# ================ Calculating Distributions and Sampling ================
# --- Split Distribution ---
# mutable struct SplitDist
#     MvNorms
#     DiscNorm

#     function SplitDist(norms, disc)
#         return new(norms, bi)
#     end
# end

# function sample(D::SplitDist, n, B_FIX=-1)
#     M = length(D.MvNorms)
#     x = rand(D.MvNorms[1], n)
#     for k = 2:M
#         m = rand(D.MvNorms[k], n)
#         x = vcat(x, m)
#     end
#     b = transpose(rand(D.DiscNorm,n))
#     x = vcat(x, b)
#     x = correctMatrix(x, B_FIX)
#     return x
# end

# function recalc!(D::SplitDist, mat)
#     M = Int((size(mat)[1]-1)/3)
#     mvnorms = []
#     for i = 1:M
#         mean = calcMeanMat(mat[3*i-2:3*i,:])
#         cov  = calcCovMat(mat[3*i-2:3*i,:], mean)
#         dist = MvNormal(mean, cov)
#         push!(mvnorms, dist)
#     end

#     mean = calcMean(mat[end,:])
#     std  = sqrt(calcVar(mat[end,:], mean))
#     discnorm   = Normal(mean, std)

#     D.MvNorms  = mvnorms
#     D.DiscNorm = discnorm
# end;

# --- Multivariate Normal Distribution ---
# mutable struct MultiNormal
#     Dist

#     function MultiNormal(dist_in)
#         new(dist_in)
#     end
# end

# function sample(D::MultiNormal, n, B_FIX=-1)
#     dist = D.Dist
#     x = rand(dist, n)
#     return correctMatrix(x, B_FIX)
# end

# function recalc!(D::MultiNormal, mat)
#     mean = calcMeanMat(mat)
#     cov  = calcCovMat(mat, mean)
#     dist = MvNormal(mean, cov)

#     D.Dist = dist
# end;

# --- Creates an initial distribution of a multivariate normal ---
#     with the mean directly in the center, with the capacities 
#     slpit four ways.
# function initialDistribution(D::DemandSystem, M)
#     x_min, x_max = D.x_min, D.x_max
#     y_min, y_max = D.y_min, D.y_max
#     demand_max = D.demand_max

#     x_mean = 0.5*(x_max + x_min)
#     y_mean = 0.5*(y_max + y_min)

#     optimal_cap = 0
#     Nodes, r = D.Nodes, D.r
#     N = Int(length(Nodes)/2)
#     for i = 1:N
#         nu, nv = nCrd(D,i)
#         Fi     = nDmd(D,i)
#         d_ik   = i_kDist(nu, nv, x_mean, y_mean)
#         coef   = exp(-r*d_ik)
#         capacity = Fi/(M*coef)

#         if (capacity > optimal_cap)
#             optimal_cap = capacity
#         end
#     end

#     vec = Float64[]
#     for k = 1:M
#         push!(vec, x_mean)
#         push!(vec, y_mean)
#         push!(vec, demand_max)
#     end
#     push!(vec, M)

#     return MultivariateNormal(vec, diagm(vec))
# end