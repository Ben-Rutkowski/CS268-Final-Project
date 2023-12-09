# include("new_objective_function.jl")

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
    
    return x_elite[:,1]
end