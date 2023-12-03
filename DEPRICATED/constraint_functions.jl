include("objective_function.jl")

# ==== Equalitiy Constraints ====
# ====       h(x) = 0        ====

# x = {0,1} for all i (6c)
# For each i, find the value of si = x^2 - x
# and return the value of sum of si^2 for each i
function hEq_1(M::Model, x)
    sum = 0
    term = NaN
    for i=1:M.n
        term = x[i]*x[i] - x[i]
        sum += term*term
    end
    return sum
end

# Sum of flow y_jk^i = ... (6f)
# Calculate sum y_jk^i - xixk + sum y_kl^i for each i,k in N
# return the sum of squares of each term
function hEq_2CalcTerm(M::Model, x, i, k)
    left_term, right_term = NaN, NaN
    n = M.n
    for j = 1:n
        
    end
end

function hEq_2(M::Model, x)
    n = M.n
    for i=1:n
        for j=1:n

        end
    end
end