# import("objective_function.jl")

# # ================ Cross Entropy Method ================
# using Distributions
# function crossEntropy(func, x, P, k_max, m=100, m_elite=10)
#     # for k = 1:k_max
#     #     X = rand(P, m)
#     # end
# end


# # ================ Cost Function with barrier ================
# function costWithBarrier(M::DemandSystem, x, p)
#     N = Int(length(M.Nodes)/2)
#     barrier = 0

#     for i = 1:N
#         barrier -= 1/nodeConstraint_i(M,x,i)
#     end
#     barrier *= 1/p

#     return costFunc(M,x) + barrier
# end