using Random

mutable struct Constraint
    # a1z1 + ... + adzd <= b

    a # the list of d coeffecients for the d variables
    b # the constant

    function Constraint(a_in, b_in)
        return new(a_in, b_in)
    end
end

mutable struct Objective
    # c + f1z1 + ... + fdzd

    f # the list of d coeffecients for the d variables
    c # the constant

    function Objective(f_in, c_in)
        return new(f_in, c_in)
    end
end

# --- Generates the initial guess for the optimum ---
function initialGuess(func::Objective)
    f = func.f
    d, ans = length(f), Float64[]
    
    for k = 1:d
        if f[k] < 0.0
            push!(ans, Inf)
        # elseif f[k] == NaN
        #     push!(ans, NaN)
        else
            push!(ans, 0.0)
        end
    end

    return ans
end

# --- Solves a guess in a subspace, given by the function and target ---
function subGuess(func::Objective, target::Constraint)

end

# --- Returns a bounding constraint keeping values positive ---
#     The value d is dimension, and k is the varaible.
function boundConstraint(d, k)
    a = zeros(Float64, d)
    a[k] = -1.0
    return Constraint(a, 0.0)
end

# --- Tests whether guess satisfies constraint ---
#     If every guess value is finite, then return whether or not 
#     the inequality holds. If there are infinities present, and if
#     there is one term that goes to negative inifinte return true, 
#     else return false.
function testConstraint(guess, constraint::Constraint)
    a, b, d = constraint.a, constraint.b, length(guess)
    value, pos_inf, neg_inf = 0.0, false, false
    for i = 1:d
        if isinf(guess[i]) && a[i] != 0.0
            if guess[i]*a[i] > 0.0
                pos_inf = true
            else
                neg_inf = true
            end
        elseif a[i] != 0.0
            value += guess[i]*a[i]
        end
    end

    if neg_inf
        return true
    elseif pos_inf
        return false
    else 
        return value <= b
    end
end

# --- Solves a constraint given a guess and a index k ---
#     If right hand side is finite, set to xk, if has one 
function solveConstraint(guess, constraint::Constraint, k)
    a, b = constraint.a, constraint.b
    ak, xk, pos_inf, neg_inf = a[k], 0.0, false, false
    for i = 1:Int(length(a))
        coef = -a[i]/ak
        if isinf(guess[i]) && coef != 0.0
            if guess[i]*coef > 0.0
                pos_inf = true
            else
                neg_inf = true
            end
        elseif coef != 0.0
            xk += guess[i]*coef
        end
    end

    # ===DEBUGGING===
    println("-- Solving Guess: ", guess, " --")
    print("-- On Constraint: ")
    printDebug(constraint)
    print(" --\n")
    # ===DEBUGGING===

    if neg_inf && !pos_inf
        guess[k] = NaN
    else
        guess[k] = xk + b/ak
    end

    # ===DEBUGGING===
    println("-- New Guess: ", guess, " --")
    # ===DEBUGGING===

    return guess
end

# --- Merges new guess and previous guess (in the case where ---
#     xk is indeterminate).
# function mergeGuess(guess, new_guess)
#     for i = 1:Int(length(guess))
#         if new_guess[i] == NaN
#             new_guess[i] = guess[i]
#         end
#     end
#     return new_guess
# end

# --- Finds the first k with a non-zero ak ---
function findK(constraint::Constraint)
    a, k = constraint.a, 0
    for i = 1:Int(length(a))
        if a[i] != 0.0
            k = i
            break
        end
    end
    return k
end

# --- Takes a target constraint equality, a1x1 + ... + akxk = b ---
#     and solves the function to reduce the dimension. The k is the 
#     index of the first non-zero ak.
function reduceFunction(func::Objective, targetCon::Constraint, k)
    # ===DEBUGGING===
    print("-- Reducing: ")
    printDebug(func)
    println(" --")
    print("-- With target: ")
    printDebug(targetCon)
    println(" --")
    # ===DEBUGGING===

    a, b = targetCon.a, targetCon.b
    f, c = func.f, func.c
    ak, fk = a[k], f[k]
    f_red = Float64[]
    for i = 1:Int(length(a))
        if i == k
            push!(f_red, 0.0)
        else
            push!(f_red, f[i] - fk*(a[i]/ak) )
        end
    end

    new_func = Objective(f_red, c + fk*(b/ak) )

    # ===DEBUGGING===
    print("-- New Function: ")
    printDebug(new_func)
    println(" --")
    # ===DEBUGGING===

    return new_func
end

# --- Takes a target constraint equaility and solves the other constratint ---
#     reducing the dimension. The k is the first index with a non-zero ak.
function reduceConstraint(constraint::Constraint, target::Constraint)
    a, b = target.a, target.b
    c, e = constraint.a, constraint.b
    a_red = Float64[]
    ak, ck = a[k], c[k]
    for i = 1:Int(length(a))
        if i == k
            push!(a_red, 0.0)
        else
            push!(a_red, c[i] - ck*(a[i]/ak) )
        end
    end

    new_con = Constraint(a_red, e - ck*(b/ak) )

    # ===DEBUGGING===
    print("-- New Constraint: ")
    printDebug(new_con)
    println(" --")
    # ===DEBUGGING===

    return new_con
end

# --- Reduces a list of constraints ---
function reduceConstraintList(cons, target::Constraint)
    new_cons = []
    for j = 1:Int(length(cons))
        push!(new_cons, reduceConstraint(cons[j], target))
    end
    return new_cons
end

# --- Subroutine for Seidel ---
function recurseSeidel(guess, func::Objective, cons)
    # idx = randperm(Int(length(cons)))
    idx = [i for i in 1:Int(length(cons))]

    for i = 1:Int(length(cons))
        success = testConstraint(guess, cons[idx[i]])
        if !success
            target = cons[idx[i]]
            k = findK(target)

            # ===DEBUGGING===
            println(" === GUESS: ", guess, " ===")
            print(" === VIOLATED CONSTRAINT: ")
            printDebug(target)
            println(" ===\n")
            println("   CALCULATING NEW PROBLEM")
            # ===DEBUGGING===

            new_guess, new_func, new_cons = initSubSeidel(
                func, cons[idx[1:(i-1)]], target, k
            )

            # ===DEBUGGING===
            println("    ++++ Entering Recursion ++++")
            # ===DEBUGGING===
            new_guess = recurseSeidel(new_guess, new_func, new_cons)
            # ===DEBUGGING===
            println("    ++++ Leaving Recursion ++++")
            # ===DEBUGGING===

            new_guess = solveConstraint(new_guess, target, k)

            if isnan(new_guess[k])
                println("    -=-=- Setting lower bound and solving seidel again -=-=-")
                target = boundConstraint(length(func.f), k)
                new_guess, new_func, new_cons = initSubSeidel(
                    func, cons[idx[1:i]], target, k
                )

                # ===DEBUGGING===
                println("    ++++ Entering Recursion ++++")
                # ===DEBUGGING===
                new_guess = recurseSeidel(new_guess, new_func, new_cons)
                # ===DEBUGGING===
                println("    ++++ Leaving Recursion ++++")
                # ===DEBUGGING===
                println("    -=-=- Finished lower bound and solve -=-=-")
            end
            guess = new_guess
        else
            # ===DEBUGGING===
            println(" === GUESS: ", guess, " ===")
            print(" === SATISFIED CONSTRAINT: ")
            printDebug(cons[idx[i]])
            println(" ===\n")
            # ===DEBUGGING===
        end
    end

    # ===DEBUGGING===
    println("Finshed subSeidel with guess: ", guess)
    # ===DEBUGGING===
    return guess
end

# --- Subroutine for Seidel. Takes a function, constraints and target ---
#     And returns the reduced parameters for the recursion.
function initSubSeidel(func::Objective, cons, target::Constraint, k)
    println("Init Seidel:")
    new_func  = reduceFunction(func, target, k)
    new_cons  = reduceConstraintList(cons, target)
    new_guess = initialGuess(new_func)

    return new_guess, new_func, new_cons
end


# ================ DEBUGGING ================

function printDebug(c::Constraint)
    a = c.a
    str = join(string.(a), "  ")
    print(str, " | ", c.b)
end

function printDebug(o::Objective)
    f = o.f
    str = join(string.(f), "  ")
    print(o.c, " | ", str)
end

# --- Solves a constraint given a guess and a index k ---
#     If right hand side is finite, set to xk, if has one 
# function solveConstraint(guess, constraint::Constraint, func::Objective, k)
    # a, b = constraint.a, constraint.b
    # ak, xk, pos_inf, neg_inf = a[k], 0.0, false, false
    # for i = 1:Int(length(a))
    #     coef = -a[i]/ak
    #     if isinf(guess[i]) && coef != 0.0
    #         if guess[i]*coef > 0.0
    #             pos_inf = true
    #         else
    #             neg_inf = true
    #         end
    #     elseif coef != 0.0
    #         xk += guess[i]*coef
    #     end
    # end

    # println("-- Solving Guess: ", guess, " --")
    # print("-- On Constraint: ")
    # printDebug(constraint)
    # print(" --\n")

    # if pos_inf
    #     guess[k] = func.a[k] < 0.0 ? Inf : 0.0
    # elseif neg_inf
    #     guess[k] = 0.0
    # else
    #     guess[k] = xk
    # end

    # println("-- New Guess: ", guess, " --")

    # return guess
# end

# --- Returns a bounding constraint keeping values positive ---
#     The value d is dimension, and k is the varaible.
# function boundConstraint(d, k)
#     a = zeros(Float64, d)
#     a[k] = -1.0
#     return Constraint(a, 0.0)
# end