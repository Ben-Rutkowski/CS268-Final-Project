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

# function printDebug(c::Constraint)
#     a = c.a
#     str = join(string.(a), "  ")
#     print(str, " | ", c.b)
# end

# function printDebug(o::Objective)
#     f = o.f
#     str = join(string.(f), "  ")
#     print(o.c, " | ", str)
# end

# --- Generates the initial guess for the optimum ---
function initialGuess(func::Objective)
    d = length(func.f)
    ans = Float64[]
    
    for k = 1:d
        if func.f[k] == 0.0
            push!(ans, 0.0)
        elseif func.f[k] < 0.0
            push!(ans, Inf)
        else
            push!(ans, -Inf)
        end
    end
    
    return ans
end

# --- Tests whether guess satisfies constraint ---
#     If every guess value is finite, then return whether or not 
#     the inequality holds. If there are infinities present, and if
#     there is one term that goes to negative inifinte return true, 
#     else return false.
function testConstraint(guess, constraint::Constraint)
    d = length(guess)
    value = 0
    pos_inf = false
    neg_inf = false
    for k = 1:d
        if isinf(guess[k]) && constraint.a[k] != 0.0
            if guess[k]*constraint.a[k] > 0.0
                pos_inf = true
            else
                neg_inf = true
            end
        elseif constraint.a[k] != 0.0
            value += constraint.a[k]*guess[k]
        end
    end

    if neg_inf
        return true
    elseif pos_inf
        return false
    else
        return value <= constraint.b
    end
end

# --- Returns a new linear program of one dimension less ---
function reduceDim(guess, func::Objective, constraints)
    last = constraints[end]
    f, a = func.f, last.a
    
    # --- Find first nonzero a in last constraint ---
    k = 0
    for j = 1:length(last.a)
        if last.a[j] != 0
            k = j
            break
        end
    end
    if k == 0
        return guess, func, constraints
    end
    ak = a[k]

    # --- Recalculate function ---
    new_fs = []
    for j = 1:length(last.a)
        if j == k
            push!(new_fs, 0.0)
        else
            push!(new_fs, f[j] - (f[k]*a[j])/ak ) 
        end
    end
    new_c = func.c + (f[k]*last.b)/ak
    new_func = Objective(new_fs, new_c)

    # print("\n- New Function: ")
    # printDebug(new_func)
    # print("\n")

    # --- Recalculate Guess ---
    new_guess = initialGuess(new_func)
    sum = 0.0
    neg_inf = false
    pos_inf = false
    for j = 1:length(last.a)
        if j != k
            coef = -a[j]/ak
            if isinf(new_guess[j]) && coef != 0.0
                if coef*new_guess[j] > 0.0
                    pos_inf = true
                else
                    neg_inf = true
                end
            elseif coef != 0.0
                sum += coef*new_guess[j]
            end
        end
    end
    if neg_inf
        if !pos_inf
            new_guess[k] = -Inf
        end
    else
        new_guess[k] = sum + last.b/ak
    end

    # print("- New Guess: ")
    # println(new_guess)

    # --- Recalculate constraints ---
    new_cons = []
    for i = 1:(length(constraints)-1)

        new_as = []
        cur_con = constraints[i]
        for j = 1:length(cur_con.a)
            if j == k
                push!(new_as, 0.0)
            else 
                term = cur_con.a[j] - (cur_con.a[k]*a[j])/ak
                push!(new_as, term)
            end
        end
        new_b = cur_con.b - (cur_con.a[k]*last.b)/ak
        new_constraint = Constraint(new_as, new_b)
        push!(new_cons, new_constraint)

        # print("- New Constraint ", i, ":\n    ")
        # printDebug(new_constraint)
        # print("\n")
        
    end

    return new_guess, new_func, new_cons
end

# --- Subroutine for Seidel ---
function subSeidel(guess, func::Objective, constraints)
    con_len = length(constraints)
    idx     = randperm(con_len)
    # idx     = 1:con_len

    # --- Test each constraint to see if it is satisfied ---
    for i = 1:con_len
        success = testConstraint(guess, constraints[idx[i]])
        # println("== CURRENT GUESS: ", guess)

        if !success
            # --- DEBUGGING ---
            # println("== VIOLATED CONSTRAINT: ", idx[i])
            # print("==    ")
            # printDebug(constraints[idx[i]])

            # println("\n\n---RECALCULATING---")
            # print("- Old Function: ")
            # printDebug(func)
            # print("\n")
            # for j = 1:i
            #     println("- Old Constraint ", j, ":")
            #     print("    ")
            #     printDebug(constraints[idx[j]])
            #     print("\n")
            # end
            # --- DEBUGGING ---
            
            new_cons = constraints[idx[1:i]]
            new_guess, new_func, new_cons = reduceDim(guess, func, new_cons)
            # println("---FINISHED CALCULATING---\n") 
            
            # println("+++ CALLING RECURSION ON ABOVE PROBLEM +++")
            guess = subSeidel(new_guess, new_func, new_cons)
            # println("+++ LEAVING RECURSION +++")
        else
            # println("success on constraint: ", i)
        end
    end

    return guess
end