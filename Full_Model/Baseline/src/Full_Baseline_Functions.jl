## Author:  Pearson
## Date: April 2020
## Julia version: 1.3.1
## Purpose: Functions for EKNR(2016) Full Model - Baseline Equilibrium

#=
    Subroutine 1
    - This function solves for the price of final products using equation (A.30) and (A.31).
    - Note that this function takes changes of factor price of labor and capital as given.
=#
function factor_price_fixpoint!(
    res_fixpoint::AbstractArray,
    guess_fixpoint::AbstractArray, # Product price, dim = NC*NS
    exos_fixpoint::NamedTuple,
    params_fixpoint::NamedTuple,
    )

    # Unpack exogenous variables and parameters
    @unpack ŵ, r̂, π, d̂, T̂ = exos_fixpoint
    @unpack NC, NS, NK, β̃ᴸ, β̃ᴷ, β̃ᴹ, θ = params_fixpoint

    # Pre-allocate memory
    p̂ = similar(guess_fixpoint)
    b̂ = similar(guess_fixpoint)
    p̂′ = zeros(eltype(guess_fixpoint),size(guess_fixpoint))

    # Resolve guess
    p̂[1:NC,1:NS] = guess_fixpoint

    # Generate other variables
    for n in 1:NC
        for l in 1:NS
            b̂[n,l] = ŵ[n]^β̃ᴸ[n,l]
            for k in 1:NK
                b̂[n,l] *= r̂[n,k]^β̃ᴷ[n,l,k]
            end
            for j in 1:NS
                #=if p̂[n,j]<0
                    println(guess_fixpoint)
                end=#
                b̂[n,l] *= p̂[n,j]^β̃ᴹ[n,l,j]
            end
        end
    end

    # Equations
    for n in 1:NC
        for j in 1:NS
            for i = 1:NC
                p̂′[n,j] += π[n,i,j]*(b̂[i,j]*d̂[n,i,j]/T̂[i,j])^-θ
            end
            p̂′[n,j] = p̂′[n,j]^(-1/θ)
            res_fixpoint[n,j] = p̂′[n,j]/p̂[n,j]-1
        end
    end

    return res_fixpoint, p̂, b̂
end

#=
    Subroutine 2
    - This function solves the changes of GDP in sector 3, Semidurable(S) in appendix or
    Nondurable(N) in the paper.
    - Note that this function takes changes of GDP in sector 1 and 2, Construction(C) and
    Durables(D) as given.
=#
function static_problem!(
    res_static::AbstractArray,
    guess_static::AbstractArray, # Product price, dim = NC*NS
    exos_static::NamedTuple,
    params_static::NamedTuple,
    )
    # Unpack exogenous variables and parameters
    @unpack π, Ŷᴷ, Y, Xᶠ, Dᴿ, wL, L̂, rK, K̂, d̂, T̂ = exos_static
    @unpack NC, NS, NK, β̃ᴸ, β̃ᴷ, ψ, θ, β̃ᴹ = params_static

    # Pre-allocate memory, note that Ωᵣ⋆={C,D,S}, Ωₖ={C,D}
    Ŷ = zeros(eltype(guess_static),NC,NS) # changes of sectoral GDP, (𝒩,Ωᵣ⋆)
    ŵ = zeros(eltype(guess_static),NC) # changes of labor wage, (𝒩)
    r̂ = zeros(eltype(guess_static),NC,NK) # changes of capital rental rate, (𝒩,Ωₖ)
    guess_fixpoint = zeros(eltype(guess_static),NC,NS) # goods price guess, (𝒩,Ωᵣ⋆)
    π̂ = zeros(eltype(guess_static),size(π)) # changes of trade share, (𝒩,𝒩,Ωᵣ⋆)
    π′ = similar(π̂) # level of trade share in the following period, (𝒩,𝒩,Ωᵣ⋆)
    Y′ = zeros(eltype(guess_static),NC,NS) # level of sectoral GDP, (𝒩,Ωᵣ⋆)
    Xˢ = zeros(eltype(guess_static),NC) # level of final demand for Semidurable(S), (𝒩)
    RHS = zeros(eltype(guess_static),NC) # Right hand side of Step 7, (𝒩)

    # Step 1
    # Resolve guess
    Ŷ[1:NC,1:2] = Ŷᴷ
    Ŷ[1:NC,3] = guess_static

    # Step 2 and 3
    # Generate other variables
    for n = 1:NC
        for j = 1:NS
            ŵ[n] += β̃ᴸ[n,j]*Y[n,j]*Ŷ[n,j]
        end
        ŵ[n] += β̃ᴸ[n,4]*(Xᶠ[n,4]-Dᴿ[n])
        ŵ[n] /= wL[n]*L̂[n]

        for k = 1:NK
            for j = 1:NS
                r̂[n,k] += β̃ᴷ[n,j,k]*Y[n,j]*Ŷ[n,j]
            end
            r̂[n,k] += β̃ᴷ[n,4,k]*(Xᶠ[n,4]-Dᴿ[n])+(Xᶠ[n,3]+Xᶠ[n,4])*ψ[k]/ψ[3]
            r̂[n,k] /= rK[n,k]*K̂[n,k]
        end
    end

    # Step 4: Solve for the goods price
    # Form the guess
    for n = 1:NC
        for l = 1:NS
            guess_fixpoint[n,l] = ŵ[n]^β̃ᴸ[n,l]
            for k = 1:NK
                guess_fixpoint[n,l] *= r̂[n,k]^β̃ᴷ[n,l,k]
            end
        end
    end

    # Pack exogenous variables and parameters for solve the price
    myexos_fixpoint = @with_kw (
        ŵ = ŵ[1:NC], r̂ = r̂[1:NC,1:NK], π = π[1:NC,1:NC,1:NS], d̂ = d̂[1:NC,1:NC,1:NS], T̂ = T̂[1:NC,1:NS]
    )
    myparams_fixpoint = @with_kw (
        NC = NC, NS = NS, NK = NK, β̃ᴸ = β̃ᴸ, β̃ᴷ = β̃ᴷ, β̃ᴹ = β̃ᴹ, θ = θ
    )
    exos_fixpoint = myexos_fixpoint()
    params_fixpoint = myparams_fixpoint()

    # Solve the fix point problem
    # println("Start to solve the fix point problem.")
    #println("Run time and memory cost:")
    #@time results_fixpoint =
    #    try
            results_fixpoint = nlsolve(
                (res_fixpoint, guess_fixpoint) -> factor_price_fixpoint!(
                    res_fixpoint, guess_fixpoint, exos_fixpoint, params_fixpoint),
                guess_fixpoint,
                ftol=1e-12,
                method=:anderson,
                show_trace=false,
            )
    #    catch err
    #        if isa(err, DomainError)
    #            error("Failed to solve the fix point problem, please try again.")
    #        end
    #    end

    # Check Convergence
    #converged(results_fixpoint) || error("Failed to converge in $(results_fixpoint.iterations) iterations.")
    #println("Successfully solved the fix point problem.\n")

    # Catch Solutions
    res_fixpoint = similar(results_fixpoint.zero)
    res_fixpoint, p̂, b̂ = factor_price_fixpoint!(
        res_fixpoint, results_fixpoint.zero, exos_fixpoint, params_fixpoint)

    # Step 5 and 6
    # Get changes of trade shares for Semidurable(S) in the appendix or Nondurable(N) in the paper
    # Form the trade share matrix at t+1
    for n = 1:NC
        for i = 1:NC
            for j = 2:NS
                π̂[n,i,j] = (b̂[i,j]*d̂[n,i,j]/T̂[i,j]/p̂[n,j])^-θ
                π′[n,i,j] = π[n,i,j]*π̂[n,i,j]
            end
            π̂[n,i,1] = 1.0
            π′[n,i,1] = 0.0
        end
        π′[n,n,1] = 1.0
    end

    # Step 7: back out Ŷ, here goes the equation
    for n = 1:NC
        for j = 1:NS
            Y′[n,j] = Ŷ[n,j]*Y[n,j]
        end
    end
    Xˢ[1:NC] = π′[1:NC,1:NC,3]'\(Y′[1:NC,3])

    nf = 0
    RHS[1:NC] = Xᶠ[1:NC,3]
    for n = 1:NC
        for j = 1:NS
            RHS[n] += β̃ᴹ[n,j,3]*Y′[n,j]
        end
        RHS[n] += β̃ᴹ[n,4,3]*(Xᶠ[n,4]-Dᴿ[n])
        res_static[nf+1] = Xˢ[n]/RHS[n]-1
        nf += 1
    end

    return res_static, Ŷ[1:NC,3], π′, ŵ, r̂, p̂
end

#=
    Subroutine 3
    - This function solves the whole dynamic problem.
    - Note that this function calls subroutine 2 multiple times, which calls subroutine 1
    multiple times.
=#
function dynamic_problem!(
    res_dynamic::AbstractArray,
    guess_dynamic::AbstractArray, # K̂[1:NC,1:NK,1] and Ŷ[1:NC,1:NS,1:T-1]
    init_dynamic::NamedTuple,
    exos_dynamic::NamedTuple,
    params_dynamic::NamedTuple,
    )
    # Unpack exogenous variables and parameters
    @unpack π₁, Y₁, Xᶠ₁, wL₁, rK₁ = init_dynamic
    @unpack Dᴿ, L̂, d̂, T̂ = exos_dynamic
    @unpack T, NC, NS, NK, β̃ᴸ, β̃ᴷ, ψ, θ, β̃ᴹ, ρ, δ, α = params_dynamic

    # Pre-allocate memory
    π = zeros(eltype(guess_dynamic),NC,NC,NS,T)
    Y = zeros(eltype(guess_dynamic),NC,NS,T)
    Xᶠ = zeros(eltype(guess_dynamic),NC,NS+1,T)
    wL = zeros(eltype(guess_dynamic),NC,T)
    rK = zeros(eltype(guess_dynamic),NC,NK,T)
    X = zeros(eltype(guess_dynamic),NC,NS,T)

    K̂ = ones(eltype(guess_dynamic),size(rK))
    Ŷ = ones(eltype(guess_dynamic),size(Y))
    X̂ᶠ = similar(Xᶠ)

    Ŷ_static = similar(Y)
    ŵ = similar(wL)
    r̂ = similar(rK)
    p̂ = similar(Y)

    # Assign initial conditions
    π[1:NC,1:NC,1:NS,1] = π₁
    Y[1:NC,1:NS,1] = Y₁[1:NC,1:NS]
    Xᶠ[1:NC,1:NS+1,1] = Xᶠ₁
    wL[1:NC,1] = wL₁
    rK[1:NC,1:NK,1] = rK₁
    for s in 1:NS
        X[1:NC,s,1] = π[1:NC,1:NC,s,1]'\(Y[1:NC,s,1])
    end

    # Resolve guess
    # guess_dynamic = reshape(guess_dynamic,NC,NK,T)
    K̂[1:NC,1:NK,1] = guess_dynamic[1:NC,1:NK,1]
    Ŷ[1:NC,1:NK,1:T-1] = guess_dynamic[1:NC,1:NK,2:T]

    # Evaluate Euler
    # println("Start to evaluate euler residuals.")
    for t = 1:T-1
        # Step 1
        # Calls subroutine 2
        myexos_static = @with_kw (
            π = π[1:NC,1:NC,1:NS,t],
            Ŷᴷ = Ŷ[1:NC,1:NK,t],
            Y = Y[1:NC,1:NS,t],
            Xᶠ = Xᶠ[1:NC,1:NS+1,t],
            Dᴿ = Dᴿ[1:NC,t+1],
            wL = wL[1:NC,t],
            L̂ = L̂[1:NC,t],
            rK = rK[1:NC,1:NK,t],
            K̂ = K̂[1:NC,1:NK,t],
            d̂ = d̂[1:NC,1:NC,1:NS,t],
            T̂ = T̂[1:NC,1:NS,t],
        )
        myparams_static = @with_kw (
            NC = NC, NS = NS, NK = NK, β̃ᴸ = β̃ᴸ, β̃ᴷ = β̃ᴷ, ψ = ψ, θ = θ, β̃ᴹ = β̃ᴹ
        )
        exos_static = myexos_static()
        params_static = myparams_static()
        guess_static = Ŷ[1:NC,3,t]

        # Solve the static problem
        # println(string("Start to solve the static problem in period ",t,"."))
        # println("Run time and memory cost:")
        #@time results_static =
        #    try
                results_static = nlsolve(
                    (res_static, guess_static) -> static_problem!(
                        res_static, guess_static, exos_static, params_static),
                    guess_static,
                    ftol=1e-12,
                    method=:anderson,
                    show_trace=false,
                )
        #    catch err
        #        if isa(err, DomainError)
        #            error("Failed to solve the static problem, please try again.")
        #        end
        #    end

        # Check Convergence
        # converged(results_static) || error("Failed to converge in $(results_static.iterations) iterations.")
        # println(string("Successfully solved the static problem in period ",t,".\n"))

        # Catch Solutions
        res_static = similar(results_static.zero)
        res_static, Ŷ_static[1:NC,3,t], π[1:NC,1:NC,1:NS,t+1], ŵ[1:NC,t], r̂[1:NC,1:NK,t], p̂[1:NC,1:NS,t] = static_problem!(
            res_static, results_static.zero, exos_static, params_static)

        # Update level variables Part I
        Ŷ[1:NC,3,t] = Ŷ_static[1:NC,3,t]
        for n = 1:NC
            wL[n,t+1] = wL[n,t]*ŵ[n,t]*L̂[n,t]
            for j = 1:NS
                Y[n,j,t+1] = Y[n,j,t]*Ŷ[n,j,t]
            end
            for k = 1:NK
                rK[n,k,t+1] = rK[n,k,t]*r̂[n,k,t]*K̂[n,k,t]
            end
        end
        X[1:NC,3,t+1] = π[1:NC,1:NC,3,t+1]'\(Y[1:NC,3,t+1])

        # Step 2
        # Solve for X̂ᶠ[:,1,t]
        X[1:NC,1,t+1] = Y[1:NC,1,t+1]
        for n = 1:NC
            X̂ᶠ[n,1,t] = Y[n,1,t+1]
            for j = 1:NS
                X̂ᶠ[n,1,t] -= β̃ᴹ[n,j,1]*Y[n,j,t+1]
            end
            X̂ᶠ[n,1,t] -= β̃ᴹ[n,4,1]*(Xᶠ[n,4,t]-Dᴿ[n,t+1])
            X̂ᶠ[n,1,t] /= Xᶠ[n,1,t]
        end

        # Step 3
        # Form Π[:,:,2,t+1] and get X[:,2,t+1]
        X[1:NC,2,t+1] = π[1:NC,1:NC,2,t+1]'\(Y[1:NC,2,t+1])

        # Step 4
        # Solve for X̂ᶠ[:,2,t]
        for n = 1:NC
            X̂ᶠ[n,2,t] = X[n,2,t+1]
            for j = 1:NS
                X̂ᶠ[n,2,t] -= β̃ᴹ[n,j,2]*Y[n,j,t+1]
            end
            X̂ᶠ[n,2,t] -= β̃ᴹ[n,4,2]*(Xᶠ[n,4,t]-Dᴿ[n,t+1])
            X̂ᶠ[n,2,t] /= Xᶠ[n,2,t]
        end

        # Step 5
        # Update level variables
            #= Update(d) levels(I)
                π[:,:,t+1] = π[:,:,t].*π̂[:,:,t]
                Y[:,:,t+1] = Y[:,:,t].*Ŷ[:,:,t]
                wL[:,t+1] = wL[:,t].*ŵ[:,t].*L̂[:,t]
                rK[:,:,t+1] = rK[:,:,t].*r̂[:,:,t].*K̂[:,:,t]
            =#
        # Update level variables(II)
        for n = 1:NC
            for j = 1:NK
                Xᶠ[n,j,t+1] = Xᶠ[n,j,t]*X̂ᶠ[n,j,t]
            end
            Xᶠ[n,3,t+1] = Xᶠ[n,3,t] #   *ϕ̂[n,t]*ψ̂[n,3,t]
            Xᶠ[n,4,t+1] = Xᶠ[n,4,t] #   *ϕ̂[n,t]*ψ̂[n,4,t]
        end

        # Step 6, 7
        # Use results in step 2 & 4 to evaluate Euler equation
        for n = 1:NC
            for k = 1:NK
                # Euler equation
                #=
                if p̂[n,k,t]*K̂[n,k,t]/X̂ᶠ[n,k,t]<0
                    println(n)
                    println(k)
                    println(t)
                    println(p̂[n,k,t])
                    println(K̂[n,k,t])
                    println(X̂ᶠ[n,k,t])
                    debug_t = t+1
                    colnames = [[string("X[1:NC,",i,",",debug_t,"]") for i in 1:3]; [string("Y[1:NC,",i,",",debug_t,"]") for i in 1:3];[string("XF[1:NC,",i,",",debug_t,"]") for i in 1:4]]
                    colnames = [Symbol(names) for names in colnames]
                    π_colnames = [[string("pi[1:NC,",i,",1,",debug_t,"]") for i in 1:21];[string("pi[1:NC,",i,",2,",debug_t,"]") for i in 1:21];[string("pi[1:NC,",i,",3,",debug_t,"]") for i in 1:21]]
                    π_colnames = [Symbol(names) for names in π_colnames]
                    CSV.write(string("debug_t",debug_t,".csv"),Tables.table([X[1:NC,1:NS,debug_t] Y[1:NC,1:NS,debug_t] Xᶠ[1:NC,1:NS+1,debug_t]];header=colnames))
                    CSV.write(string("pi_debug_t",debug_t,".csv"),Tables.table([π[1:NC,1:NC,1,debug_t] π[1:NC,1:NC,2,debug_t] π[1:NC,1:NC,3,debug_t]];header=π_colnames))
                end
                =#
                res_dynamic[n,k,t] = K̂[n,k,t]/(K̂[n,k,t]-(1-δ[k]))/ρ/
                    (α[k]*rK[n,k,t+1]/Xᶠ[n,k,t]+X̂ᶠ[n,k,t]*((1-α[k])+(p̂[n,k,t]*K̂[n,k,t]/X̂ᶠ[n,k,t])^α[k]*((1-δ[k])/(K̂[n,k,t]-(1-δ[k])))/χ̂[n,k,t]))-1.0
                # Update K̂ᵏₜ₊₁
                K̂[n,k,t+1] = χ̂[n,k,t]*(X̂ᶠ[n,k,t]/p̂[n,k,t]/K̂[n,k,t])^α[k]*(K̂[n,k,t]-(1-δ[k]))+(1-δ[k])
            end
        end

        # Step 8: Next iteration, from t -> T-1
    end
    # println("Successfully evaluate euler residuals.")

    # Step 8
    # Evaluate terminal conditions
    # println("Start to evaluate terminal conditions.")
    for n = 1:NC
        for k = 1:NK
            res_dynamic[n,k,T] = K̂[n,k,T]/1-1
        end
    end
    # println("Successfully evaluate terminal conditions.")
    #=
    colnames = [[string("Y_hat[1:NC,",i,",",T-1,"]") for i in 1:3]; [string("K_hat[1:NC,",i,",",T,"]") for i in 1:2]]
    colnames = [Symbol(names) for names in colnames]
    CSV.write(string("dynamic",T,".csv"),Tables.table([Ŷ[1:NC,1:NS,T-1] K̂[1:NC,1:NK,T]];header=colnames))
    =#
    return res_dynamic, K̂, Ŷ
end
