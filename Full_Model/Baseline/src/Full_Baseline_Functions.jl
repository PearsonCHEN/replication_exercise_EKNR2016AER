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
    @unpack π, T̂, d̂, ŵ, r̂ = exos_fixpoint
    @unpack β̃ᴸ, β̃ᴷ, β̃ᴹ, θ, NC, NS, NK = params_fixpoint

    # Pre-allocate memory
    p̂ = similar(guess_fixpoint)
    b̂ = similar(guess_fixpoint)
    p̂′ = zero(size(guess_fixpoint))

    # Resolve guess
    p̂[:,:] = guess_fixpoint

    # Generate other variables
    for n in 1:NC
        for l in 1:NS
            b̂[n,l] = ŵ[n,l]^β̃ᴸ[n,l]
            for k in 1:NK
                b̂[n,l] *= r̂[n,l,k]^β̃ᴷ[n,l,k]
            end
            for j in 1:NS
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
    @unpack Ŷᴷ, Y, Xᶠ, Dᴿ, wL, L̂, rK, K̂, d̂, T̂, π = exos_static
    @unpack NC, NS, NK, β̃ᴸ, β̃ᴷ, ψ, θ, β̃ᴹ = params_static

    # Pre-allocate memory
    Ŷ = zeros(NC,NS)
    ŵ = zeros(NC)
    r̂ = zeros(NC,NK)
    K̂ = zeros(NC,NK)
    guess_fixpoint = zeros(NC,NS)
    π̂ = similar(π)
    Π = similar(π)
    Xˢ = zeros(NC)
    Y′ = zeros(NC,NS)
    RHS = zeros(NC)

    # Step 1
    # Resolve guess
    Ŷ[:,1:2] = Ŷᴷ
    Ŷ[:,3] = guess_static

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
            r̂[n,k] += β̃ᴷ[n,4,k]*(Xᶠ[n,4]-Dᴿ[n])+(Xᶠ[n,3]+Xᶠ[n,4])*ψ[n,k]/ψ[n,3]
            r̂[n,k] /= rK[n,k]*K̂[n,k]
        end
    end

    # Step 4: Solve for the goods price
    # Form the guess
    for n = 1:NC
        for l = 1:NS
            guess_fixpoint[n,l] = ŵ[n,l]^β̃ᴸ[n,l]
            for k = 1:NK
                guess_fixpoint[n,l] *= r̂[n,l,k]^β̃ᴷ[n,l,k]
            end
        end
    end

    # Pack exogenous variables and parameters for solve the price
    myexos_fixpoint = @with_kw(

    )
    myparams_fixpoint = @with_kw(

    )
    exos_fixpoint = myexos_fixpoint()
    params_fixpoint = myparams_fixpoint()

    # Solve the fix point problem
    println("Start to solve the fix point problem.")
    println("Run time and memory cost:")
    @time results_fixpoint =
        try
            results_fixpoint = nlsolve(
                (res_fixpoint, guess_fixpoint) -> factor_price_fixpoint!(
                    res_fixpoint, guess_fixpoint, exos_fixpoint, params_fixpoint),
                guess_fixpoint,
                ftol=1e-6,
                method=:newton,
                autodiff=:forward,
                show_trace=false,
            )
        catch err
            if isa(err, DomainError)
                error("Failed to solve the fix point problem, please try again.")
            end
        end

    # Check Convergence
    converged(results_fixpoint) || error("Failed to converge in $(results_fixpoint.iterations) iterations.")
    println("Successfully solved the fix point problem.\n")

    # Catch Solutions
    res_fixpoint = similar(results_fixpoint.zero)
    res_fixpoint, p̂, b̂ = factor_price_fixpoint!(
        res_fixpoint, results_fixpoint.zero, exos_fixpoint, params_fixpoint)

    # Step 5 and 6
    # Get changes of trade shares for Semidurable(S) in the appendix or Nondurable(N) in the paper
    # Form the trade share matrix at t+1
    for n = 1:NC
        for i = 1:NC
            for j = 1:NS
                π̂[n,i,j] = (b̂[i,j]*d̂[n,i,j]/T̂[i,j]/p̂[n,j])^-θ
                Π[n,i,j] = π[n,i,j]*π̂[n,i,j]
            end
        end
    end

    # Step 7: back out Ŷ, here goes the equation
    for n = 1:NC
        for j = 1:NS
            Y′[n,j] = Ŷ[n,j]*Y[n,j]
        end
    end
    Xˢ = Π[:,:,3]'\(Y′[:,3])

    nf = 0
    RHS[:] = Xᶠ[:,3]
    for n = 1:NC
        for j = 1:NS
            RHS[n] += β̃ᴹ[n,j,3]*Y′[n,j]
        end
        RHS[n] += β̃ᴹ[n,j,4]*(Xᶠ[n,4]-Dᴿ[n])
        res_static[nf+1] = Xˢ[n]/RHS[n]-1
        nf += 1
    end

    return res_static, Ŷ,
end
