## Author:  Pearson
## Date: April 2020
## Julia version: 1.3.1
## Purpose: Functions for EKNR(2016) Full Model - Baseline Equilibrium

function factor_price_fixpoint!(
    res_fixpoint::AbstractArray,
    guess_fixpoint::AbstractArray, # Product price, dim = NC*NS
    exos_fixpoint::NamedTuple,
    params_fixpoint::NamedTuple,
    )
    # Unpack exogenous variables and parameters
    @unpack π, T̂, d̂, ŵ, r̂ = exos_fixpoint
    @unpack β̃ᴸ, β̃ᴷ, β̃ᴹ, θ, = params_fixpoint

    # Pre-allocate memory
    p̂ = similar(guess_fixpoint)
    b̂ = similar(guess_fixpoint)
    p̂′ = zero(size(guess_fixpoint))

    # Load guess, form other variables
    p̂[:,:] = guess_fixpoint
    for n in 1:1:size(p̂,1)
        for l in 1:1:size(p̂,2)
            b̂[n,l] = ŵ[n,l]^β̃ᴸ[n,l]
            for k in 1:1:size(r̂,3)
                b̂[n,l] *= r̂[n,l,k]^β̃ᴷ[n,l,k]
            end
            for j in 1:1:size(p̂,2)
                b̂[n,l] *= p̂[n,j]^β̃ᴹ[n,l,j]
            end
        end
    end

    # Equations
    for n in 1:1:size(p̂,1)
        for j in 1:1:size(p̂,2)
            for i = 1:1:size(p̂,1)
                p̂′[n,j] += π[n,i,j]*(b̂[i,j]*d̂[n,i,j]/T̂[i,j])^-θ
            end
            p̂′[n,j] = p̂′[n,j]^(-1/θ)
            res_fixpoint[n,j] = p̂′[n,j]/p̂[n,j]-1
        end
    end

    return res_fixpoint, p̂, b̂
end
