## Author: Pearson
## Date: April 2020
## Purpose: Functions for Simulation of EKNR one city case
## Julia version: 1.3.1

function Fun_1Cty_Levels!(
    res_lev::AbstractArray,
    guess_lev::AbstractArray,
    init_lev::NamedTuple,
    exos_lev::NamedTuple,
    params_lev::NamedTuple,
)
    @unpack K₁ = init_lev
    @unpack χ, Aᴰ, L, Kₛₛ = exos_lev
    @unpack βᴷ, βᴸ, δ, ρ, B, T = params_lev
    Y = similar(guess_lev)
    K = similar(guess_lev)
    Y[:] = guess_lev
    K[1] = K₁
    for tt = 2:T
        K[tt] = χ[tt-1]*Aᴰ[tt-1]*B*K[tt-1]*(K[tt-1]/L[tt-1])^(βᴷ-1)*(Y[tt-1]-1)/Y[tt-1]+
            (1-δ)*K[tt-1]
        res_lev[tt-1] = abs(Y[tt]/(Y[tt-1]/(χ[tt-1]*Aᴰ[tt-1]*B*K[tt-1])*(K[tt-1]/L[tt-1])^βᴸ/ρ*
            K[tt]/((1-δ)/(χ[tt]*Aᴰ[tt]*B)*(K[tt]/L[tt])^βᴸ+βᴷ))-1)
    end
    res_lev[T] = abs(K[T]/Kₛₛ-1)
    return res_lev, K, Y
end

function Fun_1Cty_Changes!(
    res_hat::AbstractArray,
    guess_hat::AbstractArray,
    init_hat::NamedTuple,
    exos_hat::NamedTuple,
    params_hat::NamedTuple,
)
    @unpack Y₁ = init_hat
    @unpack χ̂, Âᴰ, L̂ = exos_hat
    @unpack βᴷ, βᴸ, δ, ρ, B, T = params_hat
    K̂=similar(guess_hat)
    Ŷ=similar(guess_hat)
    Y=similar(guess_hat)
    K̂[1] = guess_hat[1]
    Ŷ[1:T-1] = guess_hat[2:T]
    Y[1] = Y₁
    for tt = 2:T
        Y[tt] = Y[tt-1]*Ŷ[tt-1]
    end
    for tt = 2:T
        K̂[tt] = χ̂[tt-1]*Âᴰ[tt-1]*(K̂[tt-1]/L̂[tt-1])^(βᴷ-1)*(Y[tt-1]*Ŷ[tt-1]-1)/
            ((Y[tt-1]-1)*Ŷ[tt-1])*(K̂[tt-1]-(1-δ))+(1-δ)
        res_hat[tt-1] = abs(Ŷ[tt-1]/(K̂[tt-1]/ρ/((1-δ)/(χ̂[tt-1]*Âᴰ[tt-1])*(K̂[tt-1]/L̂[tt-1])^βᴸ+
            βᴷ*(K̂[tt-1]-(1-δ))*(Y[tt-1]/(Y[tt-1]-1))))-1)
    end
    res_hat[T] = abs(K̂[T]/1-1)
    return res_hat, K̂, Ŷ
end
