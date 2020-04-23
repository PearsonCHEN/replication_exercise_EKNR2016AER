## Author: Pearson
## Date: April 2020
## Purpose: Functions for Simulation of EKNR one city case
## Julia version: 1.3.1
function Pkgcheck()
    deps = Pkg.dependencies()
    installs = Dict{String, VersionNumber}()
    for (uuid, dep) in deps
        dep.is_direct_dep || continue
        dep.version === nothing && continue
        installs[dep.name] = dep.version
    end
    return installs
end

function Fun_1Cty_Levels!(
    err::AbstractArray,
    Y_guess::AbstractArray,
    exos_1Cty_lev::NamedTuple,
    params_1Cty_lev::NamedTuple,
)
    @unpack χ, Aᴰ, L, K₁, Kₛₛ = exos_1Cty_lev
    @unpack βᴷ, βᴸ, δ, ρ, B, T = params_1Cty_lev
    Y = Y_guess
    K = similar(Y_guess)
    K[1] = K₁
    for tt = 2:T
        K[tt] = χ[tt-1]*Aᴰ[tt-1]*B*K[tt-1]*(K[tt-1]/L[tt-1])^(βᴷ-1)*(Y[tt-1]-1)/Y[tt-1]+
            (1-δ)*K[tt-1]
        err[tt-1] = abs(Y[tt]/(Y[tt-1]/(χ[tt-1]*Aᴰ[tt-1]*B*K[tt-1])*(K[tt-1]/L[tt-1])^βᴸ/ρ*
            K[tt]/((1-δ)/(χ[tt]*Aᴰ[tt]*B)*(K[tt]/L[tt])^βᴸ+βᴷ))-1)
    end
    err[T] = abs(K[T]/Kₛₛ-1)
    return err, K
end

function Fun_1Cty_Changes!(
    err::AbstractArray,
    hat_guess::AbstractArray,
    exos_1Cty_hat::NamedTuple,
    params_1Cty_hat::NamedTuple,
)
    @unpack χ̂, Âᴰ, L̂, Y₁ = exos_1Cty_hat
    @unpack βᴷ, βᴸ, δ, ρ, B, T = params_1Cty_hat
    K̂=similar(hat_guess)
    Ŷ=similar(hat_guess)
    Y=similar(hat_guess)
    K̂[1] = hat_guess[1]
    Ŷ[1:T-1] = hat_guess[2:T]
    Y[1] = Y₁
    for tt = 2:T
        Y[tt] = Y[tt-1]*Ŷ[tt-1]
    end
    for tt = 2:T
        K̂[tt] = χ̂[tt-1]*Âᴰ[tt-1]*(K̂[tt-1]/L̂[tt-1])^(βᴷ-1)*(Y[tt-1]*Ŷ[tt-1]-1)/
            ((Y[tt-1]-1)*Ŷ[tt-1])*(K̂[tt-1]-(1-δ))+(1-δ)
        err[tt-1] = abs(Ŷ[tt-1]/(K̂[tt-1]/ρ/((1-δ)/(χ̂[tt-1]*Âᴰ[tt-1])*(K̂[tt-1]/L̂[tt-1])^βᴸ+
            βᴷ*(K̂[tt-1]-(1-δ))*(Y[tt-1]/(Y[tt-1]-1))))-1)
    end
    err[T] = abs(K̂[T]/1-1)
    return err, [K̂ Ŷ]
end
