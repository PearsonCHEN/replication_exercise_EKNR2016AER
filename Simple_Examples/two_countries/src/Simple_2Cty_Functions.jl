## Author: Pearson
## Date: April 2020
## Purpose: Simulation of EKNR tow cities case
## Julia version: 1.3.1

#=
This function governs the equations for solving the steady state.
=#
function Fun_2Cty_SS!(
    res_ss::AbstractArray, # residuals
    guess_ss::AbstractArray, # guess
    exos_ss::NamedTuple, # exogenous variables
    params_ss::NamedTuple # parameters
)
    # Unpack exogenous variables and parameters
    @unpack χ, Aᴰ, L, dₙᵢ, Yˢₛₛ = exos_ss
    @unpack α, βᴸ, βᴷ, δ, ρ, θ, B = params_ss

    # Initialize vectors
    Yₛₛ = similar(guess_ss, 2)
    Kₛₛ = similar(guess_ss, 2)
    wₛₛ = similar(guess_ss, 2)
    rₛₛ = similar(guess_ss, 2)
    bₛₛ = similar(guess_ss, 2)
    pᴰₛₛ = similar(guess_ss, 2)
    πₙᵢₛₛ = similar(guess_ss, 2, 2)

    # Assign guess
    Yₛₛ[:] = guess_ss[:,1]
    Kₛₛ[:] = guess_ss[:,2]

    # Derive other variables in terms of guess
    Yᴰₛₛ = Yₛₛ-Yˢₛₛ
    for n = 1:size(guess_ss,1)
        wₛₛ[n] = βᴸ*Yₛₛ[n]/L[n]
        rₛₛ[n] = βᴷ*Yₛₛ[n]/Kₛₛ[n]
        bₛₛ[n] = wₛₛ[n]^βᴸ*rₛₛ[n]^βᴷ
    end

    for n = 1:size(guess_ss,1)
        pᴰₛₛ[n] = ((bₛₛ[1]*dₙᵢ[n,1]/Aᴰ[1])^-θ+(bₛₛ[2]*dₙᵢ[n,2]/Aᴰ[2])^-θ)^(-1/θ)
        for i = 1:size(guess_ss,1)
            πₙᵢₛₛ[n,i] = (bₛₛ[i]*dₙᵢ[n,i]/pᴰₛₛ[n]/Aᴰ[i])^-θ
        end
    end

    Xᴰₛₛ = πₙᵢₛₛ'\Yᴰₛₛ

    # Equations
    nf = 0
    for n = 1:size(guess_ss,1)
        res_ss[nf+1] = abs(Xᴰₛₛ[n]/pᴰₛₛ[n]/Kₛₛ[n]/(δ/χ[n])^(1/α)-1) # capital accumulation
        nf += 1
        res_ss[nf+1] = abs((1-ρ*(1-δ))*pᴰₛₛ[n]/ρ/χ[n]/(Xᴰₛₛ[n]/pᴰₛₛ[n]/Kₛₛ[n])^α/pᴰₛₛ[n]/((1-α)+
            α*βᴷ*Yₛₛ[n]/Xᴰₛₛ[n])-1) # euler equation
        nf += 1
    end

    return res_ss, Yₛₛ, Kₛₛ
end

#=
This function governs the equations for solving the problem in levels.
=#
function Fun_2Cty_Levels!(
    res_lev::AbstractArray, # residuals
    guess_lev::AbstractArray, # guess
    init_lev::NamedTuple, # initial conditions
    exos_lev::NamedTuple, # exogenous variables
    params_lev::NamedTuple # parameters
)
    # Unpack initial conditions, exogenous variables and parameters
    @unpack K₁ = init_lev
    @unpack χ, Aᴰ, L, dₙᵢ, Kₛₛ, Yˢ = exos_lev
    @unpack α, βᴷ, βᴸ, δ, ρ, θ, T = params_lev

    # Initialize vectors
    Y = similar(guess_lev)
    K = similar(guess_lev)
    w = similar(guess_lev)
    r = similar(guess_lev)
    b = similar(guess_lev)
    pᴰ = similar(guess_lev)
    πₙᵢ = similar(guess_lev, size(dₙᵢ))
    Yᴰ = similar(guess_lev)
    Xᴰ = similar(guess_lev)

    # Assign guess
    Y[:,:] = guess_lev

    # Equations and other variables
    Yᴰ = Y-Yˢ
    K[:,1] = K₁

    for n = 1:size(guess_lev,1)
        w[n,1] = βᴸ*Y[n,1]/L[n,1]
        r[n,1] = βᴷ*Y[n,1]/K[n,1]
        b[n,1] = w[n,1]^βᴸ*r[n,1]^βᴷ
    end

    for n = 1:size(guess_lev,1)
        pᴰ[n,1] = ((b[1,1]*dₙᵢ[n,1,1]/Aᴰ[1,1])^-θ+(b[2,1]*dₙᵢ[n,2,1]/Aᴰ[2,1])^-θ)^(-1/θ)
        for i = 1:size(guess_lev,1)
            πₙᵢ[n,i,1] = (b[i,1]*dₙᵢ[n,i,1]/pᴰ[n,1]/Aᴰ[i,1])^-θ
        end
    end

    Xᴰ[:,1] = (πₙᵢ[:,:,1]')\Yᴰ[:,1]

    for tt = 2:T
        for n = 1:size(guess_lev,1)
            K[n,tt] = χ[n,tt-1]*(Xᴰ[n,tt-1]/pᴰ[n,tt-1])^α*K[n,tt-1]^(1-α)+(1-δ)*K[n,tt-1]
            w[n,tt] = βᴸ*Y[n,tt]/L[n,tt]
            r[n,tt] = βᴷ*Y[n,tt]/K[n,tt]
            b[n,tt] = w[n,tt]^βᴸ*r[n,tt]^βᴷ
        end

        for n = 1:size(guess_lev,1)
            pᴰ[n,tt] = ((b[1,tt]*dₙᵢ[n,1,tt]/Aᴰ[1,tt])^-θ+(b[2,tt]*dₙᵢ[n,2,tt]/Aᴰ[2,tt])^-θ)^(-1/θ)
            for i = 1:size(guess_lev,1)
                πₙᵢ[n,i,tt] = (b[i,tt]*dₙᵢ[n,i,tt]/pᴰ[n,tt]/Aᴰ[i,tt])^-θ
            end
        end

        Xᴰ[:,tt] = (πₙᵢ[:,:,tt]')\Yᴰ[:,tt]

        for n = 1:size(guess_lev,1)
            res_lev[n,tt-1] = abs(
                (pᴰ[n,tt-1]/ρ/χ[n,tt-1]*(Xᴰ[n,tt-1]/pᴰ[n,tt-1]/K[n,tt-1])^(1-α))/
                (pᴰ[n,tt]/χ[n,tt]*(Xᴰ[n,tt]/pᴰ[n,tt]/K[n,tt])^(1-α)*
                (χ[n,tt]*(1-α)*(Xᴰ[n,tt]/pᴰ[n,tt]/K[n,tt])^α+(1-δ))+α*r[n,tt])-1
            ) # euler equations(T-1)
        end
    end

    for n = 1:size(guess_lev,1)
        res_lev[n,T] = abs((χ[n,T]*(Xᴰ[n,T]/pᴰ[n,T])^α*K[n,T]^(1-α)+(1-δ)*K[n,T])/Kₛₛ[n]-1) # terminal condition(1)
    end

    return res_lev, K, Y, Yᴰ, Xᴰ, πₙᵢ
end

function Fun_2Cty_Changes!(
    res_hat::AbstractArray, # residuals
    guess_hat::AbstractArray, # guess
    init_hat::NamedTuple, # initial conditions
    exos_hat::NamedTuple, # exogenous variables
    params_hat::NamedTuple # parameters
)
    # Unpack initial conditions, exogenous variables and parameters
    @unpack Y₁, Yᴰ₁, πₙᵢ₁ = init_hat
    @unpack χ̂, Âᴰ, L̂, d̂ₙᵢ, ϕ̂ = exos_hat
    @unpack α, βᴷ, βᴸ, δ, ρ, θ, T = params_hat

    # Initialize vectors for changes
    K̂ = similar(guess_hat)
    Ŷ = similar(guess_hat)
    ŵ = similar(guess_hat)
    r̂ = similar(guess_hat)
    b̂ = similar(guess_hat)
    p̂ᴰ = similar(guess_hat)
    π̂ₙᵢ = similar(guess_hat, size(d̂ₙᵢ))
    Ŷᴰ = similar(guess_hat)
    X̂ᴰ = similar(guess_hat)

    # Initialize vectors for levels of flow
    Y = similar(guess_hat)
    Yᴰ = similar(guess_hat)
    πₙᵢ = similar(guess_hat, size(d̂ₙᵢ))
    Xᴰ = similar(guess_hat)
    Yˢ = similar(guess_hat)

    # Assign guess and initial conditions
    K̂[:,1] = guess_hat[:,1]
    Ŷ[:,1:T-1] = guess_hat[:,2:T]
    Y[:,1] = Y₁
    Yᴰ[:,1] = Yᴰ₁
    πₙᵢ[:,:,1] = πₙᵢ₁
    Xᴰ[:,1] = (πₙᵢ[:,:,1]')\Yᴰ[:,1]
    Yˢ[:,1] = Y[:,1]-Yᴰ[:,1]
    for n = 1:1:size(guess_hat,1)
        for ts = 2:T
            Yˢ[n,ts] = Yˢ[n,ts-1].*ϕ̂[n,ts-1]
        end
    end

    # Equations and other variables
    for tt = 2:T
        for n = 1:1:size(guess_hat,1)
        # other variables and equations
            ŵ[n,tt-1] = Ŷ[n,tt-1]/L̂[n,tt-1]
            r̂[n,tt-1] = Ŷ[n,tt-1]/K̂[n,tt-1]
            b̂[n,tt-1] = ŵ[n,tt-1]^βᴸ*r̂[n,tt-1]^βᴷ
        end

        for n = 1:1:size(guess_hat,1)
            p̂ᴰ[n,tt-1] = (πₙᵢ[n,1,tt-1]*(b̂[1,tt-1]*d̂ₙᵢ[n,1,tt-1]/Âᴰ[1,tt-1])^-θ+
                πₙᵢ[n,2,tt-1]*(b̂[2,tt-1]*d̂ₙᵢ[n,2,tt-1]/Âᴰ[2,tt-1])^-θ)^(-1/θ)
        end

        for n = 1:1:size(guess_hat,1)
            for i = 1:1:size(guess_hat,1)
                π̂ₙᵢ[n,i,tt-1] = (b̂[i,tt-1]*d̂ₙᵢ[n,i,tt-1]/p̂ᴰ[n,tt-1]/Âᴰ[i,tt-1])^-θ
                πₙᵢ[n,i,tt] = π̂ₙᵢ[n,i,tt-1]*πₙᵢ[n,i,tt-1]
            end
            Y[n,tt] = Ŷ[n,tt-1]*Y[n,tt-1]
            Yᴰ[n,tt] = Y[n,tt]-Yˢ[n,tt]
        end

        Xᴰ[:,tt] = (πₙᵢ[:,:,tt]')\(Yᴰ[:,tt])

        for n = 1:1:size(guess_hat,1)
            X̂ᴰ[n,tt-1] = Xᴰ[n,tt]/Xᴰ[n,tt-1]
            K̂[n,tt] = χ̂[n,tt-1]*(X̂ᴰ[n,tt-1]/p̂ᴰ[n,tt-1]/K̂[n,tt-1])^α*(K̂[n,tt-1]-(1-δ))+(1-δ)
            res_hat[n,tt-1] = K̂[n,tt-1]/(K̂[n,tt-1]-(1-δ))/ρ/(X̂ᴰ[n,tt-1]*((1-α)+(1-δ)*
            (K̂[n,tt-1]*p̂ᴰ[n,tt-1]/X̂ᴰ[n,tt-1])^α/χ̂[n,tt-1]/(K̂[n,tt-1]-(1-δ)))+
            α*βᴷ*Y[n,tt]/Xᴰ[n,tt-1])-1 # euler equations(T-1)
        end
    end

    for n = 1:1:size(guess_hat,1)
        res_hat[n,T] = abs(K̂[n,T]/1-1) # terminal condition(1)
    end

    return res_hat, K̂, Ŷ, Ŷᴰ, X̂ᴰ, π̂ₙᵢ
end
