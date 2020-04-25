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
    wₛₛ = βᴸ.*Yₛₛ./L
    rₛₛ = βᴷ.*Yₛₛ./Kₛₛ
    bₛₛ = wₛₛ.^βᴸ.*rₛₛ.^βᴷ
    pᴰₛₛ[:] = ((bₛₛ[1].*dₙᵢ[:,1]./Aᴰ[1]).^-θ+(bₛₛ[2].*dₙᵢ[:,2]./Aᴰ[2]).^-θ).^(-1/θ)
    πₙᵢₛₛ[:,1] = (bₛₛ[1].*dₙᵢ[:,1]./pᴰₛₛ[:]./Aᴰ[1]).^-θ
    πₙᵢₛₛ[:,2] = (bₛₛ[2].*dₙᵢ[:,2]./pᴰₛₛ[:]./Aᴰ[2]).^-θ
    Xᴰₛₛ = πₙᵢₛₛ'\Yᴰₛₛ

    # Equations
    res_ss[1:2] = abs.(Xᴰₛₛ./pᴰₛₛ./Kₛₛ./(δ./χ).^(1/α).-1)' # capital accumulation
    res_ss[3:4] = abs.((1-ρ*(1-δ))*pᴰₛₛ./ρ./χ./(Xᴰₛₛ./pᴰₛₛ./Kₛₛ).^α./pᴰₛₛ./((1-α).+
        α*βᴷ.*Yₛₛ./Xᴰₛₛ).-1)' # euler equation
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
    w[:,1] = βᴸ.*Y[:,1]./L[:,1]
    r[:,1] = βᴷ.*Y[:,1]./K[:,1]
    b[:,1] = w[:,1].^βᴸ.*r[:,1].^βᴷ
    pᴰ[:,1] = ((b[1,1].*dₙᵢ[:,1,1]./Aᴰ[1,1]).^-θ.+(b[2,1].*dₙᵢ[:,2,1]./Aᴰ[2,1]).^-θ).^
        (-1/θ)
    πₙᵢ[:,1,1] = (b[1,1].*dₙᵢ[:,1,1]./pᴰ[:,1]./Aᴰ[1,1]).^-θ
    πₙᵢ[:,2,1] = (b[2,1].*dₙᵢ[:,2,1]./pᴰ[:,1]./Aᴰ[2,1]).^-θ
    Xᴰ[:,1] = (πₙᵢ[:,:,1]')\Yᴰ[:,1]
    for tt = 2:T
        K[:,tt] = χ[:,tt-1].*(Xᴰ[:,tt-1]./pᴰ[:,tt-1]).^α.*K[:,tt-1].^(1-α).+(1-δ).*K[:,tt-1]
        w[:,tt] = βᴸ.*Y[:,tt]./L[:,tt]
        r[:,tt] = βᴷ.*Y[:,tt]./K[:,tt]
        b[:,tt] = w[:,tt].^βᴸ.*r[:,tt].^βᴷ
        pᴰ[:,tt] = ((b[1,tt].*dₙᵢ[:,1,tt]./Aᴰ[1,tt]).^-θ.+
            (b[2,tt].*dₙᵢ[:,2,tt]./Aᴰ[2,tt]).^-θ).^(-1/θ)
        πₙᵢ[:,1,tt] = (b[1,tt].*dₙᵢ[:,1,tt]./pᴰ[:,tt]./Aᴰ[1,tt]).^-θ
        πₙᵢ[:,2,tt] = (b[2,tt].*dₙᵢ[:,2,tt]./pᴰ[:,tt]./Aᴰ[2,tt]).^-θ
        Xᴰ[:,tt] = (πₙᵢ[:,:,tt]')\Yᴰ[:,tt]
        res_lev[:,tt-1] = abs.(
            (pᴰ[:,tt-1]./ρ./χ[:,tt-1].*(Xᴰ[:,tt-1]./pᴰ[:,tt-1]./K[:,tt-1]).^(1-α))./
            (pᴰ[:,tt]./χ[:,tt].*(Xᴰ[:,tt]./pᴰ[:,tt]./K[:,tt]).^(1-α).*
            (χ[:,tt].*(1-α).*(Xᴰ[:,tt]./pᴰ[:,tt]./K[:,tt]).^α.+(1-δ))+α.*r[:,tt]).-1
        ) # euler equations(T-1)
    end
    res_lev[:,T] = abs.((χ[:,T].*(Xᴰ[:,T]./pᴰ[:,T]).^α.*K[:,T].^(1-α)+(1-δ).*K[:,T])./
        Kₛₛ[:].-1) # terminal condition(1)
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
    for ts = 2:T
        Yˢ[:,ts] = Yˢ[:,ts-1].*ϕ̂[:,ts-1]
    end

    # Equations and other variables
    for tt = 2:T
        # other variables and equations
        ŵ[:,tt-1] = Ŷ[:,tt-1]./L̂[:,tt-1]
        r̂[:,tt-1] = Ŷ[:,tt-1]./K̂[:,tt-1]
        b̂[:,tt-1] = ŵ[:,tt-1].^βᴸ.*r̂[:,tt-1].^βᴷ
        p̂ᴰ[:,tt-1] = (πₙᵢ[:,1,tt-1].*(b̂[1,tt-1].*d̂ₙᵢ[:,1,tt-1]./Âᴰ[1,tt-1]).^-θ.+
            πₙᵢ[:,2,tt-1].*(b̂[2,tt-1].*d̂ₙᵢ[:,2,tt-1]./Âᴰ[2,tt-1]).^-θ).^(-1/θ)
        π̂ₙᵢ[:,1,tt-1] = (b̂[1,tt-1].*d̂ₙᵢ[:,1,tt-1]./p̂ᴰ[:,tt-1]./Âᴰ[1,tt-1]).^-θ
        π̂ₙᵢ[:,2,tt-1] = (b̂[2,tt-1].*d̂ₙᵢ[:,2,tt-1]./p̂ᴰ[:,tt-1]./Âᴰ[2,tt-1]).^-θ

        # update intertemporal variables
        πₙᵢ[:,:,tt] = π̂ₙᵢ[:,:,tt-1].*πₙᵢ[:,:,tt-1]
        Y[:,tt] = Ŷ[:,tt-1].*Y[:,tt-1]
        Yᴰ[:,tt] = Y[:,tt]-Yˢ[:,tt]
        Xᴰ[:,tt] = (πₙᵢ[:,:,tt]')\(Yᴰ[:,tt])
        X̂ᴰ[:,tt-1] = Xᴰ[:,tt]./Xᴰ[:,tt-1]
        K̂[:,tt] = χ̂[:,tt-1].*(X̂ᴰ[:,tt-1]./p̂ᴰ[:,tt-1]./K̂[:,tt-1]).^α.*(K̂[:,tt-1].-(1-δ)).+(1-δ)

        # euler equations
        res_hat[:,tt-1] = K̂[:,tt-1]./(K̂[:,tt-1].-(1-δ))./ρ./(X̂ᴰ[:,tt-1].*((1-α).+(1-δ).*
            (K̂[:,tt-1].*p̂ᴰ[:,tt-1]./X̂ᴰ[:,tt-1]).^α./χ̂[:,tt-1]./(K̂[:,tt-1].-(1-δ))).+
            α*βᴷ.*Y[:,tt]./Xᴰ[:,tt-1]).-1
            # euler equations(T-1)
    end
    res_hat[:,T] = abs.(K̂[:,T]./1 .-1) # terminal condition(1)
    return res_hat, K̂, Ŷ, Ŷᴰ, X̂ᴰ, π̂ₙᵢ
end
