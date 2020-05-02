## Author:  Pearson
## Date: April 2020
## Julia version: 1.3.1
## Purpose: Solve EKNR(2016) Full Model - Baseline Equilibrium

# Precompile Packages
using Parameters
using CSV
using NLsolve
#using Plots

include("Full_Baseline_Functions.jl")
## Read data
Data = CSV.read(joinpath(@__DIR__, "..", "input", "part1_21cty.csv"))
π_data = CSV.read(joinpath(@__DIR__, "..", "input", "part2_21cty.csv"))

###################################################################
# Parameters
###################################################################
# Program Parameters
NC = 21 # Number of countries
NS = 3 # Number of sectors
NK = 2 # Number of capital types(C and D, corresponding to sector 1 and 2)
T₁ = 211
T_data = 0
T_tail = 150
T = T_data+T_tail

# Pick the final period
Data = Data[Data.date.==T₁, 3:end]
π_data = π_data[π_data.date.==T₁, 5:end]

# Preference
ρ = Data.rho[1] # Intertemporal elasticity of substitution
ψ = [Data.psiC[1] Data.psiD[1] Data.psiN[1]]

# Production
α = [Data.alphaC[1] Data.alphaD[1] 0.5] # Investment efficiency
δ = [Data.deltaC[1] Data.deltaD[1]] # Depreciation rete
θ = 2 # Trade elasticity for sector D and S
β̃ᴸ = ones(NC, NS+1)
β̃ᴷ = ones(NC, NS+1, NK)
β̃ᴹ = ones(NC, NS+1, NS)
β̃ᴸ = [Data.betaTilLC Data.betaTilLD Data.betaTilLS Data.betaTilLR] # Labor's income share
β̃ᴷ_temp = [Data.betaTilKCC Data.betaTilKDC Data.betaTilKSC Data.betaTilKRC Data.betaTilKCD Data.betaTilKDD Data.betaTilKSD Data.betaTilKRD]
β̃ᴹ_temp = [Data.betaTilICC Data.betaTilIDC Data.betaTilISC Data.betaTilIRC Data.betaTilICD Data.betaTilIDD Data.betaTilISD Data.betaTilIRD Data.betaTilICS Data.betaTilIDS Data.betaTilISS Data.betaTilIRS]
for country in eachindex(β̃ᴷ_temp[:,1])
    β̃ᴷ[country, :, :] = reshape(β̃ᴷ_temp[country, :], NS+1, NK) # Capital's income share
    β̃ᴹ[country, :, :] = reshape(β̃ᴹ_temp[country, :], NS+1, NS) # Intermediate's income share
end

###################################################################
# Guess
###################################################################
K̂ = ones(NC, NK, T)
Ŷ = ones(NC, NK, T)
guess_dynamic = ones(NC, NK, T)
guess_dynamic[1:NC,1:NK,1] = K̂[1:NC,1:NK,1]
guess_dynamic[1:NC,1:NK,2:T] = Ŷ[1:NC,1:NK,1:T-1]

###################################################################
# Preallocate memory
###################################################################
# Level variables
π = zeros(NC, NC, NS, T)
Y = zeros(NC, NS+1, T)
# X = zeros(NC, NS+1, T)
Xᶠ = zeros(NC, NS+1, T)
rK = zeros(NC, NK, T)
wL = zeros(NC, T)

# Exogenous variables
D = zeros(NC, NS+1)

###################################################################
# Initial Conditions
###################################################################
for country in 1:NC
    π[country,1:NC,1:NS,1] = Matrix(π_data[(country-1)*NC+1:country*NC, 2:4])
end
Y[1:NC,1:NS+1,1] = [Data.yC1 Data.yD1 Data.yS1 Data.yR1]
# X[1:NC,1:NS+1,1] = [Data.xC1 Data.xD1 Data.xS1 Data.xR1]
Xᶠ[1:NC,1:NS+1,1] = [Data.xFC1 Data.xFD1 Data.xFS1 Data.xFR1]
rK[1:NC,1:NK,1] = [Data.capCinc1 Data.capDinc1]
wL[1:NC,1] = Data.labinc1

###################################################################
# Exogenous Variables
###################################################################
# Time-varying variables(Hold constant in this procedure)
Â = ones(NC, NS, T)
χ̂ = ones(NC, NS, T)
d̂ = ones(NC, NC, NS, T)
L̂ = ones(NC, T)
ϕ̂ = ones(NC, T)

D[1:NC,1:NS+1,1] = [Data.DC1 Data.DD1 Data.DS1 Data.DR1]
D = repeat(D,1,1,T)
T̂ = ones(NC,NS,T)
#ψ̂ˢ = ones(NC, T)
#ϕ̂ᵂ = ones(NC*NS, T)

## Pack init conditions, exogenous variables and parameters
myinit_dynamic = @with_kw (
    π₁ = π[1:NC,1:NC,1:NS,1], Y₁ = Y[1:NC,1:NS+1,1], Xᶠ₁ = Xᶠ[1:NC,1:NS+1,1], wL₁ = wL[1:NC,1], rK₁ = rK[1:NC,1:NK,1])
myexos_dynamic = @with_kw (
    Dᴿ = D[1:NC,4,1:T], L̂ = L̂[1:NC,1:T], d̂ = d̂[1:NC,1:NC,1:NS,1:T], T̂ = T̂[1:NC,1:NS,1:T])
myparams_dynamic = @with_kw (
    T = T, NC = NC, NS = NS, NK = NK, β̃ᴸ = β̃ᴸ, β̃ᴷ = β̃ᴷ, ψ = ψ, θ = θ, β̃ᴹ = β̃ᴹ, ρ = ρ, δ = δ, α = α)
init_dynamic = myinit_dynamic()
exos_dynamic = myexos_dynamic()
params_dynamic = myparams_dynamic()

## Solve the dynamic problem
# Solve the fix point problem
println("Start to solve the dynamic problem.")
println("Run time and memory cost:")
@time results_dynamic =
    try
        results_dynamic = nlsolve(
            (res_dynamic, guess_dynamic) -> dynamic_problem!(
                res_dynamic, guess_dynamic, init_dynamic, exos_dynamic, params_dynamic),
            guess_dynamic,
            ftol=1e-6,
            method=:newton,
            autodiff=:forward,
            show_trace=false,
        )
    catch err
        if isa(err, DomainError)
            error("Failed to solve the dynamic problem, please try again.")
        end
    end

# Check Convergence
converged(results_dynamic) || error("Failed to converge in $(results_dynamic.iterations) iterations.")
println("Successfully solved the dynamic problem.\n")

# Catch Solutions
res_dynamic = similar(results_dynamic.zero)
res_dynamic, K̂, Ŷ = dynamic_problem!(
    res_dynamic, results_dynamic.zero, exos_dynamic, params_dynamic)

## Config use
show(Data[1,:], allcols=:true)
