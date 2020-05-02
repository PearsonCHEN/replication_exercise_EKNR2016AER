## Author:  Pearson
## Date: April 2020
## Julia version: 1.3.1
## Purpose: Solve EKNR(2016) Full Model - Baseline Equilibrium

# Precompile Packages
using Parameters
using CSV
using NLsolve
#using Plots

# Program Parameters
T₀ = 211
T_data = 0
T_tail = 150
T = T_data+T_tail

# Read data
Data = CSV.read(joinpath(@__DIR__, "..", "input", "part1_21cty.csv"))
π_data = CSV.read(joinpath(@__DIR__, "..", "input", "part2_21cty.csv"))
Data = Data[Data.date.==T₀, 3:end]
π_data = π_data[π_data.date.==T₀, 5:end]

###################################################################
# Parameters
###################################################################
NC = size(Data, 1) # Number of countries
NS = 3 # Number of sectors
NK = 2 # Number of capital types(C and D, corresponding to sector 1 and 2)

# Preference
ρ = Data.rho[1] # Intertemporal elasticity of substitution

# Production
α = [Data.alphaC[1] Data.alphaD[1] 0.5] # Investment efficiency
δ = [Data.deltaC[1] Data.deltaD[1]] # Depreciation rete
θ = [2 2] # Trade elasticity for sector D and S
β̃ᴸ = ones(NC, NS+1)
β̃ᴷ = ones(NC, NS+1, NS-1)
β̃ᴹ = ones(NC, NS, NS)
β̃ᴸ = [Data.betaTilLC Data.betaTilLD Data.betaTilLS Data.betaTilLR] # Labor's income share
β̃ᴷ_temp = [Data.betaTilKCC Data.betaTilKDC Data.betaTilKSC Data.betaTilKRC Data.betaTilKCD Data.betaTilKDD Data.betaTilKSD Data.betaTilKRD]
β̃ᴹ_temp = [Data.betaTilICC Data.betaTilICD Data.betaTilICS Data.betaTilIDC Data.betaTilIDD Data.betaTilIDS Data.betaTilISC Data.betaTilISD Data.betaTilISS]
for country in eachindex(β̃ᴷ_temp[:,1])
    β̃ᴷ[country, :, :] = reshape(β̃ᴷ_temp[country, :], NS+1, NS-1) # Capital's income share
    β̃ᴹ[country, :, :] = reshape(β̃ᴹ_temp[country, :], NS, NS) # Intermediate's income share
end

###################################################################
# Exogenous Variables
###################################################################
# Time-invariant in this numerical procedure
Â = ones(NC, NS, T)
χ̂ = ones(NC, NS, T)
d̂ = ones(NC, NC, NS, T)
L̂ = ones(NC, T)
ϕ̂ = ones(NC, T)
#ψ̂ˢ = ones(NC, T)
#ϕ̂ᵂ = ones(NC*NS, T)

# Data for Sector R
β̃ᴷᴿ = [Data.betaTilKRC Data.betaTilKRD]
β̃ᴹᴹ = [Data.betaTilIRC Data.betaTilIRD Data.betaTilIRS]

###################################################################
# Guess
###################################################################
K̂ = ones(NC, NK, T)
Ŷ = ones(NC, NK, T)
guess = ones(NC, NK, T)
guess[1:NC,1:NK,1] = K̂[1:NC,1:NK,1]
guess[1:NC,1:NK,2:T] = Ŷ[1:NC,1:NK,1:T-1]

###################################################################
# Preallocate memory
###################################################################
π = zeros(NC, NC, NS, T)
D = zeros(NC, NS+1, T)
Y = zeros(NC, NS+1, T)
X = zeros(NC, NS+1, T)
Xᶠ = zeros(NC, NS+1, T)
rK = zeros(NC, NK, T)
wL = zeros(NC, T)

###################################################################
# Initial Conditions
###################################################################
for country in 1:NC
    π[country,1:NC,1:NS,1] = Matrix(π_data[(country-1)*NC+1:country*NC, 2:4])
end
D[1:NC,1:NS+1,1] = [Data.DC1 Data.DD1 Data.DS1 Data.DR1]
Y[1:NC,1:NS+1,1] = [Data.yC1 Data.yD1 Data.yS1 Data.yR1]
X[1:NC,1:NS+1,1] = [Data.xC1 Data.xD1 Data.xS1 Data.xR1]
Xᶠ[1:NC,1:NS+1,1] = [Data.xFC1 Data.xFD1 Data.xFS1 Data.xFR1]
rK[1:NC,1:NK,1] = [Data.capCinc1 Data.capDinc1]
wL[1:NC,1] = Data.labinc1

# Pack init conditions, exogenous variables and parameters
myinit_dynamic = @with_kw (
    π₁ = π[1:NC,1:NC,1:NS,1], Y₁ = Y[1:NC,1:NS+1,1], Xᶠ₁ = Xᶠ[1:NC,1:NS+1,1], wL₁ = wL[1:NC,1], rK₁ = rK[1:NC,1:NK,1])
myexos_dynamic = @with_kw (Dᴿ = Dᴿ, L̂ = L̂, d̂ = d̂, T̂ = T̂)
myparams_dynamic = @with_kw (
    T = T, NC = NC, NS = NS, NK = NK, β̃ᴸ = β̃ᴸ, β̃ᴷ = β̃ᴷ, ψ = ψ, θ = θ, β̃ᴹ = β̃ᴹ, ρ = ρ, δ = δ, α = α)

# Config use
show(Data[1,:], allcols=:true)
