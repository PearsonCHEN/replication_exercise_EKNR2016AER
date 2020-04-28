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
K̂ = ones(NC, NS, T)
Ŷ = ones(NC, NS, T)
guess = ones(NC, NS, T); guess[:,:,1] = K̂[:,:,1]; guess[:,:,2:T] = Ŷ[:,:,1:T-1]

###################################################################
# Initial Conditions
###################################################################
π = ones(NC, NC, NS)
for country in 1:1:NC
    π[country, :, :] = Matrix(π_data[(country-1)*NC+1:country*NC, 2:4])
end
D = [Data.DC1 Data.DD1 Data.DS1 Data.DR1]
Y₀ = [Data.yC1 Data.yD1 Data.yS1 Data.yR1]
X₀ = [Data.xC1 Data.xD1 Data.xS1 Data.xR1]
Xᶠ = [Data.xFC1 Data.xFD1 Data.xFS1 Data.xFR1]
rK = [Data.capCinc1 Data.capDinc1]
wL = Data.labinc1

# Config use
show(Data[1,:], allcols=:true)
