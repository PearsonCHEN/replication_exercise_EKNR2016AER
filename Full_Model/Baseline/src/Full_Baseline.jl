## Author:  Pearson
## Date: April 2020
## Julia version: 1.3.1
## Purpose: Solve EKNR(2016) Full Model - Baseline Equilibrium

# Precompile Packages
using Parameters
using CSV
using NLsolve
#using Plots

Data = CSV.read(joinpath(@__DIR__, "..", "input", "part1_21cty.csv"))
π_data = CSV.read(joinpath(@__DIR__, "..", "input", "part2_21cty.csv"))
Data = Data[Data.date.==211, 3:end]
π_data = π_data[π_data.date.==211, 5:end]
