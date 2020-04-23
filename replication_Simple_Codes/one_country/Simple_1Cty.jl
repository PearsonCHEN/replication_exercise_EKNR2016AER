## Author:  Pearson
## Julia version: 1.3.1
## Purpose: Solve EKNR(2016) Simple 1 City Case
using Pkg
using Parameters

# Version Control
include("Simple_1Cty_Functions.jl")

for package in ["MATLAB", "NLsolve"]
  if Pkgcheck()[package]==nothing
    Pkg.add(package)
  end
end

#Pin packages to version used in the paper
if Pkgcheck()["MATLAB"]!=v"0.7.3"
    Pkg.pin(PackageSpec(name="MATLAB", version="0.7.3"))
end

if Pkgcheck()["NLsolve"]!=v"1.0.1"
    Pkg.pin(PackageSpec(name="NLsolve", version="4.3.0"))
end

#Compile packages
using MATLAB
using NLsolve

# Program Parameters
T_data = 20
T_tail = 500
T = T_data+T_tail

# Initial Conditions and Time Invariant Parameters
K_init_rel_SS = 1/2
βᴸ = 2/3
βᴷ = 1-βᴸ
δ = 0.1
ρ = 0.95
B = βᴸ^-βᴸ*βᴷ^-βᴷ

# Read data
fname = MatFile(joinpath(@__DIR__,"EKNR_Simple_1Cty_Shocks.mat"))
χ_data = get_variable(fname,"chi_data")
Aᴰ_data = get_variable(fname,"AD_data")
L_data = get_variable(fname,"L_data")

# Hypothetical Data on Shocks(Assumed Constant During T_tail)
χ_tail = χ_data[T_data]*ones(T_tail,1)
Aᴰ_tail = Aᴰ_data[T_data]*ones(T_tail,1)
L_tail = L_data[T_data]*ones(T_tail,1)
χ = [χ_data'; χ_tail]
Aᴰ = [Aᴰ_data'; Aᴰ_tail]
L = [L_data'; L_tail]

# Steady State Values (FROM ANALYTICAL SOLUTION)
Kₛₛ = L[T]*βᴷ/βᴸ*(ρ*χ[T]*Aᴰ[T]/(1-ρ*(1-δ)))^(1/βᴸ)
Yₛₛ = (1-ρ*(1-δ))/(1-ρ*(1-δ)-δ*ρ*βᴷ)

#=
    Solve for Solution to Problem in Levels
=#

# Initial Conditions and Guess
K₁ = Kₛₛ*K_init_rel_SS
lev_guess = ones(T,1)*Yₛₛ

# Wrap All Exogenous Variables and Parameters
myexos_1Cty_lev = @with_kw (χ = χ, Aᴰ = Aᴰ, L = L, K₁ = K₁, Kₛₛ = Kₛₛ)
myparams_1Cty_lev = @with_kw (βᴷ = βᴷ, βᴸ = βᴸ, δ = δ, ρ = ρ, B = B, T = T)
exos_1Cty_lev = myexos_1Cty_lev()
params_1Cty_lev = myparams_1Cty_lev()

# Solve for the Problem in levels
@time results_lev = nlsolve(
    (err, lev_guess) -> Fun_1Cty_Levels!(err, lev_guess, exos_1Cty_lev, params_1Cty_lev),
    lev_guess,
    ftol=1e-6,
    autodiff=:forward,
    method=:newton,
    show_trace=true,
)

# Check Convergence
converged(results_lev) || error("Failed to converge in $(results_lev.iterations) iterations")
println("Successfully solved problem in levels.")

# Catch Solutions
Y_levsolution = results_lev.zero
err = similar(Y_levsolution)
residual_lev, K_levsolution = Fun_1Cty_Levels!(
    err,
    results_lev.zero,
    exos_1Cty_lev,
    params_1Cty_lev,
)

# Hypothetical Data on Shocks in Change Form
χ̂ = [χ[2:T]./χ[1:T-1]; 1]
Âᴰ = [Aᴰ[2:T]./Aᴰ[1:T-1]; 1]
L̂ = [L[2:T]./L[1:T-1]; 1]

#=
    Solve for Solution to Problem in changes
=#

# Wrap All Exogenous Variables and Parameters
myexos_1Cty_hat = @with_kw (χ̂ = χ̂, Âᴰ = Âᴰ, L̂ = L̂, Y₁ = Y_levsolution[1])
myparams_1Cty_hat = @with_kw (βᴷ = βᴷ, βᴸ = βᴸ, δ = δ, ρ = ρ, B = B, T = T)
exos_1Cty_hat = myexos_1Cty_hat()
params_1Cty_hat = myparams_1Cty_hat()

# Solve for Solution to Problem in Changes
hat_guess = ones(T,1)
@time results_hat = nlsolve(
    (err, hat_guess) -> Fun_1Cty_Changes!(err, hat_guess, exos_1Cty_hat, params_1Cty_hat),
    hat_guess,
    ftol=1e-6,
    autodiff=:forward,
    method=:newton,
    show_trace=true,
)

# Check Convergence
converged(results_hat) || error("Failed to converge in $(results_hat.iterations) iterations")
println("Successfully solved problem in changes.")

# Catch Solutions
err = similar(results_hat.zero)
residual_hat, hat_solutions = Fun_1Cty_Changes!(
    err,
    results_hat.zero,
    exos_1Cty_hat,
    params_1Cty_hat,
)
K_hatsolution = hat_solutions[1:T,1]
Y_hatsolution = hat_solutions[1:T,2]

#=
    Plot Figures
=#
K_hatsolution_lev = similar(K_hatsolution)
Y_hatsolution_lev = similar(Y_hatsolution)
K_hatsolution_lev[1] = K₁
Y_hatsolution_lev[1] = Y_levsolution[1]

for tt in 2:T
    K_hatsolution_lev[tt] = K_hatsolution_lev[tt-1]*K_hatsolution[tt-1]
    Y_hatsolution_lev[tt] = Y_hatsolution_lev[tt-1]*Y_hatsolution[tt-1]
end

using Plots
l = @layout [a; b]
p1 = Plots.plot(K_levsolution[1:75], line=(:solid, 2));
    Plots.plot!(K_hatsolution_lev[1:75], line=(:dash, 2))
p2 = Plots.plot(Y_levsolution[1:75], line=(:solid, 2));
    Plots.plot!(Y_hatsolution_lev[1:75], line=(:dash, 2))
plot(p1, p2, layout = l)
