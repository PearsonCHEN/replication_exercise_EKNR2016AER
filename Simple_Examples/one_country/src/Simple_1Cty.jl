## Author:  Pearson
## Date: April 2020
## Julia version: 1.3.1
## Purpose: Solve EKNR(2016) Simple 1 City Case

# Precompile Packages
using Parameters
using MATLAB
using NLsolve
using Plots

# Introduce Functions
include("Simple_1Cty_Functions.jl")

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
fname = MatFile(joinpath(@__DIR__, "..", "input", "EKNR_Simple_1Cty_Shocks.mat"))
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

# Hypothetical Data on Shocks in Change Form
χ̂ = [χ[2:T]./χ[1:T-1]; 1]
Âᴰ = [Aᴰ[2:T]./Aᴰ[1:T-1]; 1]
L̂ = [L[2:T]./L[1:T-1]; 1]

# Steady State Values (FROM ANALYTICAL SOLUTION)
Kₛₛ = L[T]*βᴷ/βᴸ*(ρ*χ[T]*Aᴰ[T]/(1-ρ*(1-δ)))^(1/βᴸ)
Yₛₛ = (1-ρ*(1-δ))/(1-ρ*(1-δ)-δ*ρ*βᴷ)

###########################################################################################
#                         Solve for Solution to the Problem in Levels
###########################################################################################
# Initial Conditions and Guess
K₁ = Kₛₛ*K_init_rel_SS
guess_lev = ones(eltype(Yₛₛ),T,1).*Yₛₛ

# Wrap All Exogenous Variables and Parameters
myinit_lev = @with_kw (K₁ = K₁,)
myexos_lev = @with_kw (χ = χ, Aᴰ = Aᴰ, L = L, K₁ = K₁, Kₛₛ = Kₛₛ)
myparams_lev = @with_kw (βᴷ = βᴷ, βᴸ = βᴸ, δ = δ, ρ = ρ, B = B, T = T)
init_lev = myinit_lev()
exos_lev = myexos_lev()
params_lev = myparams_lev()

# Solve for the Problem in levels
println("\n\nStart to solve the problem in levels.")
println("Run time and memory cost:")
@time results_lev =
    try
        results_lev = nlsolve(
        (res_lev, guess_lev) -> Fun_1Cty_Levels!(res_lev,
            guess_lev, init_lev, exos_lev, params_lev),
        guess_lev,
        ftol=1e-6,
        autodiff=:forward,
        method=:newton,
        show_trace=false,
    )
    catch err
        if isa(err, DomainError)
            error("Failed to solve the problem in levels, please try again.")
        end
    end

# Check Convergence
converged(results_lev) || error("Failed to converge in $(results_lev.iterations) iterations")
println("Successfully solved problem in levels.\n")

# Catch Solutions
res_lev = similar(results_lev.zero)
res_lev, K_lev, Y_lev = Fun_1Cty_Levels!(
    res_lev, results_lev.zero, init_lev, exos_lev, params_lev)

###########################################################################################
#                         Solve for Solution to the Problem in Changes
###########################################################################################
# Guess
guess_hat = ones(T,1)

# Wrap All Exogenous Variables and Parameters
myinit_hat = @with_kw (Y₁ = Y_lev[1],)
myexos_hat = @with_kw (χ̂ = χ̂, Âᴰ = Âᴰ, L̂ = L̂)
myparams_hat = @with_kw (βᴷ = βᴷ, βᴸ = βᴸ, δ = δ, ρ = ρ, B = B, T = T)
init_hat = myinit_hat()
exos_hat = myexos_hat()
params_hat = myparams_hat()

# Solve for Solution to Problem in Changes
println("Start to solve the problem in changes.")
println("Run time and memory cost:")
@time results_hat =
    try
        results_hat = nlsolve(
            (res_hat, guess_hat) -> Fun_1Cty_Changes!(
                res_hat, guess_hat, init_hat, exos_hat, params_hat),
        guess_hat,
        ftol=1e-6,
        autodiff=:forward,
        method=:newton,
        show_trace=false,
    )
    catch err
        if isa(err, DomainError)
            error("Failed to solve the problem in changes, please try again.")
        end
    end

# Check Convergence
converged(results_hat) || error("Failed to converge in $(results_hat.iterations) iterations")
println("Successfully solved problem in changes.\n")

# Catch Solutions
res_hat = similar(results_hat.zero)
res_hat, K̂, Ŷ = Fun_1Cty_Changes!(
    res_hat, results_hat.zero, init_hat, exos_hat, params_hat)
K_hat = similar(K̂)
Y_hat = similar(Ŷ)
K_hat[1] = K₁
Y_hat[1] = Y_lev[1]
for tt in 2:T
    K_hat[tt] = K_hat[tt-1]*K̂[tt-1]
    Y_hat[tt] = Y_hat[tt-1]*Ŷ[tt-1]
end

###########################################################################################
#                                      Plot Figures
###########################################################################################
println("Plotting the figure.")
l = @layout [a; b]
p1 = Plots.plot(K_lev[1:75], line=(:solid, 2), label="Level Solution");
    Plots.plot!(K_hat[1:75], line=(:dash, 2), label="Change Solution")
p2 = Plots.plot(Y_lev[1:75], line=(:solid, 2), label="Level Solution");
    Plots.plot!(Y_hat[1:75], line=(:dash, 2), label="Change Solution")
figures_all = plot(p1, p2, layout = l)
display(figures_all)
savefig(figures_all, joinpath(@__DIR__, "..", "output", "EKNR_1Cty.pdf"))
println("The figure is saved.")
