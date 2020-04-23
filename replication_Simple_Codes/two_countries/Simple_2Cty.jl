## Author:  Pearson
## Date: April 2020
## Julia version: 1.3.1
## Purpose: Solve EKNR(2016) Simple 2 Cities Case

# Precompile Packages
using Parameters
using MATLAB
using NLsolve
using Plots

# Introduce Functions
include("Simple_2Cty_Functions.jl")

# Program Parameters
T_data = 20
T_tail = 250
T = T_data+T_tail

# Initial Conditions and Time Invariant Parameters
K_init_rel_SS=rand(2,1)
α = 0.55
βᴸ = 2/3
βᴷ = 1-βᴸ
δ = 0.1
ρ = 0.95
θ = 2
ω = [1/2, 1/2]
B = βᴸ^-βᴸ*βᴷ^-βᴷ

# Read data
fname = MatFile(joinpath(@__DIR__,"EKNR_Simple_2Cty_Shocks.mat"))
χ_data = get_variable(fname,"chi_data")
Aᴰ_data = get_variable(fname,"AD_data")
ϕ_data = get_variable(fname,"phi_data")
L_data = get_variable(fname,"L_data")
dₙᵢ_data = get_variable(fname,"dni_data")

# Hypothetical Data on Shocks(Assumed Constant During T_tail)
χ_tail = χ_data[:,T_data]*ones(1,T_tail)
Aᴰ_tail = Aᴰ_data[:,T_data]*ones(1,T_tail)
ϕ_tail = ϕ_data[:,T_data]*ones(1,T_tail)
L_tail = L_data[:,T_data]*ones(1,T_tail)
dₙᵢ = repeat(dₙᵢ_data[:,:,T_data],1,1,T)
χ = [χ_data χ_tail]
Aᴰ = [Aᴰ_data Aᴰ_tail]
ϕ = [ϕ_data ϕ_tail]
L = [L_data L_tail]
dₙᵢ[:,:,1:T_data] = dₙᵢ_data

# Hypothetical Data on Shocks in Change Form
χ̂ = [χ[:,2:end]./χ[:,1:end-1] ones(eltype(χ),size(χ[:,end]))]
Âᴰ = [Aᴰ[:,2:end]./Aᴰ[:,1:end-1] ones(eltype(Aᴰ),size(Aᴰ[:,end]))]
ϕ̂ = [ϕ[:,2:end]./ϕ[:,1:end-1] ones(eltype(ϕ),size(ϕ[:,end]))]
L̂ = [L[:,2:end]./L[:,1:end-1] ones(eltype(L),size(L[:,end]))]
d̂ₙᵢ = ones(eltype(dₙᵢ), size(dₙᵢ)); d̂ₙᵢ[:,:,1:end-1] = dₙᵢ[:,:,2:end]./dₙᵢ[:,:,1:end-1]

###########################################################################################
#                            Solve for the steady state
###########################################################################################
# Guess and Other Exogenous Variables
Yˢₛₛ = ω.*ϕ[:,T]
guess_ss = [1.5*Yˢₛₛ 200*Yˢₛₛ]

# Wrap All Exogenous Variables and Parameters
myexos_ss = @with_kw (χ = χ[:,T], Aᴰ = Aᴰ[:,T], L = L[:,T], dₙᵢ = dₙᵢ[:,:,T], Yˢₛₛ = Yˢₛₛ)
myparams_ss = @with_kw (α = α, βᴸ = βᴸ, βᴷ = βᴷ, δ = δ, ρ = ρ, θ = θ, B = B)
exos_ss = myexos_ss()
params_ss = myparams_ss()

# Solve for Steady State
println("\n\nStart to solve the steady state.")
println("Run time and memory cost:")
@time results_ss = nlsolve(
    (res_ss, guess_ss) -> Fun_2Cty_SS!(res_ss, guess_ss, exos_ss, params_ss),
    guess_ss,
    ftol=1e-6,
    method=:newton,
    autodiff=:forward,
    show_trace=false,
)

# Check Convergence
converged(results_ss) || error("Failed to converge in $(results_ss.iterations) iterations.")
println("Successfully solved steady state.\n")

# Catch Solutions for Steady State
res_ss = similar(results_ss.zero)
res_ss, Yₛₛ, Kₛₛ = Fun_2Cty_SS!(res_ss, results_ss.zero, exos_ss, params_ss)

###########################################################################################
#                         Solve for Solution to the Problem in Levels
###########################################################################################
# Initial Conditions, Guess and Other Exogenous Variables
K₁ = K_init_rel_SS.*Kₛₛ
Yˢ = ω.*ϕ
guess_lev = repeat(Yₛₛ,1,T)

# Wrap All Exogenous Variables and Parameters
myinit_lev = @with_kw (K₁ = K₁,)
myexos_lev = @with_kw (χ = χ, Aᴰ = Aᴰ, L = L, dₙᵢ = dₙᵢ, Kₛₛ = Kₛₛ, Yˢ = Yˢ)
myparams_lev = @with_kw (α = α, βᴷ = βᴷ, βᴸ = βᴸ, δ = δ, ρ = ρ, θ = θ, T = T)
init_lev = myinit_lev()
exos_lev = myexos_lev()
params_lev = myparams_lev()

# Solve the Problem in Levels
println("Start to solve the problem in levels.")
println("Run time and memory cost:")
@time results_lev =
    try
        results_lev = nlsolve(
            (res_lev, guess_lev) -> Fun_2Cty_Levels!(
                res_lev, guess_lev, init_lev, exos_lev, params_lev),
            guess_lev,
            ftol=1e-6,
            method=:newton,
            autodiff=:forward,
            show_trace=false,
        )
    catch err
        if isa(err, DomainError)
            error("Failed to solve the problem in levels, please try again.")
        end
    end

# Check Convergence
converged(results_lev) || error("Failed to converge in $(results_lev.iterations) iterations.")
println("Successfully solved the problem in levels.\n")

# Catch Solutions
res_lev = similar(results_lev.zero)
res_lev, K_lev, Y_lev, Yᴰ_lev, Xᴰ_lev, πₙᵢ_lev = Fun_2Cty_Levels!(
    res_lev, results_lev.zero, init_lev, exos_lev, params_lev)

###########################################################################################
#                         Solve for Solution to the Problem in Changes
###########################################################################################
# Guess
guess_hat = ones(2,T)

# Wrap All Exogenous Variables and Parameters
myinit_hat = @with_kw (Y₁ = Y_lev[:,1], Yᴰ₁ = Yᴰ_lev[:,1], πₙᵢ₁ = πₙᵢ_lev[:,:,1])
myexos_hat = @with_kw (χ̂ = χ̂, Âᴰ = Âᴰ, L̂ = L̂, d̂ₙᵢ = d̂ₙᵢ, ϕ̂ = ϕ̂)
myparams_hat = @with_kw (α = α, βᴷ = βᴷ, βᴸ = βᴸ, δ = δ, ρ = ρ, θ = θ, T = T)
init_hat = myinit_hat()
exos_hat = myexos_hat()
params_hat = myparams_hat()

# Solve the Problem in Changes
println("Start to solve the problem in changes.")
println("Run time and memory cost:")
@time results_hat =
    try
        results_hat = nlsolve(
            (res_hat, guess_hat) -> Fun_2Cty_Changes!(
                res_hat, guess_hat, init_hat, exos_hat, params_hat),
            guess_hat,
            ftol=1e-6,
            method=:newton,
            autodiff=:forward,
            show_trace=false,
        )
    catch err
        if isa(err, DomainError)
            error("Failed to solve the problem in changes, please try again.")
        end
    end

# Check Convergence
converged(results_hat) || error("Failed to converge in $(results_hat.iterations) iterations.")
println("Successfully solved the problem in changes.\n")

# Catch Solutions
res_hat = similar(results_hat.zero)
res_hat, K̂, Ŷ, Ŷᴰ, X̂ᴰ, π̂ₙᵢ = Fun_2Cty_Changes!(
    res_hat, results_hat.zero, init_hat, exos_hat, params_hat)
K_hat = similar(K_lev)
Y_hat = similar(Y_lev)
K_hat[:,1] = K₁
Y_hat[:,1] = Y_lev[:,1]
for t = 2:T
    K_hat[:,t] = K_hat[:,t-1].*K̂[:,t-1]
    Y_hat[:,t] = Y_hat[:,t-1].*Ŷ[:,t-1]
end

###########################################################################################
#                                      Plot Figures
###########################################################################################
println("Plotting the figure.")
l = @layout [a b; c d]
p1 = Plots.plot(K_lev[1,1:80], line=(:solid, 2), label="Level Solution");
    Plots.plot!(K_hat[1,1:80], line=(:dash, 2), label="Change Solution");
p2 = Plots.plot(K_lev[2,1:80], line=(:solid, 2), label="Level Solution");
    Plots.plot!(K_hat[2,1:80], line=(:dash, 2), label="Change Solution");
p3 = Plots.plot(Y_lev[1,1:80], line=(:solid, 2), label="Level Solution");
    Plots.plot!(Y_hat[1,1:80], line=(:dash, 2), label="Change Solution");
p4 = Plots.plot(Y_lev[2,1:80], line=(:solid, 2), label="Level Solution");
    Plots.plot!(Y_hat[2,1:80], line=(:dash, 2), label="Change Solution");
figures_all = plot(p1, p2, p3, p4, layout = l)
display(figures_all)
savefig(figures_all, "EKNR_2Cty.pdf")
println("The figure is saved.")
