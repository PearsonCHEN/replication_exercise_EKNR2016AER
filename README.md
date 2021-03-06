# Replication Exercise for Eaton et al.(2016, AER)
This project is a replication exercise for the seminal work about "hat algebra"(Eaton et al., 2016).

# Purpose
The authors of the [original paper](https://www.aeaweb.org/articles?id=10.1257/aer.20101557) have already offered replication files in AEA(American Economic Association). However, their analysis is based on Stata and MATLAB. The primary goal of this replication exercise is to learn both the programming language "Julia" and "hat algebra" method in solving economic dynamic transition.

# Julia
Julia is an outstanding young programming language for scientific computation. Julia offers excellent features in **type stability** and **multiple-dispatch**(See [JuliaDocumentation](https://docs.julialang.org/en/v1/manual/documentation/index.html) and [QuantEcon](https://julia.quantecon.org/)). Although Julia is quite similar to other script languages(MATLAB, Python, R, etc.), the code has to be well-designed to achieve peak performance. This replication exercise is a simple investigation.

# Hat algebra
Hat algebra transforms the traditional dynamic economic system to that of representation of flow variables and "changes". This method helps shy away from unclear definition of exogenous variables, or fundamentals(Caliendo et al., 2019) like capital, trade cost and productivity shocks etc.

# Reference
Eaton, Jonathan, Samuel Kortum, Brent Neiman, and John Romalis, "Trade and the Global Recession." American Economic Review, 106 (11): 3401-38, 2016. DOI: 10.1257/aer.20101557  
Caliendo, Lorenzo, Maximiliano Dvorkin, and Fernando Parro, "Trade and labor market dynamics: General equilibrium analysis of the china trade shock."  Econometrica, 87 (3), 741-835, 2019. DOI: 10.3982/ECTA13758
