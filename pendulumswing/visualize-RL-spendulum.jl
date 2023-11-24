using MAT
using JLD
using LinearAlgebra
using Distributions
using Plots; 
default(grid=false, label="", linewidth=5,margin=15Plots.pt)

results_EFE = load("pendulumswing/results/EFE.jld")
results_GPRL = matread("pendulumswing/results/GPRL.mat")

T      = results_EFE["T"]
tsteps = results_EFE["tsteps"]
ulims  = results_EFE["sys_ulims"]
goal_m = results_EFE["goals_m"]
goal_v = results_EFE["goals_v"]
y_GPRL = results_GPRL["y_GPRL"]
u_GPRL = results_GPRL["u_GPRL"]
N = size(y_GPRL,1)

"Plot trajectories"

p1 = plot(ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), guidefontsize=12, legendfontsize=10, tickfontsize=10)
plot!(tsteps, y_GPRL[:,1], color="green", label="trial 1")
hline!([goal_m[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="", xlabel="time [s]", ylabel="angle (ϑ)", guidefontsize=12, legendfontsize=10, tickfontsize=10)

p2 = plot(ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), guidefontsize=12, legendfontsize=10, tickfontsize=10)
plot!(tsteps, y_GPRL[:,2], yaxis=false, color="green", label="trial 2")
hline!([goal_m[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="", xlabel="time [s]", ylabel="", guidefontsize=12, legendfontsize=10, tickfontsize=10)

p3 = plot(ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), guidefontsize=12, legendfontsize=10, tickfontsize=10)
plot!(tsteps, y_GPRL[:,3], yaxis=false, color="green", label="trial 3")
hline!([goal_m[1]], color="black", linestyle=:dash, alpha=0.5, label="", xlabel="time [s]", ylabel="", guidefontsize=12, legendfontsize=10, tickfontsize=10)

p4 = plot(ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), guidefontsize=12, legendfontsize=10, tickfontsize=10)
plot!(tsteps, y_GPRL[:,4], yaxis=false, color="green", label="trial 4")
hline!([goal_m[1]], color="black", linestyle=:dash, alpha=0.5, label="", xlabel="time [s]", ylabel="", guidefontsize=12, legendfontsize=10, tickfontsize=10)

plot(p1,p2,p3,p4, layout=(1,4), size=(900,200), right_margin=[-20Plots.pt -20Plots.pt -20Plots.pt 10Plots.pt])
savefig("pendulumswing/figures/GPRL-trajectories.png")
