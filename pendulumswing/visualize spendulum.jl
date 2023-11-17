using MAT
using JLD
using LinearAlgebra
using Distributions
using Plots; 
default(grid=false, label="", linewidth=3,margin=20Plots.pt)

results_EFE = load("pendulumswing/results/EFE.jld")
results_MSE = load("pendulumswing/results/MSE.jld")
results_NMPC = matread("pendulumswing/results/NMPC.mat")

T      = results_EFE["T"]
tsteps = results_EFE["tsteps"]
ulims  = results_EFE["sys_ulims"]
goal_m = results_EFE["goals_m"]
goal_v = results_EFE["goals_v"]
y_EFE  = results_EFE["y_EFE"]
y_MSE  = results_MSE[ "y_MSE"]
y_NMPC = results_NMPC["y_NMPC"]
u_EFE  = results_EFE["u_EFE"]
u_MSE  = results_MSE[ "u_MSE"]
u_NMPC = results_NMPC["u_NMPC"]

"Plot trajectories"

p1 = hline([goals[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="", xlabel="", ylabel="angle (ϑ)", guidefontsize=12, legendfontsize=10, tickfontsize=10)
plot!(tsteps, y_EFE[1:end], linewidth=5, color=palette(:tab10)[1], label="EFE")

p2 = hline([goals[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="", xlabel="", ylabel="angle (ϑ)", guidefontsize=12, legendfontsize=10, tickfontsize=10)
plot!(tsteps, y_MSE[1:end], linewidth=5, color=palette(:tab10)[2], label="MSE")

p3 = hline([goals[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="", xlabel="time [s]", ylabel="angle (ϑ)", guidefontsize=12, legendfontsize=10, tickfontsize=10)
plot!(tsteps, y_NMPC[2:end], linewidth=5, color=palette(:tab10)[3], label="known dynamics")

plot(p1,p2,p3, layout=(3,1), size=(600,500), margin=10Plots.pt)

savefig("pendulumswing/figures/angle-trajectories.png")

"Plot controls"

p1 = hline([ulims[1]], color="black", linestyle=:dash, alpha=0.5)
hline!([ulims[2]], color="black", linestyle=:dash, alpha=0.5)
plot!(tsteps, u_EFE[1:end-T], linewidth=3, color=palette(:tab10)[1], label="EFE", ylims=sys_ulims.*1.3, xlabel="", ylabel="control (u)", guidefontsize=12, legendfontsize=10, tickfontsize=10)

p2 = hline([ulims[1]], color="black", linestyle=:dash, alpha=0.5)
hline!([ulims[2]], color="black", linestyle=:dash, alpha=0.5)
plot!(tsteps, u_MSE[1:end-T], linewidth=3, color=palette(:tab10)[2], label="MSE", ylims=sys_ulims.*1.3, xlabel="", ylabel="control (u)", guidefontsize=12, legendfontsize=10, tickfontsize=10)

p3 = hline([ulims[1]], color="black", linestyle=:dash, alpha=0.5)
hline!([ulims[2]], color="black", linestyle=:dash, alpha=0.5)
plot!(tsteps, u_NMPC[2:end], linewidth=3, color=palette(:tab10)[3], label="known dynamics", ylims=sys_ulims.*1.3, xlabel="time [s]", ylabel="control (u)", guidefontsize=12, legendfontsize=10, tickfontsize=10)

plot(p1,p2,p3, layout=grid(3,1), size=(600,500), margin=10Plots.pt)

savefig("pendulumswing/figures/control-trajectories.png")