using JLD
using LinearAlgebra
using Distributions
using Plots; 
default(grid=false, label="", linewidth=3,margin=20Plots.pt)

results_EFE = load("pendulumswing/results/EFE.jld")
results_QCR = load("pendulumswing/results/QCR.jld")

T      = results_EFE["T"]
tsteps = results_EFE["tsteps"]
ulims  = results_EFE["sys_ulims"]
goal_m = results_EFE["goals_m"]
goal_v = results_EFE["goals_v"]
y_EFE  = results_EFE["y_EFE"]
y_QCR  = results_QCR["y_QCR"]
u_EFE  = results_EFE["u_EFE"]
u_QCR  = results_QCR["u_QCR"]

"Plot trajectories"

p1 = plot(xlabel="time [s]", ylabel="angle (ϑ)", guidefontsize=12, legendfontsize=12, tickfontsize=12)
hline!([goal_m[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="")
plot!(tsteps, y_EFE[1:end], linewidth=5, color="red", label="EFE")

p2 = plot(xlabel="time [s]", ylabel="angle (ϑ)", guidefontsize=12, legendfontsize=12, tickfontsize=12)
hline!([goal_m[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="")
plot!(tsteps, y_QCR[1:end], linewidth=5, color="blue", label="QCR")

plot(p2,p1, layout=(1,2), size=(900,250), margin=20Plots.pt)

savefig("pendulumswing/figures/angle-trajectories.png")

p1 = plot(xlabel="time [s]", ylabel="angle (ϑ)", guidefontsize=12, legendfontsize=10, tickfontsize=10)
hline!([goal_m[1]], color="black", linestyle=:dash, alpha=0.5, ylims=(-2, 8), yticks=([0, 3.1415, 6.2830],[0, "π",  "2π"]), label="")
plot!(tsteps, y_EFE[1:end], linewidth=5, color="red", label="EFE")
plot!(tsteps, y_QCR[1:end], linewidth=5, color="blue", label="QCR")
plot!(size=(900,250), margin=20Plots.pt)

"Plot controls"

p1 = hline([ulims[1]], color="black", linestyle=:dash, alpha=0.5, xlabel="time [s]")
hline!([ulims[2]], color="black", linestyle=:dash, alpha=0.5)
plot!(tsteps, u_EFE[1:end-T], linewidth=3, color="red", label="EFE", ylims=ulims.*1.3, xlabel="", ylabel="control (u)", guidefontsize=12, legendfontsize=10, tickfontsize=10)

p2 = hline([ulims[1]], color="black", linestyle=:dash, alpha=0.5, xlabel="time [s]")
hline!([ulims[2]], color="black", linestyle=:dash, alpha=0.5)
plot!(tsteps, u_QCR[1:end-T], linewidth=3, color="blue", label="QCR", ylims=ulims.*1.3, xlabel="", ylabel="control (u)", guidefontsize=12, legendfontsize=10, tickfontsize=10)

plot(p1,p2, layout=grid(1,2), size=(900,200), margin=20Plots.pt)

savefig("pendulumswing/figures/control-trajectories.png")

p1 = plot(ylabel="control (u)", xlabel="time [s]", guidefontsize=12, legendfontsize=10, tickfontsize=10)
hline!([ulims[1]], color="black", linestyle=:dash, alpha=0.5)
hline!([ulims[2]], color="black", linestyle=:dash, alpha=0.5)
plot!(tsteps, u_EFE[1:end-T], linewidth=3, color="red", label="EFE", ylims=ulims.*1.3)
plot!(tsteps, u_QCR[1:end-T], linewidth=3, color="blue", label="QCR", ylims=ulims.*1.3)
plot!(size=(900,250), margin=20Plots.pt)
