using Revise
using JLD
using Statistics
using LaTeXStrings
using Plots
# pyplot()
default(label="", linewidth=5, margin=12Plots.pt)

includet("../mv_normal_gamma.jl"); using .mv_normal_gamma

num_exps = 100
N = 100

id_rnd = zeros(num_exps, N)
id_MSE = zeros(num_exps, N)
id_EFE = zeros(num_exps, N)
id_sin = zeros(num_exps, N)

for nn in 1:num_exps

    sys_settings = load("./experiments/data/system-$nn.jld")
    sys_theta = sys_settings["sys_theta"]
    sys_sd = sys_settings["sys_sd"]

    DD = load("./experiments/results/learningrate-rnd-$nn.jld")
    id_rnd[nn,:] = [mv_normal_gamma.logpdf(MvNormalGamma(DD["mu"][k], DD["Lambda"][k], DD["alpha"][k], DD["beta"][k]), sys_theta, inv(sys_sd^2)) for k in 1:N]
    
    DD = load("./experiments/results/learningrate-MSE-$nn.jld")
    id_MSE[nn,:] = [mv_normal_gamma.logpdf(MvNormalGamma(DD["mu"][k], DD["Lambda"][k], DD["alpha"][k], DD["beta"][k]), sys_theta, inv(sys_sd^2)) for k in 1:N]

    DD = load("./experiments/results/learningrate-EFE-$nn.jld")
    id_EFE[nn,:] = [mv_normal_gamma.logpdf(MvNormalGamma(DD["mu"][k], DD["Lambda"][k], DD["alpha"][k], DD["beta"][k]), sys_theta, inv(sys_sd^2)) for k in 1:N]
    
    DD = load("./experiments/results/learningrate-sin-$nn.jld")
    id_sin[nn,:] = [mv_normal_gamma.logpdf(MvNormalGamma(DD["mu"][k], DD["Lambda"][k], DD["alpha"][k], DD["beta"][k]), sys_theta, inv(sys_sd^2)) for k in 1:N]
end

rate_rnd = mean(id_rnd, dims=1)'
rate_sin = mean(id_sin, dims=1)'
rate_MSE = mean(id_MSE, dims=1)'
rate_EFE = mean(id_EFE, dims=1)'

index = 1:100
mint = round(minimum(rate_rnd[index])-14)
maxt = round(maximum(rate_rnd[index]))

plot(xlabel="time (#steps)", 
    grid=true, 
    guidefontsize=12, 
    tickfontsize=12, 
    legendfontsize=12, 
    ylabel="log-probability", 
    size=(600,300), 
    ylims=(mint, maxt),
    yticks=[round(mint), round(mean([mint,maxt])), round(maxt)],
    # yscale=:log10,
    legend=:bottomright)

# plot!(rate_rnd, linewidth=5, label="rnd")
plot!(rate_sin[index], linewidth=4, label="sine")
plot!(rate_MSE[index], linewidth=4, label="MSE")
plot!(rate_EFE[index], linewidth=4, label="EFE")

savefig("./experiments/figures/compare-learningrates.png")