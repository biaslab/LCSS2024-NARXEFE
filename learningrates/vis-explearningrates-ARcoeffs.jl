using Revise
using JLD
using LinearAlgebra
using Statistics
using Distributions
using LaTeXStrings
using Plots
default(label="", linewidth=5, margin=15Plots.pt)

# includet("../mv_normal_gamma.jl"); using .mv_normal_gamma

num_exps = 100
N = 100

ll_rnd = zeros(num_exps, N)
ll_MSE = zeros(num_exps, N)
ll_EFE = zeros(num_exps, N)
ll_sin = zeros(num_exps, N)

n2_rnd = zeros(num_exps, N)
n2_MSE = zeros(num_exps, N)
n2_EFE = zeros(num_exps, N)
n2_sin = zeros(num_exps, N)

for nn in 1:num_exps

    sys_settings = load("./experiments/data/system-$nn.jld")
    sys_theta = sys_settings["sys_theta"]

    DD = load("./experiments/results/learningrate-rnd-$nn.jld")
    ll_rnd[nn,:] = [logpdf(MvNormal(DD["mu"][k], inv(Hermitian(DD["Lambda"][k]))), sys_theta) for k in 1:N]
    n2_rnd[nn,:] = [norm(sys_theta - DD["mu"][k]) for k in 1:N]
    
    DD = load("./experiments/results/learningrate-MSE-$nn.jld")
    ll_MSE[nn,:] = [logpdf(MvNormal(DD["mu"][k], inv(Hermitian(DD["Lambda"][k]))), sys_theta) for k in 1:N]
    n2_MSE[nn,:] = [norm(sys_theta - DD["mu"][k], 2) for k in 1:N]

    DD = load("./experiments/results/learningrate-EFE-$nn.jld")
    ll_EFE[nn,:] = [logpdf(MvNormal(DD["mu"][k], inv(Hermitian(DD["Lambda"][k]))), sys_theta) for k in 1:N]
    n2_EFE[nn,:] = [norm(sys_theta - DD["mu"][k], 2) for k in 1:N]
    
    DD = load("./experiments/results/learningrate-sin-$nn.jld")
    ll_sin[nn,:] = [logpdf(MvNormal(DD["mu"][k], inv(Hermitian(DD["Lambda"][k]))), sys_theta) for k in 1:N]
    n2_sin[nn,:] = [norm(sys_theta - DD["mu"][k], 2) for k in 1:N]
end

mll_rnd = mean(ll_rnd, dims=1)'
mll_sin = mean(ll_sin, dims=1)'
mll_MSE = mean(ll_MSE, dims=1)'
mll_EFE = mean(ll_EFE, dims=1)'
sll_rnd = std( ll_rnd, dims=1)' ./ num_exps
sll_sin = std( ll_sin, dims=1)' ./ num_exps
sll_MSE = std( ll_MSE, dims=1)' ./ num_exps
sll_EFE = std( ll_EFE, dims=1)' ./ num_exps

index = 1:100

plot(xlabel="time (#steps)", 
    grid=true, 
    guidefontsize=12, 
    tickfontsize=12, 
    legendfontsize=12, 
    ylabel=L"p(\theta^{*} | \mathcal{D}_k)", 
    size=(600,200), 
    xscale=:log10,
    legend=:topleft)

# plot!(rate_rnd, linewidth=5, label="rnd")
plot!(mll_sin[index], ribbon=sll_sin[index], linewidth=4, label="sine")
plot!(mll_MSE[index], ribbon=sll_MSE[index], linewidth=4, label="MSE")
plot!(mll_EFE[index], ribbon=sll_EFE[index], linewidth=4, label="EFE")

savefig("./experiments/figures/probpost-idARcoeffs.png")
savefig("./experiments/figures/probpost-idARcoeffs.pdf")

mn2_rnd = mean(n2_rnd, dims=1)'
mn2_sin = mean(n2_sin, dims=1)'
mn2_MSE = mean(n2_MSE, dims=1)'
mn2_EFE = mean(n2_EFE, dims=1)'
sn2_rnd = std( n2_rnd, dims=1)' ./ num_exps
sn2_sin = std( n2_sin, dims=1)' ./ num_exps
sn2_MSE = std( n2_MSE, dims=1)' ./ num_exps
sn2_EFE = std( n2_EFE, dims=1)' ./ num_exps

plot(xlabel="time (#steps)", 
    grid=true, 
    guidefontsize=12, 
    tickfontsize=12, 
    legendfontsize=12, 
    ylabel=L"|| \theta^{*} - \mu_k ||_2", 
    size=(600,200), 
    xscale=:log10,
    legend=:bottomleft)

plot!(mn2_sin[index], ribbon=sn2_sin[index], linewidth=4, label="sine")
plot!(mn2_MSE[index], ribbon=sn2_MSE[index], linewidth=4, label="MSE")
plot!(mn2_EFE[index], ribbon=sn2_EFE[index], linewidth=4, label="EFE")

savefig("./experiments/figures/l2dist-idARcoeffs.png")