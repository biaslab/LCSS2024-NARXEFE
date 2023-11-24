using Revise
using JLD
using LinearAlgebra
using Statistics
using Distributions
using LaTeXStrings
using Plots
default(label="", linewidth=5, margin=15Plots.pt)

includet("../location_scale_tdist.jl"); using .location_scale_tdist

num_exps = 100
N = 100
H = 2
Mu = 1
My = 1

ll_rnd = zeros(num_exps, N)
ll_QCR = zeros(num_exps, N)
ll_EFE = zeros(num_exps, N)
ll_sin = zeros(num_exps, N)

n2_rnd = zeros(num_exps, N)
n2_QCR = zeros(num_exps, N)
n2_EFE = zeros(num_exps, N)
n2_sin = zeros(num_exps, N)

for nn in 1:num_exps

    sys_settings = load("learningrates/data/system-pol$H-My$My-Mu$Mu-$nn.jld")
    sys_theta = sys_settings["sys_theta"]

    # DD = load("learningrates/results/learningrate-pol$H-My$My-Mu$Mu-rnd-$nn.jld")
    # ll_rnd[nn,:] = [location_scale_tdist.pdf(MvLocationScaleT(2DD["alpha"][k], DD["mu"][k], inv(DD["alpha"][k]/DD["beta"][k]*DD["Lambda"][k])), sys_theta) for k in 1:N]
    # n2_rnd[nn,:] = [norm(sys_theta - DD["mu"][k]) for k in 1:N]
    
    DD = load("learningrates/results/learningrate-pol$H-My$My-Mu$Mu-QCR-$nn.jld")
    # pθD = MvLocationScaleT(2DD["alpha"][k], DD["mu"][k], inv(DD["alpha"][k]/DD["beta"][k]*DD["Lambda"][k]))
    ll_QCR[nn,:] = [location_scale_tdist.pdf(MvLocationScaleT(2DD["alpha"][k], DD["mu"][k], inv(DD["alpha"][k]/DD["beta"][k]*DD["Lambda"][k])), sys_theta) for k in 1:N]
    n2_QCR[nn,:] = [norm(sys_theta - DD["mu"][k], 2) for k in 1:N]

    DD = load("learningrates/results/learningrate-pol$H-My$My-Mu$Mu-EFE-$nn.jld")
    # pθD = MvLocationScaleT(2DD["alpha"][k], DD["mu"][k], inv(DD["alpha"][k]/DD["beta"][k]*DD["Lambda"][k]))
    ll_EFE[nn,:] = [location_scale_tdist.pdf(MvLocationScaleT(2DD["alpha"][k], DD["mu"][k], inv(DD["alpha"][k]/DD["beta"][k]*DD["Lambda"][k])), sys_theta) for k in 1:N]
    n2_EFE[nn,:] = [norm(sys_theta - DD["mu"][k], 2) for k in 1:N]
    
    # DD = load("learningrates/results/learningrate-pol$H-My$My-Mu$Mu-sin-$nn.jld")
    # # pθD = MvLocationScaleT(2DD["alpha"][k], DD["mu"][k], inv(DD["alpha"][k]/DD["beta"][k]*DD["Lambda"][k]))
    # ll_sin[nn,:] = [location_scale_tdist.pdf(MvLocationScaleT(2DD["alpha"][k], DD["mu"][k], inv(DD["alpha"][k]/DD["beta"][k]*DD["Lambda"][k])), sys_theta) for k in 1:N]
    # n2_sin[nn,:] = [norm(sys_theta - DD["mu"][k], 2) for k in 1:N]
end

mll_rnd = mean(ll_rnd, dims=1)'
mll_sin = mean(ll_sin, dims=1)'
mll_QCR = mean(ll_QCR, dims=1)'
mll_EFE = mean(ll_EFE, dims=1)'
sll_rnd = std( ll_rnd, dims=1)' ./ num_exps
sll_sin = std( ll_sin, dims=1)' ./ num_exps
sll_QCR = std( ll_QCR, dims=1)' ./ num_exps
sll_EFE = std( ll_EFE, dims=1)' ./ num_exps

index = 1:100

plot(xlabel="time (#steps)", 
    grid=true, 
    guidefontsize=12, 
    tickfontsize=12, 
    legendfontsize=12, 
    ylabel=L"p(\theta^{*} | \mathcal{D}_k)", 
    size=(450,200), 
    xscale=:log10,
    yscale=:log10,
    legend=:bottomright)

# plot!(rate_rnd, linewidth=5, label="rnd")
# plot!(mll_sin[index], ribbon=sll_sin[index], linewidth=4, color="orange", label="sine")
plot!(mll_QCR[index], ribbon=sll_QCR[index], linewidth=4, color="blue", label="")
plot!(mll_EFE[index], ribbon=sll_EFE[index], linewidth=4, color="red", label="")

savefig("learningrates/figures/probpost-pol$H-My$My-Mu$Mu-idARcoeffs.png")
savefig("learningrates/figures/probpost-pol$H-My$My-Mu$Mu-idARcoeffs.pdf")

mn2_rnd = mean(n2_rnd, dims=1)'
mn2_sin = mean(n2_sin, dims=1)'
mn2_QCR = mean(n2_QCR, dims=1)'
mn2_EFE = mean(n2_EFE, dims=1)'
sn2_rnd = std( n2_rnd, dims=1)' ./ num_exps
sn2_sin = std( n2_sin, dims=1)' ./ num_exps
sn2_QCR = std( n2_QCR, dims=1)' ./ num_exps
sn2_EFE = std( n2_EFE, dims=1)' ./ num_exps

plot(xlabel="time (#steps)", 
    grid=true, 
    guidefontsize=12, 
    tickfontsize=12, 
    legendfontsize=12, 
    ylabel=L"|| \theta^{*} - \mu_k ||_2", 
    size=(600,200), 
    # xscale=:log10,
    legend=:bottomleft)

plot!(mn2_sin[index], ribbon=sn2_sin[index], linewidth=4, label="sine")
plot!(mn2_QCR[index], ribbon=sn2_QCR[index], linewidth=4, label="QCR")
plot!(mn2_EFE[index], ribbon=sn2_EFE[index], linewidth=4, label="EFE")

savefig("./experiments/figures/l2dist-idARcoeffs.png")