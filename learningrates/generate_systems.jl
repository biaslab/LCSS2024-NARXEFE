using JLD
using DSP
using Distributions
using ProgressMeter
using LinearAlgebra

# Experimental parameters
num_exps = 100
N = 100
fs = 20 # Hertz
Δt = 1/fs

# Basis
H = 2
sys_basis(x) = cat([1.0; [x.^d for d in 1:H]]...,dims=1)
Mu = 1
My = 1
M = size(sys_basis(zeros(My + 1 + Mu)),1)

# Randomized variables
pσ = Gamma(3.0, 1/300.)
pf = Uniform(0.5, 2.0)

ppp = Progress(num_exps)
for nn in 1:num_exps

    # Define system parameters
    sys_mnoise_sd = rand(pσ)
    sys_lowpass_f = rand(pf)
    df = digitalfilter(Lowpass(sys_lowpass_f; fs=fs), Butterworth(maximum([Mu, My])))
    sys_coefficients = [0.0; sys_basis([coefb(df)[2:My+1]; coefa(df)[1:Mu+1]])[2:end]]

    save("learningrates/data/system-pol$H-My$My-Mu$Mu-$nn.jld", "sys_sd", sys_mnoise_sd, "sys_theta", sys_coefficients)

    next!(ppp)
end
finish!(ppp)