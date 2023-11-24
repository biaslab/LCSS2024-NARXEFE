using Revise
using FFTW
using DSP
using SpecialFunctions
using LinearAlgebra
using ProgressMeter
using Distributions
using JLD

includet("../NARXAgents.jl"); using .NARXAgents
includet("../NARXsystem.jl"); using .NARXsystem
includet("../mv_normal_gamma.jl"); using .mv_normal_gamma 

# Experimental parameters
num_exps = 100
N = 100
fs = 20 # Hertz
Δt = 1/fs
tsteps = collect(range(0.0, step=Δt, length=N))
T = 10

# Basis
H = 2
sys_basis(x) = cat([1.0; [x.^d for d in 1:H]]...,dims=1)
Mu = 2
My = 2
M = size(sys_basis(zeros(My + 1 + Mu)),1)

# Control parameters
input_lims = (-1.,1.)

# Specify prior distributions
α0 = 10.0
β0 = 1e-1
μ0 = 1e-8*ones(M)
Λ0 = 1e-1diagm(ones(M))
goals = [Normal(0.5, 1.0) for t in 1:T]

# ppp = Progress(num_exps)
@showprogress for nn in 1:num_exps

    # Load system parameters
    sys_settings = load("learningrates/data/system-pol$H-My$My-Mu$Mu-$nn.jld")
    sys_mnoise_sd = sys_settings["sys_sd"]
    sys_coefficients = sys_settings["sys_theta"]

    # Sinusoidal controls
    Ω  = rand(5)*10
    controls = mean([sin.(ω.*tsteps) for ω in Ω]);

    # Start system
    system = NARXsys(sys_coefficients, 
                     sys_basis, 
                     sys_mnoise_sd, 
                     order_outputs=My, 
                     order_inputs=Mu, 
                     input_lims=input_lims)

    py = []
    μ = [μ0]
    Λ = [Λ0]
    α = [α0]
    β = [β0]
    FE = zeros(N)

    agent = NARXAgent(μ0, Λ0, α0, β0,
                      goal_prior=goals, 
                      delay_inp=Mu, 
                      delay_out=My, 
                      pol_degree=H)

    outputs = zeros(N)
    inputs = zeros(N)
    inputs_ = [inputs; zeros(T)]

    for k in 1:N

        # Evolve system
        NARXsystem.update!(system, controls[k])
        outputs[k] = system.observation
        inputs[k] = system.input_buffer[1]
        
        # Make predictions
        push!(py, predictions(agent, inputs_[k:k+T], time_horizon=T))
        
        # Update beliefs
        NARXAgents.update!(agent, outputs[k], inputs[k])
        
        push!( μ, agent.μ )
        push!( Λ, agent.Λ )
        push!( α, agent.α )
        push!( β, agent.β )

        FE[k] = agent.free_energy
    end

    save("learningrates/results/learningrate-pol$H-My$My-Mu$Mu-sin-$nn.jld", "py", py, "y_sin", outputs, "u_sin", inputs, "FE", FE,
        "mu", μ, "Lambda", Λ, "alpha", α, "beta", β, 
        "mu0", μ0, "Lambda0", Λ0, "alpha0", α0, "beta0", β0, 
        "goal", goals, "u_lims", input_lims, "thorizon", T)

    # next!(ppp)
end
# finish!(ppp)
