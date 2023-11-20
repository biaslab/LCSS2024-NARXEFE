module location_scale_tdist

using LinearAlgebra
using Distributions
using SpecialFunctions

export LocationScaleT, MvLocationScaleT, pdf, logpdf, params, dims

mutable struct LocationScaleT <: ContinuousUnivariateDistribution
 
    ν ::Real
    μ ::Real
    σ ::Real

    function LocationScaleT(ν::Float64, μ::Float64, σ::Float64)
        
        if ν <= 0.0; error("Degrees of freedom parameter must be positive."); end
        if σ <= 0.0; error("Scale parameter must be positive."); end

        return new(ν, μ, σ)
    end
end

mutable struct MvLocationScaleT <: ContinuousMultivariateDistribution
 
    ν ::Real
    μ ::Vector
    Σ ::Matrix

    function MvLocationScaleT(ν::Float64, μ::Vector{Float64}, Σ::Matrix{Float64})
        
        if ν <= 0.0; error("Degrees of freedom parameter must be positive."); end
        if length(μ) !== size(Σ,1); error("Dimensionalities of mean and covariance matrix don't match."); end

        return new(ν, μ, Σ)
    end
end

params(p::LocationScaleT) = (p.ν, p.μ, p.σ)
params(p::MvLocationScaleT) = (p.ν, p.μ, p.Σ)

dims(p::MvLocationScaleT) = length(p.μ)

function pdf(p::LocationScaleT, x)
    ν, μ, σ = params(p)
    return gamma( (ν+1)/2 ) / ( gamma(ν/2) *sqrt(π*ν)*σ ) * ( 1 + (x-μ)^2/(ν*σ^2) )^( -(ν+1)/2 )
end

function pdf(p::MvLocationScaleT, x)
    d = dims(p)
    ν, μ, Σ = params(p)
    return sqrt(1/( (ν*π)^d*det(Σ) )) * gamma((ν+d)/2)/gamma(ν/2) * (1 + 1/ν*(x-μ)'*inv(Σ)*(x-μ))^(-(ν+d)/2)
end

function logpdf(p::LocationScaleT, x)
    ν, μ, σ = params(p)
    return loggamma( (ν+1)/2 ) - loggamma(ν/2) - 1/2*log(πν) - log(σ) + ( -(ν+1)/2 )*log( 1 + (x-μ)^2/(ν*σ^2) )
end

function logpdf(p::MvLocationScaleT, x)
    d = dims(p)
    ν, μ, Σ = params(p)
    return -d/2*log(ν*π) - 1/2*logdet(Σ) +loggamma((ν+d)/2) -loggamma(ν/2) -(ν+d)/2*log(1 + 1/ν*(x-μ)'*inv(Σ)*(x-μ))
end

end