using Distributions
using AbstractMCMC
using MCMCChains: Chains
using StatsFuns
using StatsBase

@testset "Bundles" begin
    logl(x::AbstractVector{T}) where T =  exp(-x[1]^2 / 2) / √(2π)
    priors = [Uniform(-1, 1)]
    model = NestedModel(logl, priors)
    spl = Nested(10)
    chain = sample(model, spl; dlogz = 0.2, param_names = ["x"], chain_type = Chains)
    samples = sample(model, spl; dlogz = 0.2, chain_type = Array)
end

# @testset "Flat" begin
#     logl(x::AbstractVector{T}) where T = zero(T)
#     priors = [Uniform(0, 1)]
#     model = NestedModel(logl, priors)

#     for method in [:single, :multi]
#         spl = Nested(4, method = method)
#         chain = sample(model, spl, dlogz = 0.2, chain_type = Array)

#         @test spl.logz ≈ 0 atol = 1e-9 # TODO
#         @test spl.h ≈ 0 atol = 1e-9 # TODO
#     end
# end

@testset "Gaussian" begin
    σ = 0.1
    μ1 = ones(2)
    μ2 = -ones(2)
    inv_σ = diagm(0 => fill(1 / σ^2, 2))

    function logl(x)
        dx1 = x .- μ1
        dx2 = x .- μ2
        f1 = -dx1' * (inv_σ * dx1) / 2
        f2 = -dx2' * (inv_σ * dx2) / 2
        return logaddexp(f1, f2)
    end

    priors = [Uniform(-5, 5), Uniform(-5, 5)]
    model = NestedModel(logl, priors)
    
    analytic_logz = log(2 * 2π * σ^2 / 100)

    for method in [:single, :multi]
        spl = Nested(100, method = method)
        chain = sample(model, spl, dlogz = 0.1, chain_type = Array)

        @test spl.logz ≈ analytic_logz atol = 2sqrt(spl.h / spl.nactive) # within 2sigma
        @test sort!(findpeaks(chain[:, 1, 1])[1:2]) ≈ [-1, 1] rtol = 3e-2
        @test sort!(findpeaks(chain[:, 2, 1])[1:2]) ≈ [-1, 1] rtol = 3e-2
    end
end