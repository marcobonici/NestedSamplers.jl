# Interface Implementations for running full sampling loops
using Printf

StatsBase.sample(rng::AbstractRNG, model, sampler::Nested; kwargs...) =
    mcmcsample(rng, model, sampler; kwargs...)

StatsBase.sample(model, sampler::Nested; kwargs...) =
    StatsBase.sample(GLOBAL_RNG, model, sampler; kwargs...)

function mcmcsample(
    rng,
    model,
    sampler::Nested;
    progress=true,
    progressname="Nested Sampling",
    discard_initial=0,
    thinning=1,
    chain_type=Any,
    kwargs...
)
    # Obtain the initial sample and state.
    sample, state = step(rng, model, sampler; kwargs...)

    # Discard initial samples.
    for _ in 2:discard_initial
        # Obtain the next sample and state.
        sample, state = step(rng, model, sampler, state; kwargs...)
    end

    # Save the sample.
    samples = AbstractMCMC.samples(sample, model, sampler; kwargs...)
    samples = save!!(samples, sample, 1, model, sampler; kwargs...)

    # Step through the sampler until stopping.
    i = 2
        
    @ifwithprogress progress name=progressname begin
        while !nested_isdone(rng, model, sampler, samples, state, i; progress=progress, kwargs...)
            # Discard thinned samples.
            for _ in 1:(thinning - 1)
                # Obtain the next sample and state.
                _, state = step(rng, model, sampler, state; kwargs...)
            end
            # Obtain the next sample and state.
            sample, state = step(rng, model, sampler, state; kwargs...)

            # Save the sample.
            samples = save!!(samples, sample, i, model, sampler; kwargs...)
            # Increment iteration counter.
            i += 1
        end
        logz_remain = maximum(state.logl) + state.logvol
        delta_logz = logaddexp(state.logz, logz_remain) - state.logz
    end

    # Wrap the samples up.
    return bundle_samples(samples, model, sampler, state, chain_type; kwargs...), state
end

function nested_isdone(rng, model, smapler, samples, state, i; progress=true, maxiter=Inf, maxcall=Inf, dlogz=0.5, maxlogl=Inf, kwargs...)
    # 1) iterations exceeds maxiter
    done_sampling = state.it > maxiter
    # 2) number of loglike calls has been exceeded
    done_sampling |= state.ncall > maxcall
    # 3) remaining fractional log-evidence below threshold
    logz_remain = maximum(state.logl) + state.logvol
    delta_logz = logaddexp(state.logz, logz_remain) - state.logz
    done_sampling |= delta_logz < dlogz
    # 4) last dead point loglikelihood exceeds threshold
    done_sampling |= state.logl_dead > maxlogl
    # 5) number of effective samples
    # TODO

    if progress
        str = @sprintf "\33[2K\rsampling: iter=%d\tncall=%d\tΔlogz=%.2g\tlogl=%.2g\tlogz=%.2g" i state.ncall delta_logz state.logl_dead state.logz
        print(str)
    end

    return done_sampling
end
