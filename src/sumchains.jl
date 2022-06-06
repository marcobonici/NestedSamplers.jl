function merge_chains(chains_array,
    sampler::Nested,
    check_wsum=true;
    param_names=missing,
    kwargs...)
    nactive = sampler.nactive
    #merge all chains together
    num_chains = length(chains_array)
    cat_chains = Array(chains_array[1].value)
    for i in 2:num_chains
        cat_chains = cat(cat_chains, Array(chains_array[i].value), dims = 1)
    end
    #sort chains according to logl
    cat_chains = cat_chains[sortperm(cat_chains[:, end]), :]
    #evaluate logwt and logz
    logvol = log1mexp(-1 / (nactive*num_chains))
    logl = cat_chains[1,end,1]
    logwt = logl + logvol
    cat_chains[1,end-1,1] = logwt
    logz = logwt
    logvol -= 1 / (nactive*num_chains)
    for i in 2:length(cat_chains[:,1,1])
        logwt = logvol + cat_chains[i,end,1]
        cat_chains[i, end-1, 1] = logwt
        logz = logaddexp(logz, logwt)
        logvol -= 1 / (nactive*num_chains)
    end

    #normalize logwt, such that their sum is 1.0
    for i in 1:length(cat_chains[:,1,1])
        logwt = exp(cat_chains[i,end-1,1] - logz)
        cat_chains[i,end-1,1] = logwt
    end

    if check_wsum
        wsum = sum(cat_chains[:, end-1, 1])
        err = 1e-4 #add something similar to what is done in the standard approach
        isapprox(wsum, 1, atol=err) || @warn "Weights sum to $wsum instead of 1; possible bug"
    end

    names = deepcopy(param_names) #this to avoid problems when doing more runs in parallel
    # Parameter names
    if names === missing
        names = ["Parameter $i" for i in 1:length(cat_chains[1, :]) - 2]
    end
    push!(names, "weights", "logl")

    return Chains(cat_chains, names, Dict(:internals => ["weights", "logl"]), evidence=logz)

end
