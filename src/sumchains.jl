function original_logw(logw_array, logZ)
    for i in 1:length(logw_array)
        logw_array[i] = log(logw_array[i]*exp(logZ))
    end
    return logw_array
end

function eval_logz(logw_array)
    myz = logw_array[1]
    for i in 2:length(logw_array)
        myz = logaddexp(myz, logw_array[i])
    end
    return myz
end

function eval_logz(chain::Chains, nactive)
    return eval_logz(eval_logw(original_logl(original_logw(Array(chain.value[:,end,1]),
    chain.logevidence), nactive), nactive))
end

function eval_logw(logl_array, nactive)
    logw_array = similar(logl_array)
    logvol = log1mexp(-1 / nactive)
    for i in 1:length(logw_array)
        logw_array[i] = logl_array[i]+ logvol
        logvol -= 1 / nactive
    end
    return logw_array
end

function original_logl(logw_array, nactive)
    logl_array = similar(logw_array)
    logvol = log1mexp(-1 / nactive)
    for i in 1:length(logw_array)
        logl_array[i] = logw_array[i]-logvol
        logvol -= 1 / nactive
    end

    return logl_array
end

function eval_logl(chain::Chains, logl)
    input_array = Array(chain.value)[:, 1:end-1, 1]
    len = length(input_array[:,1])
    logl_array = zeros(len)
    for i in 1:len
        logl_array[i] = logl(input_array[i,:])
    end
    return logl_array
end

function merge_chains(chains_array,
    sampler::Nested,
    check_wsum=true,
    kwargs...)
    num_chains = length(chains_array)
    cat_chains = Array(chains_array[1].value)
    for i in 2:num_chains
        cat_chains = cat(cat_chains, Array(chains_array[i].value), dims = 1)
    end
    cat_chains = cat_chains[sortperm(cat_chains[:, end]), :]
    logvol = log1mexp(-1 / (sampler.nactive*num_chains))
    logl = cat_chains[1,end,1]
    logwt = logl + logvol
    cat_chains[1,end-1,1] = logwt
    logz = logwt
    logvol -= 1 / (sampler.nactive*num_chains)
    for i in 2:length(cat_chains[:,1,1])
        logwt = logvol + cat_chains[i,end,1]
        cat_chains[i, end-1, 1] = logwt
        logz = logaddexp(logz, logwt)
        logvol -= 1 / (sampler.nactive*num_chains)
    end
    result = 0.
    n_live = sampler.nactive*num_chains
    for i in 1:length(cat_chains[:,1,1])
        h = exp(-(i-1)/n_live)-exp(-i/n_live)
        result += exp(cat_chains[i,end,1])*h
    end
    return log(result)

end
