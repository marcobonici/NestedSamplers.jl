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
    input_array = Array(chain.value)[:,1:end-1,1]
    len = length(input_array[:,1])
    logl_array = zeros(len)
    for i in 1:len
        logl_array[i] = logl(input_array[i,:])
    end
    return logl_array
end
