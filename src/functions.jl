#This calculates the photons density
photons(E::Float64; λ::Int64 = 525) = ( E / (6.626e-34 * 3e8/(λ*10^-9)))*10e-8

####################These functions are for filtering and adjusting the traces################

function drift_cancel(t, x_data::AbstractArray; return_fit = false)
    pfit = Polynomials.fit(t, x_data, 1)
    x_lin = x_data - pfit.(t)
    if return_fit
        return x_lin, pfit
    else
        return x_lin
    end
end

function subtract_baseline(x_data::AbstractArray, subtract_range::Tuple{Int64, Int64}; return_val = false)
    rng_begin, rng_end = subtract_range
    baseline_adjust = sum(x_data[rng_begin:rng_end])/length(x_data[rng_begin:rng_end])
    x_adj = x_data .- baseline_adjust
    if return_val
        return x_adj, baseline_adjust
    else
        return x_adj
    end
end

function normalize(x_data; negative = true, return_val = true)
    if negative 
        norm_factor = minimum(x_data)
    else
        norm_factor = maximum(x_data)
    end
    x_norm = x_data / norm_factor
    if return_val
        return x_norm, norm_factor
    else
        return x_norm
    end
end

function cwt_filter(x_data; wave = WT.dog2, periods = 1:9, return_cwt = true)
    y = cwt(x_data, wavelet(wave))
    x_cwt = sum(real.(y[:,periods]), dims = 2);
    x_cwt ./= minimum(x_cwt);
    if return_cwt
        return vec(x_cwt), y
    else
        return vec(x_cwt)
    end
end

function fft_spectrum(t, data::Array{T, 1}) where T <: Real
    #FFTW filtering
    x_fft = fft(data) |> fftshift
    dt = t[2] - t[1]
    freqs = FFTW.fftfreq(length(t), 1.0/dt) |> fftshift
    over_0 = findall(freqs .> 0);
    return freqs[over_0], x_fft[over_0] 
end
