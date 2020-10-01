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

"""
Filter functions should be in the form (t,x) -> func(t,x)

The concatenated file, the sweeps are removed and replaced with traces. 
If the traces are different length, they are padded with 0's. 
The kwarg pad controls whether or not the padding is added to the beginning (:pre)
or the end (:post)

Prestim_time sets the amount of time (in seconds) before the END of the stimulus. This sets it so the effective time is always the prestim time

T_cutoff truncates the data to the time (in seconds)
"""
function concat(path_arr; t_cutoff = 3.5, t_eff = 0.5, filter_func = clean_data, sweep_avg = true, pad = :post)
    abfs = map(p -> extract_abf(p)[1:2], path_arr)
    n_traces = length(path_arr)
    
    dt = abfs[1][1][2] - abfs[1][1][1]
    t = collect(0.0:dt:(t_cutoff+t_eff))
    concatenated_trace = zeros(n_traces, length(t), 3)
    #Average multiple traces
    for (i, (t, raw_data)) in enumerate(abfs)
        print(i)
        if sweep_avg
            data = sum(raw_data, dims = 1)/size(raw_data,1)
        else
            data = raw_data
        end
        if filter_func == nothing
            x_ch1, x_ch2, x_stim = data[1,:,:] 
            x_stim = x_stim .> 0.2
        else
            x_ch1, x_ch2, x_stim = filter_func(t, data)
        end
                
        t_stim_end = findall(x -> x == true, x_stim)[end]
        t_start = t_stim_end - (t_eff/dt) |> Int64
        t_end = t_stim_end + (t_cutoff/dt) |> Int64
        
        concatenated_trace[i, :, 1] = x_ch1[t_start:t_end] 
        concatenated_trace[i, :, 2] = x_ch2[t_start:t_end] 
        concatenated_trace[i, :, 3] = x_stim[t_start:t_end] 
    end 
    t, concatenated_trace
end

#########################################Everything Below here is for Pepperburg analysis
"""
This function normalizes data and sets the minimum of the highest intensity nose component
"""
function norm_nose(data)
    n_swps, n_points, n_chs = size(data)
    norm_data = zeros(n_swps, n_points, n_chs)
    for i_ch in 1:n_chs
        norm_data[:,:,i_ch] = data[:,:,i_ch] ./ -minimum(data[:,:,i_ch])
    end
    return norm_data
end

"""
This function reorders data based on the minimum value in the sweeps
"""
function reorder_data(data; dim = 2, rev = false)
    min_sweeps = minimum(data, dims = 2)[:, 1, 1]
    if rev
        return data[sortperm(min_sweeps), :, :];
    else
        return data[sortperm(min_sweeps)|>reverse, :, :];
    end
end

"""
This function uses a histogram method to find the Rmax. 
"""
function peak_finder(X::AbstractArray; rmax_thresh = -0.50, precision = 500)
    if minimum(X) < rmax_thresh
        bins = LinRange(minimum(X[X .< rmax_thresh]), maximum(X[X .< rmax_thresh]), precision)
        h = fit(Histogram, X[X .< rmax_thresh], bins)
        return collect(h.edges...)[argmax(h.weights)]
    else
        #If none of the points are above the rmax_thresh, then the analysis fails and returns nothing
        return nothing
    end
end

"""
This function conducts a Pepperburg analysis on a single trace. 
"""
function pepperburg_analysis(X::AbstractArray; dt = 5.0e-5, rank = 6, graphically = false, peak_args...)
    rmax = peak_finder(X; peak_args...)
    if rmax != nothing
        #Now we need to find the values at 60% of the rmax found here (otherwise known as rank 6)
        rmax_idx = findall(x -> round(x, digits = 4) < round(rmax, digits = 4), X)[end]
        if length(rmax_idx) == 0
            #println("this is a fucked B-wave that hits the Rmax, but never returns")
            return nothing
        end
        rmax_idx = rmax_idx[end]
        rmax_rank = rmax*(rank/10)
        rridx = findall(x -> round(x, digits = 5) > round(rmax_rank, digits = 5), X)
        rmax_rank_idx = rridx[rridx .> rmax_idx][1]
        if graphically
            return rmax, rmax_idx, rmax_rank, rmax_rank_idx
        else
            rmax_time = rmax_idx * dt
            rmax_rank_time = rmax_rank_idx * dt
            return rmax_rank_time - rmax_time
        end
    else
        #return NaN
    end
end