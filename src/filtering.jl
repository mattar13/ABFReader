####################These functions are for filtering and adjusting the traces################
"""
This function adjusts the baseline, similar to how it is done in clampfit. 
    it can cancel baseline based on: 
    - The mean of a region
    - The mean of the whole trace
    - The slope of a region
    - The slope of the whole trace
"""
function baseline_cancel(trace::NeuroTrace; mode::Symbol = :mean, region = :prestim, cust_rng = (1,10))
    if isa(region, Tuple{Float64, Float64})
        rng_begin = round(Int, region[1]/dt)+1
        if region[2] > trace.t[end]
            rng_end = length(trace.t)
        else
            rng_end = round(Int, region[2]/dt)+1
        end
        rng_end_arr = nothing
    elseif isa(region, Tuple{Int64, Int64})
        rng_begin, rng_end = region
        rng_end_arr = nothing
    elseif region == :whole
        rng_begin = 1
        rng_end = length(trace)
        rng_end_arr = nothing
    elseif region == :prestim
        rng_begin = 1
        rng_end = nothing
        rng_end_arr = findstimRng(trace)[:, 1]
    end

    data = copy(trace)
    if mode == :mean
        for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                #never adjust the stim
                if ch != trace.stim_ch
                    if rng_end_arr != nothing
                        rng_end = rng_end_arr[swp]
                    end
                    baseline_adjust = sum(trace[swp, rng_begin:rng_end, ch])/(rng_end-rng_begin)
                    data[swp,:, ch] .= trace.data_array[swp,:,ch] .- baseline_adjust
                else
                    data[swp,:,ch] .= trace[swp,:,ch]
                end
            end
        end
    elseif mode == :slope
		for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                #never adjust the stim
                if ch != trace.stim_ch
                    if rng_end_arr != nothing
                        rng_end = rng_end_arr[swp]
                    end
                    pfit = Polynomials.fit(trace.t[rng_begin:rng_end], trace[swp, rng_begin:rng_end , ch], 1)
                    data[swp, :, ch] .= trace[swp, :, ch] - pfit.(trace.t)
                else
                    data[swp,:,ch] .= trace[swp,:,ch]
                end
            end
        end
	else
		data = trace.data_array
    end
    return data
end

function baseline_cancel!(trace::NeuroTrace; mode::Symbol = :mean, region = :prestim)
    if isa(region, Tuple{Float64, Float64})
        rng_begin = round(Int, region[1]/trace.dt)+1
        if region[2] > trace.t[end]
            rng_end = length(trace.t)
        else
            rng_end = round(Int, region[2]/trace.dt)+1
        end
        rng_end_arr = nothing
    elseif isa(region, Tuple{Int64, Int64})
        rng_begin, rng_end = region
        rng_end_arr = nothing
    elseif region == :whole
        rng_begin = 1
        rng_end = length(trace)
        rng_end_arr = nothing
    elseif region == :prestim
        rng_begin = 1
        rng_end = nothing
        rng_end_arr = findstimRng(trace)[:, 1]
    end
    
    if mode == :mean
        for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                #never adjust the stim
                if ch != trace.stim_ch
                    if rng_end_arr != nothing
                        rng_end = rng_end_arr[swp]
                    end
                    baseline_adjust = sum(trace[swp, rng_begin:rng_end, ch])/(rng_end-rng_begin)
                    trace.data_array[swp,:, ch] .= trace.data_array[swp,:,ch] .- baseline_adjust
                end
            end
        end
    elseif mode == :slope
		 for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                #never adjust the stim
                if ch != trace.stim_ch
                    if rng_end_arr != nothing
                        rng_end = rng_end_arr[swp]
                    end
                    pfit = Polynomials.fit(trace.t[rng_begin:rng_end], trace[swp, rng_begin:rng_end , ch], 1)
                    #println(ch + pfit.(trace.t) |> size)
                    trace.data_array[swp, :, ch] .= trace[swp, :, ch] - pfit.(trace.t)
                end
            end
        end
    end
end

"""
This function applies a n-pole lowpass filter
"""
function lowpass_filter(trace::NeuroTrace; freq = 40.0, pole = 8)
    
    responsetype = Lowpass(freq; fs =  1/trace.dt)
    designmethod = Butterworth(8)
    digital_filter = digitalfilter(responsetype, designmethod)
    data = copy(trace)
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                data[swp,:,ch] .= filt(digital_filter, trace[swp, :, ch])
            else
                data[swp,:,ch] .= trace[swp,:,ch]
            end
        end
    end
    return data
end

function lowpass_filter!(trace::NeuroTrace; freq = 40.0, pole = 8)
    
    responsetype = Lowpass(freq; fs =  1/trace.dt)
    designmethod = Butterworth(8)
    digital_filter = digitalfilter(responsetype, designmethod)
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                trace.data_array[swp,:,ch] .= filt(digital_filter, trace[swp, :, ch])
            end
        end
    end
end

lowpass_filter(trace::NeuroTrace, freq; pole = 8) = lowpass_filter(trace; freq = freq, pole = pole)

function notch_filter(trace::NeuroTrace; pole = 8, center = 60.0, std = 0.1)
    
    responsetype = Bandstop(center-std, center+std; fs = 1/trace.dt)
	designmethod = Butterworth(8)
	digital_filter = digitalfilter(responsetype, designmethod)
    data = copy(trace)
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                data[swp,:,ch] .= filt(digital_filter, trace[swp, :, ch])
            else
                data[swp,:,ch] .= trace[swp,:,ch]
            end
        end
    end
    return data
end

function notch_filter!(trace::NeuroTrace; pole = 8, center = 60.0, std = 0.1)
    
    responsetype = Bandstop(center-std, center+std; fs = 1/trace.dt)
	designmethod = Butterworth(8)
	digital_filter = digitalfilter(responsetype, designmethod)
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                trace.data_array[swp,:,ch] .= filt(digital_filter, trace[swp, :, ch])
            end
        end
    end
end

function cwt_filter(trace::NeuroTrace; wave = WT.dog2, periods = 1:9, return_cwt = true)
    
    data = similar(trace.data_array)
    stim_begin, stim_end = findstimRng(trace)
    data = copy(trace)
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                y = cwt(trace[swp, :, ch], wavelet(wave))
                data[swp,:,ch] .= sum(real.(y[:,periods]))/length(y);
            else
                data[swp,:,ch] .= trace[swp,:,ch]
            end
        end
    end
    data
end

function cwt_filter!(trace::NeuroTrace; wave = WT.dog2, periods = 1:9)
    
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                y = cwt(trace[swp, :, ch], wavelet(wave))
                trace.data_array[swp,:,ch] .= sum(real.(y[:,periods]))/length(y);
            end
        end
    end
end

"""
If the traces contain multiple runs, then this file averages the data
"""
function average_sweeps(trace::NeuroTrace)
    
    data = copy(trace)
    for ch in 1:size(trace,3)
        data[:,:,ch] .= sum(trace[:,:,ch], dims = 1)/size(trace,1)
    end
    return data
end

average_sweeps!(trace::NeuroTrace) = trace.data_array = sum(trace, dims = 1)/size(trace,1) 

function normalize(trace::NeuroTrace; rng = (-1,0))
    data = copy(trace)
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                data[swp,:,ch] .= (trace[swp,:,ch] ./ minimum(trace[swp,:,ch], dims = 2))
            end
        end
    end
    return data
end

function normalize!(trace::NeuroTrace; rng = (-1,0))
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #never adjust the stim
            if ch != trace.stim_ch
                trace[swp,:,ch] .= (trace[swp,:,ch] ./ minimum(trace[swp,:,ch], dims = 2))
            end
        end
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
This function removes the stimulus artifact. 
The estimated duration of the stimulus artifact is 0.0025s or 25us
"""
function remove_artifact(t, data; est_duration = 0.0025, artifact_val = :mean)
	dt = t[2]-t[1]
	x_ch1  = data[1,:,1] 
	x_ch2  = data[1,:,2] 
	x_stim = data[1,:,3] .> 0.2
	offset = round(Int,est_duration/dt)
	t_stim_start = findall(x -> x == true, x_stim)[1]
	t_stim_end = findall(x -> x == true, x_stim)[end]
	stim_snip_ch1 = x_ch1[1:t_stim_start] 
	stim_snip_ch2 = x_ch2[1:t_stim_start]
	
	artifact_thresh_ch1 = (sum(stim_snip_ch1)/length(stim_snip_ch1))
	artifact_thresh_ch2 = (sum(stim_snip_ch2)/length(stim_snip_ch2))
	if artifact_val == :mean
        data[1,t_stim_start:(t_stim_start+offset),1] .= artifact_thresh_ch1
        data[1,t_stim_start:(t_stim_start+offset),2] .= artifact_thresh_ch2
        data[1,t_stim_end:(t_stim_end+offset),1] .= artifact_thresh_ch1
        data[1,t_stim_end:(t_stim_end+offset),2] .= artifact_thresh_ch2
        return t, data
    else
        data[1,t_stim_start:(t_stim_start+offset),1] .= artifact_val
        data[1,t_stim_start:(t_stim_start+offset),2] .= artifact_val
        data[1,t_stim_end:(t_stim_end+offset),1] .= artifact_val
        data[1,t_stim_end:(t_stim_end+offset),2] .= artifact_val
        return t, data
    end
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