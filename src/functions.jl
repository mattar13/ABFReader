"""
This function is for computing the R-squared of a polynomial
"""
function RSQ(poly::Polynomial, x, y)
	ŷ = poly.(x)
	ȳ = sum(ŷ)/length(ŷ)
	SSE = sum((y-ŷ).^2)
	SST = sum((y.-ȳ).^2)
	1-SSE/SST
end

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
    if region == :whole
        rng_begin = 1
        rng_end = length(trace)
    elseif region == :prestim
        rng_begin = 1
        rng_end = findstimRng(trace)[1]
    elseif region == :custom
        rng_begin, rng_end = cust_rng
    end

    if mode == :mean
        baseline_adjust = sum(trace[:,rng_begin:rng_end,:]; dims = 2)/(rng_end-rng_begin)
        data = trace.data_array .- baseline_adjust
    elseif mode == :slope
		data = similar(trace.data_array)
		for (i,ch) in enumerate(eachchannel(trace))
        	pfit = Polynomials.fit(trace.t, ch, 1)
			println(ch + pfit.(trace.t) |> size)
        	data[:, :, i] = ch - pfit.(trace.t)
		end
	else
		data = trace.data_array
    end
    return NeuroTrace(
		trace.t, 
		data, #Add the data here 
		trace.date_collected,
		trace.tUnits,
		trace.dt,
		trace.chNames,
		trace.chUnits,
		trace.labels,
		trace.stim_ch,
		)
end

function baseline_cancel!(trace::NeuroTrace, t, x_data::AbstractArray; return_fit = false)
   if region == :whole
        rng_begin = 1
        rng_end = length(trace)
    elseif region == :prestim
        rng_begin = 1
        rng_end = findstimRng(trace)[1]
    elseif region == :custom
        rng_begin, rng_end = cust_rng
    end

     if mode == :mean
        baseline_adjust = sum(trace[:,rng_begin:rng_end,:]; dims = 2)/(rng_end-rng_begin)
        trace.data_array = trace.data_array .- baseline_adjust
    elseif mode == :slope
		for (i,ch) in enumerate(eachchannel(trace))
        	pfit = Polynomials.fit(trace.t, ch, 1)
			println(ch + pfit.(trace.t) |> size)
        	trace[:,rng_begin:rng_end,i] .= ch - pfit.(trace.t)
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
    data = similar(trace.data_array)
    for (i,ch) in enumerate(eachchannel(trace))
		if i != trace.stim_ch
        	data[:,:,i] = filt(digital_filter, getchannel(trace,i))
		else
			data[:,:,i] = trace[:,:,i]
		end
    end
    return NeuroTrace(
        trace.t, 
        data, #Add the data here 
        trace.date_collected,
        trace.tUnits,
        trace.dt,
        trace.chNames,
        trace.chUnits,
        trace.labels,
        trace.stim_ch,
        )	
end

function lowpass_filter!(trace::NeuroTrace; freq = 40.0, pole = 8)
    responsetype = Lowpass(freq; fs =  1/trace.dt)
    designmethod = Butterworth(8)
    digital_filter = digitalfilter(responsetype, designmethod)
    for (i,ch) in enumerate(eachchannel(trace))
        if i != trace.stim_ch
        	trace[:,:,i] = filt(digital_filter, getchannel(trace,i))
		else
			trace[:,:,i] = trace[:,:,i]
		end
    end
end

lowpass_filter(trace::NeuroTrace, freq; pole = 8) = lowpass_filter(trace; freq = freq, pole = pole)

function notch_filter(trace::NeuroTrace; pole = 8, center = 60.0, std = 0.1)
	responsetype = Bandstop(center-std, center+std; fs = 1/trace.dt)
	designmethod = Butterworth(8)
	digital_filter = digitalfilter(responsetype, designmethod)
    for (i,ch) in enumerate(eachchannel(trace))
		if i != trace.stim_ch
        	data[:,:,i] = filt(digital_filter, getchannel(trace,i))
		else
			data[:,:,i] = trace[:,:,i]
		end
    end
    return NeuroTrace(
        trace.t, 
        data, #Add the data here 
        trace.date_collected,
        trace.tUnits,
        trace.dt,
        trace.chNames,
        trace.chUnits,
        trace.labels,
        trace.stim_ch,
        )		
end

function notch_filter!(trace::NeuroTrace; pole = 8, center = 60.0, std = 0.1)
    responsetype = Bandstop(center-std, center+std; fs = 1/trace.dt)
	designmethod = Butterworth(8)
	digital_filter = digitalfilter(responsetype, designmethod)
    for (i,ch) in enumerate(eachchannel(trace))
        if i != trace.stim_ch
        	trace[:,:,i] = filt(digital_filter, getchannel(trace,i))
		else
			trace[:,:,i] = trace[:,:,i]
		end
    end
end

function cwt_filter(trace::NeuroTrace; wave = WT.dog2, periods = 1:9, return_cwt = true)
    data = similar(trace.data_array)
    for (i,ch) in enumerate(eachchannel(trace))
        y = cwt(ch, wavelet(wave))
        data[:,:,i] = sum(real.(y[:,periods]), dims = 2)/size(y, 2);
    end
    return NeuroTrace(
        trace.t, 
        data, #Add the data here 
        trace.date_collected,
        trace.tUnits,
        trace.dt,
        trace.chNames,
        trace.chUnits,
        trace.labels,
        trace.stim_ch,
        )
end

function cwt_filter!(trace::NeuroTrace; wave = WT.dog2, periods = 1:9)
    data = similar(trace.data_array)
	for (i,ch) in enumerate(eachchannel(trace))
        y = cwt(ch, wavelet(wave))
        trace[:,:,i] = sum(real.(y[:,periods]), dims = 2)/size(y, 2);
	end
end

"""
If the traces contain multiple runs, then this file averages the data
"""
function average_sweeps(nt::NeuroTrace)
    data = sum(nt, dims = 1)/size(nt,1)
    new_obj = copy(nt)
    new_obj.data_array = data
    return new_obj
end

average_sweeps!(nt::NeuroTrace) = nt.data_array = sum(nt, dims = 1)/size(nt,1)

function normalize(trace::NeuroTrace; rng = (-1,0))
    x = trace.data_array
    z = (x .- minimum(x))/(maximum(x) - minimum(x))*(rng[2] - rng[1]) .+ rng[1]
    new_obj = copy(nt)
    new_obj.data_array = z
    return new_obj
end

function normalize!(trace::NeuroTrace; rng = (-1,0))
    x = trace.data_array
    trace.data_array = (x .- minimum(x))/(maximum(x) - minimum(x))*(rng[2] - rng[1]) .+ rng[1]
end

function fft_spectrum(t, data::Array{T, 1}) where T <: Real
    #FFTW filtering
    x_fft = fft(data) |> fftshift
    dt = t[2] - t[1]
    freqs = FFTW.fftfreq(length(t), 1.0/dt) |> fftshift
    over_0 = findall(freqs .> 0);
    return freqs[over_0], x_fft[over_0] 
end

#A single filtering pipeline
function clean_data_cwt(t, data; negative = true, wave = WT.dog2, cutoff_octave = 9)
    x_ch1 = data[1,:,1]; x_ch2 = data[1,:,2]; x_stim = data[1,:,3] .> 0.2;
    #Cancelling drift
    x_lin1 = drift_cancel(t, x_ch1);
    x_lin2 = drift_cancel(t, x_ch2);
    #Baseline subtraction
    stim_idxs = findall(x -> x == true, x_stim) #Stimulus is same for both channels
    x_adj1 = subtract_baseline(x_lin1, (1, stim_idxs[1]));
    x_adj2 = subtract_baseline(x_lin2, (1, stim_idxs[1]));
    #Normalization
    x_norm1, norm_factor1 = normalize(x_adj1);
    x_norm2, norm_factor2 = normalize(x_adj2);
    #CWT filtering (Probably not ready for CWT filtering )
    x_cwt1, cwt1_raster = cwt_filter(x_norm1, periods = 1:cutoff_octave);
    x_cwt2, cwt2_raster = cwt_filter(x_norm2, periods = 1:cutoff_octave);
    return x_cwt1.*norm_factor1, x_cwt2.*norm_factor2, x_stim
    #return x_cwt1, x_cwt2, x_stim
end

#A single filtering pipeline
function clean_data(t, data; negative = true, wave = WT.dog2, cutoff_octave = 9)
    x_ch1 = data[1,:,1]; x_ch2 = data[1,:,2]; x_stim = data[1,:,3] .> 0.2;
    #Cancelling drift
    x_lin1 = drift_cancel(t, x_ch1);
    x_lin2 = drift_cancel(t, x_ch2);
    #Baseline subtraction
    stim_idxs = findall(x -> x == true, x_stim) #Stimulus is same for both channels
    x_adj1 = subtract_baseline(x_lin1, (1, stim_idxs[1]));
    x_adj2 = subtract_baseline(x_lin2, (1, stim_idxs[1]));
    #Normalization should be done to all data points after concatenation
    return x_adj1, x_adj2, x_stim
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
    if rmax !== nothing
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