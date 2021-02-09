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

function RSQ(ŷ, y)
	ȳ = sum(ŷ)/length(ŷ)
	SSE = sum((y-ŷ).^2)
	SST = sum((y.-ȳ).^2)
	1-SSE/SST
end


"""
This function calculates the min, max, mean, and std of each trace
"""
function calculate_basic_stats(data::Experiment)
    mins = minimum(data.data_array, dims = 2)[:,1,:]
    maxes = maximum(data.data_array, dims = 2)[:,1,:]
    means = zeros(size(data,1), size(data,3))
    stds = zeros(size(data,1), size(data,3))
    for swp in 1:size(data,1), ch in 1:size(data,3)
        stim_begin, stim_end = data.stim_protocol[swp].index_range
        pre_stim = data[:, 1:stim_begin, :]
        post_stim = data[:, stim_begin:size(data,2), :]
        means[swp, ch] = sum(pre_stim[swp, :, ch])/size(pre_stim,2)
        stds[swp, ch] = std(pre_stim[swp, :, ch])
    end
    return mins, maxes, means, stds
end

rolling_mean(arr::AbstractArray; radius = 5) = [sum(arr[i:i+radius])/radius for i = 1:length(arr)-radius]

"""
This function uses a histogram method to find the saturation point. 
    - In ERG traces, a short nose component is usually present in saturated values
    - Does this same function work for the Rmax of nonsaturated responses?
    - Setting the saturated threshold to infinity will completely disregard the histogram method
"""
function saturated_response(trace::Experiment{T}; saturated_thresh = :determine, polarity::Int64 = -1, precision = 500, z = 1.3, kwargs...) where T
    if isa(saturated_thresh, Symbol)
        #Figure out if the saturated threshold needs to be determined
        saturated_thresh = size(trace,1)/precision/2
    end
    #Make an empty array for recording the rmaxes
    rmaxs = T[]
    for ch in 1:size(trace,3)
        data = Float64[]
        for swp in 1:size(trace,1)
            if isempty(trace.stim_protocol)
                #This catch is here for if no stim protocol has been set
                #println("No stimulus protocol exists")
                stim_begin = 1
            else
                stim_begin = trace.stim_protocol[swp].index_range[1] #We don't want to pull values from before the stim
            end
            push!(data,  trace[:, stim_begin:size(trace,2), ch]...)
        end
        #We are going to concatenate all sweeps together into one histogram
        mean = sum(data)/length(data)
        deviation = z*std(data)
        #Here we cutoff all points after the sweep returns to the mean
        if polarity < 0
            idxs = findall(data .< (mean - deviation))
            if isempty(idxs)
                #This is a weird catch, but no points fall under the mean. 
                push!(rmaxs, minimum(data))
                continue
            end
            data = data[idxs]
            #For negative components
            bins = LinRange(minimum(data), min(0.0, mean-deviation), precision)
        elseif polarity > 0
            idxs = findlast(data .> (mean + deviation))
            if isempty(idxs)
                #This is a weird catch, but no points fall under the mean. 
                push!(rmaxsm, minimum(data))
                continue
            end
            data = data[idxs]
            #For positive components
            bins = LinRange(max(0.0, mean+deviation), maximum(data),  precision)
        else
            throw(error("Polarity incorrect"))
        end
        h = Distributions.fit(Histogram, data, bins; )
        edges = collect(h.edges...)[2:end]
        weights = h.weights./length(data)
        
        #println(maximum(weights))
        #println(saturated_thresh)
        #return edges, weights

        if maximum(weights) > saturated_thresh
            push!(rmaxs, edges[argmax(weights)])
        else
            push!(rmaxs, minimum(data))
        end
    end
    rmaxs |> vec
end

"""
This function only works on concatenated files with more than one trace
    Rmax argument should have the same number of sweeps and channels as the 
    In the rdim calculation, it is better to adjust the higher percent
    Example: no traces in 20-30% range, try 20-40%
"""
function dim_response(trace::Experiment{T}, rmaxes::Array{T, 1}; return_idx = true, polarity::Int64 = -1, rmax_lin = [0.10, 0.40]) where T <: Real
    #We need
    if size(trace,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Rdim will not work"))
    elseif size(trace,3) != size(rmaxes,1)
        throw(ErrorException("The number of rmaxes is not equal to the channels of the dataset"))
    else
        rdims = zeros(T, size(trace,3))
        dim_idx = zeros(Int64, size(trace,3))
        for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                rmax_val = rmax_lin .* rmaxes[ch]
                if rmax_val[1] > rmax_val[2]
                    rmax_val = reverse(rmax_val)
                end
                #rdim_thresh = rmaxes[ch] * 0.15
                
                if polarity < 0
                    minima = minimum(trace[swp, :, ch])
                else
                    minima = maximum(trace[swp, :, ch])
                end
                if rmax_val[1] < minima < rmax_val[2]
                    if minima < rdims[ch] && polarity < 0
                        rdims[ch] = minima
                        dim_idx[ch] = swp 
                    elseif minima > rdims[ch] && polarity > 0
                        rdims[ch] = minima
                        dim_idx[ch] = swp
                    end
                end
            end
        end
        if return_idx #In most cases, the rdim will be used to calculate the time to peak
            rdims |> vec, dim_idx |> vec
        else
            rdims |> vec
        end
    end
end

#This dispatch is for if there has been no rmax provided. 
dim_response(trace::Experiment; z = 0.0, rdim_percent = 0.15) = dim_response(trace, saturated_response(trace; z = z), rdim_percent = rdim_percent)

"""
This function calculates the time to peak using the dim response properties of the concatenated file
"""
function time_to_peak(trace::Experiment{T}, idxs::Array{Int64,1}) where T <: Real
    if size(trace,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Tpeak will not work"))
    elseif size(trace,3) != size(idxs,1)
        throw(ErrorException("The number of indexes is not equal to the channels of the dataset"))
    else
        t_peak = T[]
        for (ch, swp) in enumerate(idxs)
            if swp != 0
                t_series = trace.t[findall(trace.t .>= 0.0)]
                data = trace[swp, findall(trace.t .>= 0), ch]
                #println(argmin(data))
                push!(t_peak, t_series[argmin(data)])
            end
        end
        t_peak
    end
end

function get_response(trace::Experiment, rmaxes::Array{T,1}) where T <: Real
    responses = zeros(size(trace,1), size(trace,3))
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            minima = minimum(trace[swp, :, ch]) 
            responses[swp, ch] = minima < rmaxes[ch] ? rmaxes[ch] : minima
        end
    end
    responses
end

#Pepperburg analysis
"""
This function conducts a Pepperburg analysis on a single trace. 

    Two dispatches are available. 
    1) A rmax is provided, does not need to calculate rmaxes
    2) No rmax is provided, so one is calculated
"""
function pepperburg_analysis(trace::Experiment{T}, rmaxes::Array{T, 1}; recovery_percent = 0.60, kwargs...) where T <: Real
    if size(trace,1) == 1
        throw(error("Pepperburg will not work on single sweeps"))
    end
    r_rec = rmaxes .* recovery_percent
    #try doing this  different way
    t_dom = zeros(T, size(trace,1), size(trace,3))
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            not_recovered = findall(trace[swp, :, ch] .< r_rec[ch])
            if isempty(not_recovered)
                #println("Trace never exceeded $(recovery_percent*100)% the Rmax")
                t_dom[swp, ch] = NaN
            elseif isempty(trace.stim_protocol)
                #println("No stimulus protocol exists")
                t_dom[swp, ch] = trace.t[not_recovered[end]]
            else
                t_dom[swp, ch] = trace.t[not_recovered[end]] - trace.t[trace.stim_protocol[swp].index_range[1]]
            end
        end
    end
    t_dom
end

pepperburg_analysis(trace::Experiment{T}; kwargs...) where T <: Real= pepperburg_analysis(trace, saturated_response(trace; kwargs...); kwargs...)  


"""
The dominant time constant is calculated by fitting the normalized Rdim with the response recovery equation
"""
function recovery_tau(trace::Experiment{T}, dim_idx::Array{Int64,1}; τRec::T = 1.0) where T <: Real
    if size(trace,3) != length(dim_idx)
        throw(error("Size of dim indexes does not match channel size for trace"))
    else
        tauRec = T[]
        #This function uses the recovery model and takes t as a independent variable
        model(x,p) = map(t -> Recovery(t, p[1], p[2]), x)
        for ch in 1:size(trace,3)
            xdata = trace.t
            ydata = trace[dim_idx[ch], :, ch]
            xdata = xdata[argmin(ydata):end] .- xdata[argmin(ydata)]
            ydata = ydata[argmin(ydata):end]
            p0 = [xdata[1], τRec]
            fit = curve_fit(model, xdata, ydata, p0)
            println(sum(fit.resid))
            push!(tauRec, fit.param[2])
        end
        return tauRec
    end
end

"""
The integration time is fit by integrating the dim flash response and dividing it by the dim flash response amplitude
    - A key to note here is that the exact f(x) of the ERG trace is not completely known
    - The integral is therefore a defininte integral and a sum of the area under the curve
"""
function integration_time(trace::Experiment{T}, dim_idx::Array{Int64,1}) where T <: Real
    if size(trace,3) != length(dim_idx)
        throw(error("Size of dim indexes does not match channel size for trace"))
    else
        int_time = T[]
        for ch in 1:size(trace,3)
            dim_trace = trace[dim_idx[ch], :, ch]
            #The integral is calculated by taking the sum of all points (in μV) and dividing by the time range (in ms)
             #We have to make sure this response is in μV
            if trace.chUnits[ch] == "mV"
                rdim = minimum(dim_trace)*1000
                sum_data = sum(dim_trace.*1000)*(trace.dt*1000)
            else
                rdim = minimum(dim_trace)
                sum_data = sum(dim_trace)*trace.dt
            end
            push!(int_time, sum_data/rdim)
        end
        return int_time
    end
end