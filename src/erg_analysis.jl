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

function RSQ(ŷ::Array{T}, y::Array{T}) where T <: Real
	ȳ = sum(ŷ)/length(ŷ)
	SSE = sum((y-ŷ).^2)
	SST = sum((y.-ȳ).^2)
	1-SSE/SST
end


"""
This function calculates the min, max, mean, and std of each data
"""
function calculate_basic_stats(data::Experiment{T}) where T <: Real
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
    - In ERG datas, a short nose component is usually present in saturated values
    - Does this same function work for the Rmax of nonsaturated responses?
    - Setting the saturated threshold to infinity will completely disregard the histogram method
"""
function saturated_response(data::Experiment{T}; 
          polarity::Int64 = -1, precision::Int64 = 500, z = 4
          family = false
     ) where T <: Real
     #We want to pick the region to analyze first
     if polarity < 0
          #Pick the local minima
          first_idxs = zeros(Int64, size(data,1), size(data,3))
          rmaxes = zeros(size(data,1), size(data,3))
          minima = argmin(data, dims = 2)
          for idx in minima
               first_idxs[idx[1], idx[3]] = idx[2]
          end
          for swp in 1:size(data,1), ch in 1:size(data,3)
               first_idx = first_idxs[swp, ch]
               x_data = data.t[first_idx:end] 
	          y_data = data.data_array[swp,first_idx:end,ch]
               mean = sum(y_data)/length(y_data)
	          deviation = z*std(y_data)
	          last_idx = findall(y_data .> mean)[1]
               x_data = x_data[1:last_idx] 
	          y_data = y_data[1:last_idx]
               d1 = diff(y_data)
	          z_diff = sum(d1)/length(d1) + 4*std(d1)
	          idxs = findall(d1.>z_diff)
               if !isempty(idxs)
                    #There is a nose component
                    bins = LinRange(mean, min(0.0, mean-deviation),  500)
                    h = Distributions.fit(Histogram, y_data[y_data.<mean], bins)
                    edges = collect(h.edges...)[2:end]
                    weights = h.weights./length(y_data)
                    rmax = edges[argmax(weights)]
               else
                    rmax = minimum(y_data)
               end
               rmaxes[swp, ch] = rmax
          end
          return rmaxes
     elseif polarity > 0
          #In this case we should just return the local maxima  
          return maximum(data, dims = 2)
     end
     #First we want to find out if the nose component exists. Otherwise we can just return the minima
end

"""
This function only works on concatenated files with more than one data
    Rmax argument should have the same number of sweeps and channels as the 
    In the rdim calculation, it is better to adjust the higher percent
    Example: no datas in 20-30% range, try 20-40%
"""
function dim_response(data::Experiment{T}, rmaxes::Array{T, 1}; return_idx = true, polarity::Int64 = -1, rmax_lin = [0.10, 0.50]) where T <: Real
    #We need
    if size(data,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Rdim will not work"))
    elseif size(data,3) != size(rmaxes,1)
        throw(ErrorException("The number of rmaxes is not equal to the channels of the dataset"))
    else
        rdims = zeros(T, size(data,3))
        #rdims = fill(NaN, size(data,3))
        dim_idx = zeros(Int64, size(data,3))
        for swp in 1:size(data,1)
            for ch in 1:size(data,3)
                rmax_val = rmax_lin .* rmaxes[ch]
                #println(rmax_val)
                if rmax_val[1] > rmax_val[2]
                    rmax_val = reverse(rmax_val)
                end
                #rdim_thresh = rmaxes[ch] * 0.15
                
                if polarity < 0
                    minima = minimum(data[swp, :, ch])
                else
                    minima = maximum(data[swp, :, ch])
                end
                #println(rmax_val[1])
                #println(minima)
                #println(rmax_val[1]< minima)

                #println(rmax_val[2])
                #println(minima)
                #println(minima <= rmax_val[2])
                if rmax_val[1] < minima < rmax_val[2]
                    #println("Minima in range")
                    if minima < rdims[ch] && polarity < 0
                        rdims[ch] = minima
                        dim_idx[ch] = swp 
                    elseif minima > rdims[ch] && polarity > 0
                        rdims[ch] = minima
                        dim_idx[ch] = swp
                    end
                else
                    #println("Minima not in range")
                    #rdims[ch] = NaN
                end
            end
        end
        rdims = map(x -> x == 0.0 ? NaN : x, rdims)
        #dim_idx = map(x -> x == 0.0 ? NaN : x, dim_idx)
        #println(rdims)
        if return_idx #In most cases, the rdim will be used to calculate the time to peak
            rdims |> vec, dim_idx |> vec
        else
            rdims |> vec
        end
    end
end

"""
This function calculates the time to peak using the dim response properties of the concatenated file
"""
function time_to_peak(data::Experiment{T}, dim_idx::Array{Int64,1}) where T <: Real
    if size(data,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Tpeak will not work"))
    elseif size(data,3) != size(dim_idx,1)
        throw(ErrorException("The number of indexes is not equal to the channels of the dataset"))
    else
        t_peak = T[]
        for (ch, swp) in enumerate(dim_idx)
            if swp != 0
                t_series = data.t[findall(data.t .>= 0.0)]
                temp_record = data[swp, findall(data.t .>= 0), ch]
                #println(argmin(data))
                push!(t_peak, t_series[argmin(temp_record)])
            else
                push!(t_peak, NaN)
            end
        end
        t_peak
    end
end

function get_response(data::Experiment{T}, rmaxes::Array{T,1}; polarity = -1) where T <: Real
    responses = zeros(size(data,1), size(data,3))
    for swp in 1:size(data,1)
        for ch in 1:size(data,3)
            if polarity == -1
                minima = minimum(data[swp, :, ch]) 
                responses[swp, ch] = minima < rmaxes[ch] ? rmaxes[ch] : minima
            elseif polarity == 1
                maxima = maximum(data[swp, :, ch]) 
                responses[swp, ch] = maxima >= rmaxes[ch] ? rmaxes[ch] : maxima
            end
        end
    end
    responses
end

#Pepperburg analysis
"""
This function conducts a Pepperburg analysis on a single data. 

    Two dispatches are available. 
    1) A rmax is provided, does not need to calculate rmaxes
    2) No rmax is provided, so one is calculated
"""
function pepperburg_analysis(data::Experiment{T}, rmaxes::Array{T, 1}; 
    recovery_percent = 0.60, kwargs...
    ) where T <: Real
    if size(data,1) == 1
        throw(error("Pepperburg will not work on single sweeps"))
    end
    r_rec = rmaxes .* recovery_percent
    #try doing this  different way
    t_dom = zeros(T, size(data,1), size(data,3))
    for swp in 1:size(data,1)
        for ch in 1:size(data,3)
            not_recovered = findall(data[swp, :, ch] .< r_rec[ch])
            if isempty(not_recovered)
                #println("data never exceeded $(recovery_percent*100)% the Rmax")
                t_dom[swp, ch] = NaN
            elseif isempty(data.stim_protocol)
                #println("No stimulus protocol exists")
                t_dom[swp, ch] = data.t[not_recovered[end]]
            else
                t_dom[swp, ch] = data.t[not_recovered[end]] - data.t[data.stim_protocol[swp].index_range[1]]
            end
        end
    end
    t_dom
end

pepperburg_analysis(data::Experiment{T}; kwargs...) where T <: Real= pepperburg_analysis(data, saturated_response(data; kwargs...); kwargs...)  



"""
The integration time is fit by integrating the dim flash response and dividing it by the dim flash response amplitude
- A key to note here is that the exact f(x) of the ERG data is not completely known
- The integral is therefore a defininte integral and a sum of the area under the curve
"""

function integration_time(data::Experiment{T}, dim_idx::Array{Int64,1}) where T <: Real
    if size(data,3) != length(dim_idx)
        throw(error("Size of dim indexes does not match channel size for data"))
    else
        int_time = T[]
        for ch in 1:size(data,3)
            if dim_idx[ch] == 0
                push!(int_time, NaN)
            else
                dim_data = data[dim_idx[ch], :, ch]
                #The integral is calculated by taking the sum of all points (in μV) and dividing by the time range (in ms)
                #We have to make sure this response is in μV
                if data.chUnits[ch] == "mV"
                    rdim = minimum(dim_data)*1000
                    sum_data = sum(dim_data.*1000)*(data.dt*1000)
                else
                    rdim = minimum(dim_data)
                    sum_data = sum(dim_data)*data.dt
                end
                push!(int_time, sum_data/rdim)
            end
        end
        return int_time
    end
end

# The below functions are created by fitting a model 

"""
The dominant time constant is calculated by fitting the normalized Rdim with the response recovery equation
"""
function recovery_tau(data::Experiment{T}, dim_idx::Array{Int64,1}; τRec::T = 1.0) where T <: Real
    if size(data,3) != length(dim_idx)
        throw(error("Size of dim indexes does not match channel size for data"))
    else
        fits = []
        gofs = []
        #This function uses the recovery model and takes t as a independent variable
        model(x,p) = map(t -> REC(t, -1.0, p[2]), x)
        for ch in 1:size(data,3)
            # println(dim_idx[ch])
            if dim_idx[ch] == 0
                push!(fits, fill(NaN, 2))
                push!(gofs, NaN)
            elseif dim_idx[ch] == NaN
                push!(fits, fill(NaN, 2))
                push!(gofs, NaN)
            else
                xdata = data.t
                ydata = data[dim_idx[ch], :, ch] 
                ydata ./= minimum(ydata) #Normalize the Rdim
                #cutoff all points below -0.5 and above -1.0
                begin_rng = findall(ydata .>= 1.0)[end]
                xdata = xdata[begin_rng:end]
                ydata = ydata[begin_rng:end]

                cutoff = findall(ydata .< 0.5)              
                if isempty(cutoff)
                    #println("Exception")
                    end_rng = length(ydata)
                else
                    end_rng = cutoff[1]
                end

                xdata = xdata[1:end_rng] .- xdata[1]
                ydata = -ydata[1:end_rng]
                p0 = [ydata[1], τRec]
                fit = curve_fit(model, xdata, ydata, p0)
                #report the goodness of fit
                SSE = sum(fit.resid.^2)
                ȳ = sum(model(xdata, fit.param))/length(xdata)
                SST = sum((ydata .- ȳ).^2)
                GOF = 1- SSE/SST
                push!(fits, fit.param)
                push!(gofs, GOF)
            end
        end
        return fits, gofs
    end
end


function amplification(data::Experiment{T}, rmaxes::Array{T,1}; 
        time_cutoff = 0.1,
        lb = [0.0, 0.001], ub = [Inf, 0.040]
    ) where T <: Real
    amp = zeros(2, size(data,1), size(data,3))
    gofs = zeros(T, size(data,1), size(data,3))
    for swp in 1:size(data,1), ch in 1:size(data,3)
        model(x, p) = map(t -> AMP(t, p[1], p[2], rmaxes[ch]), x)
        idx_end = findall(data.t .>= time_cutoff)[1]
        xdata = data.t[1:idx_end]
        ydata = data[swp,1:idx_end,ch]
        p0 = [200.0, 0.002]        
        fit = curve_fit(model, xdata, ydata, p0, lower = lb, upper = ub)
        #Check Goodness of fit
        SSE = sum(fit.resid.^2)
        ȳ = sum(model(xdata, fit.param))/length(xdata)
        SST = sum((ydata .- ȳ).^2)
        GOF = 1 - SSE/SST
        amp[1, swp, ch] = fit.param[1] #Alpha amp value
        amp[2, swp, ch] = fit.param[2] #Effective time value
        gofs[swp, ch] = GOF
    end
    return amp, gofs
end

#Lets get the file we want to use first
function IR_curve(data::Experiment{T}; 
        ih = 100.0, n::Real = 2.0,  lb = [0.0, 1.0, 0.0], ub = [Inf, 4.0, Inf],
        polarity::Int64 = -1, normalize = false, ignore_oversaturation = true,
        report_GOF = false
    ) where T <: Real
    if length(data.filename) > 1
        sensitivity = zeros(size(data,3))
        fit_rmaxs = zeros(size(data,3))
        fit_ns = zeros(size(data,3))
        #The file is not a concatenation in clampfit
        intensity = Float64[]
                
        #We can choose if the Rmax polarity is negative or positive
        if polarity == -1
            resp = -minimum(data, dims = 2)[:,1,:] 
        elseif polarity == 1
            resp = maximum(data, dims = 2)[:,1,:]
        end
        rmaxs = -saturated_response(data)
        #We need to make sure that the rmax doesn't include the nose component
        
        if ignore_oversaturation == true
            for i in 1:size(resp,2)
                over_saturated = resp[:,i] .> rmaxs[i] #Find out anything over the rmax
                under_saturated = resp[:,i] .<= rmaxs[i] #Find out anything under the rmax
                resp[:, i] = (rmaxs[i] .* over_saturated) + (resp[:, i] .* under_saturated)
            end
        end


        if normalize
            #If we choose to normalize, the rmaxes will become 1 and responses
            resp ./= maximum(resp)
            rmaxs = ones(size(rmaxs))
        end
        
        for (idx,info) in enumerate(data.filename)
            t_begin, t_end = data.stim_protocol[idx].timestamps
            t_stim = (t_end - t_begin)*1000
            file_info = formatted_split(info, format_bank)
            OD = Float64(file_info.ND) |> Transferrance
            Per_Int = Float64(file_info.Intensity)
            #println("$(file_info.ND) -> $(OD), $(Per_Int), $(t_stim)")
            photons = stimulus_model([OD, Per_Int, t_stim])
            push!(intensity, photons)
        end

        model_pars = [ih, n, sum(rmaxs)/2]
        for i in 1:size(data,3)
            if data.chUnits[i]== "mV" && normalize == false
                #The data is in mV we want μV
                resp[:,i] .*= 1000
                rmaxs[i] = rmaxs[i] * 1000
            end
            #println(rmaxs)
            model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)
            #Fit the Ih curve for each channel
            fit = curve_fit(model, intensity, resp[:,i], model_pars, lower = lb, upper = ub)
            println(fit.param)
            if report_GOF
                SSE = sum(fit.resid.^2)
                ȳ = sum(model(intensity, fit.param))/length(intensity)
                SST = sum((resp[:,i] .- ȳ).^2)
                GOF = 1- SSE/SST
                println("Goodness of fit: $GOF")
            end
            sensitivity[i] = fit.param[1]
            fit_ns[i] = fit.param[2]
            fit_rmaxs[i] = fit.param[3]

        end

        return (intensity, resp, sensitivity, fit_ns, fit_rmaxs)
    else
        #The file is a preconcatenated file
    end
end