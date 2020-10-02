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
    #Normalization
    x_norm1, norm_factor1 = normalize(x_adj1);
    x_norm2, norm_factor2 = normalize(x_adj2);
    #CWT filtering (Probably not ready for CWT filtering )
    x_cwt1, cwt1_raster = cwt_filter(x_norm1, periods = 1:cutoff_octave);
    x_cwt2, cwt2_raster = cwt_filter(x_norm2, periods = 1:cutoff_octave);
    return x_cwt1.*norm_factor1, x_cwt2.*norm_factor2, x_stim
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


###############################These are the IR and Amplification models#############

"""
# Adult Intensity-Response models

## The relationship is demonstrated by 
\$R = f(I)\$ 

\$f(I) = R_{max}\\frac{I^n}{I^n_{1/2}+I^n}\$

if Response values are normalized to 1, then \$R_{max}\$ = 1 and can be cancelled out to form the equations

### Variables: 
- R: The response amplitude is the dependent variable
- I: The stimulus light intensity (I) is the independent variable
### Parameters: 
- R_max: Maximum saturating value(\$R_{max}\$)
- Ih: The flash strength required to elicit half of \$R_{max}\$: (\$I_{1/2}\$)
- n: The power of the equation
### Function usage
[IN 1]:  IR(I, Ih, n)

[OUT 1]: Response
"""
IR(I, Ih, n) = I^n / (Ih^n + I^n)

"""
# Developmental Intensity response (>P14)

## The relationship is demonstrated by 
\$R = f(I)\$ 
 where 
\$f(I) =R_{max}\\left(\\alpha(1 - e^{SI}) + (1-\\alpha)\\frac{I^n}\$

if Response values are normalized to 1, then \$R_{max}\$ = 1 and can be cancelled out to form the equations

### Variables: 
- R: The response amplitude is the dependent variable
- I: The stimulus light intensity (I) is the independent variable
### Parameters: 
- R_max: Maximum saturating value(\$R_{max}\$)
- Ih: The flash strength required to elicit half of \$R_{max}\$: (\$I_{1/2}\$)
- n: The power of the equation
- \$\\alpha\$: The temperature-dependent weighting coefficient:  
- S: he fractional sensitivity
### Function usage
[IN 1]:  IR_dev(I, Ih, n, α, SI, S)

[OUT 1]: Response_dev
"""
IR_dev(I, Ih, n, α, SI, S) = α*(1-exp(SI)) + (1-α)*(I^n / (Ih^n + S))

"""
# Amplification 

Amplification is a time series, therefore it is a function of time

## The relationship is demonstrated by
\$R = f(t)\$

\$\\\$f(t) = R_{max}(1-e^{-\\alpha(t-t_{eff})^2})\$\$

### Variables
- R: The response is the dependent variable
- t: Time is the independent variable.

### Parameters
- \$t_{eff}\$: The effective time delay is a short delay between stimulus onset and response onset indicative of the biomolecuar diffusion rates
- \$\\alpha\$: The amplification coefficient  represents the rate of the response increases from the biomolecular processes. 

### Function usage
[IN 1]:  AMP(t, α, t_eff, rmax)

[OUT 1]: Response

"""
AMP(t, α, t_eff, rmax) = t > t_eff ? rmax * (1 - exp(-α*(t-t_eff)^2)) : 0.0 