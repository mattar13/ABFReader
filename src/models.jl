"""
This function is used to calculate the photon density based on the photon energy from the calibrator

The equation for this function is as follows

E = Photon Energy
C = speed of Light
E/(C)
"""
photons(E::Float64; λ::Int64 = 525) = ( E / (6.626e-34 * 3e8/(λ*10e-9)))*10e-8


"""
This function is used to calculate the transferrance from the optical density (D)
"""
Transferrance(D) = 10^-D

#For the next equation we have 3 fit variables previously determined


"""
This function is the relationship between: 
    Independent Variables: 
        Transferrance (T) 
        LED Intensity (I)
        Stimulus Time (t_stim)
    Dependent Variable
        Photons (P)
"""
stimulus_model(x::Array{Float64,1}, p::Array{Float64,1}) = x[1]*(p[1]*x[2]^2 + p[2]*x[2] + p[3])*x[3]
stimulus_model(x::Array{Float64,2}, p::Array{Float64,1}) = [excite_model(x[i,:], p) for i in 1:size(x,1)]
stimulus_model(x::Array{Float64,1}) = stimulus_model(x, [25352.59, 43857.01, 929.56]) 
stimulus_model(x::Array{Float64,2}) = stimulus_model(x, [25352.59, 43857.01, 929.56]) 
    
##############################These are the IR and Amplification models#############

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

\$f(I) =R_{max}\\left(\\alpha(1 - e^{SI}) + (1-\\alpha)\\frac{I^n}{Ih^n + S}\$

if Response values are normalized to 1, then \$R_{max}\$ = 1 and can be cancelled out to form the equations

### Variables: 
- R: The response amplitude is the dependent variable
- I: The stimulus light intensity (I) is the independent variable
### Parameters: 
- R_max: Maximum saturating value(\$R_{max}\$)
- Ih: The flash strength required to elicit half of \$R_{max}\$: (\$I_{1/2}\$)
- n: The power of the equation
- (\$\\alpha\$): The temperature-dependent weighting coefficient:  
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

\$f(t) = R_{max}(1-e^{-\\alpha(t-t_{eff})^2})\$

### Variables
- R: The response is the dependent variable
- t: Time is the independent variable.

### Parameters
- (\$t_{eff}\$): The effective time delay is a short delay between stimulus onset and response onset indicative of the biomolecuar diffusion rates
- (\$\\alpha\$): The amplification coefficient  represents the rate of the response increases from the biomolecular processes. 

### Function usage
[IN 1]:  AMP(t, α, t_eff, rmax)

[OUT 1]: Response

"""
AMP(t, α, t_eff, rmax) = t > t_eff ? rmax * (1 - exp(-α*(t-t_eff)^2)) : 0.0 

