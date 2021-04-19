using NeuroPhys
using DataFrames, XLSX, Query
#%% 1) Plotting and analyzing the IR analysis from Cole and Dustin
#Open the data frame first 
IR_datafile = "E:\\Data\\ERG\\Gnat\\IR_analysis.xlsx"
df = DataFrame(XLSX.readtable(IR_datafile, "Pauls_IR_analysis")...)

#%%
experiments = df |> 
    @unique({_.Age, _.Genotype, _.Wavelength, _.Drugs, _.Photoreceptors, _.Rearing}) |> 
    #@filter(_.Photons > 0.0) |> 
    DataFrame

experiments
#%%

df
#%%
for (i, row) in enumerate(eachrow(experiments))
    #println(row)
    println("Summarizing data from $i / $(length(eachrow(experiments)))")
    Qi = df |>
        @filter(_.Genotype == row.Genotype) |>
        @filter(_.Age == row.Age) |>
        @filter(_.Photoreceptors == row.Photoreceptors) |> 
        @filter(_.Drugs != "b-waves") |> 
        @filter(_.Rearing == row.Rearing) |>
        @filter(_.Wavelength == row.Wavelength) |> 
        @filter(_.Photons != 0.0) |>
        @map({
                _.Path, _.Age, _.Wavelength, _.Photoreceptors, _.Rearing,
                _.Photons, _.Response, _.Rmax
            }) |> 
        DataFrame
    if any(Qi.Photons .> 0.0)

        model(x, p) = map(I -> IR(I, p[1], p[2])*p[3], x)
        println(Qi)
        model_pars = [1.0, 2.0, Qi.Rmax[1]]
        lb = [0.0, 0.2, 0.0]
        ub = [Inf, 4.0, Inf]
        xdata = Qi.Photons
        ydata = Qi.Response
        fit = curve_fit(model, xdata, ydata, model_pars, lower = lb, upper = ub)
        println(fit.param)
        SSE = sum(fit.resid.^2)
        ȳ = sum(model(xdata, fit.param))/length(xdata)
        SST = sum((ydata .- ȳ).^2)
        R2 = 1- SSE/SST
        println("Goodness of fit: $R2")
        fig_i = plot(xdata, ydata, 
            st = :scatter, xaxis = :log,
            xlabel = "Intensity (Photons/μM^2)", ylabel = "Response (μV)",
            label = "", grid = false
            )
        plot!(fig_i, x -> model(x, fit.param), LinRange(minimum(xdata), maximum(xdata), 10000), 
            label = "Fit", title = "P$(row.Age) $(row.Wavelength)nm $(row.Photoreceptors) $(row.Rearing)"
        
        )
        savefig("$(row.Age)_$(row.Wavelength)_$(row.Photoreceptors).png")
    end
end



#%% 2) Figure out the new naming convention used by Paul
-