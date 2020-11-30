#%% Test to see if the importing works
using Revise
using NeuroPhys
using Plots
println("Package properly exported")

#%% Test the exporting and filtering of .abf files
target_path = "test\\to_filter.abf"
#%%
trace = extract_abf(target_path); #Extract the data
println(trace.data_array |> size)
#%%
drift_trace = baseline_cancel(trace; mode = :slope, region = :prestim) #Cancel drift
baseline_trace = baseline_cancel(drift_trace; mode = :mean, region = :prestim) #Baseline data
filter_trace = lowpass_filter(baseline_trace) #Lowpass filter using a 40hz 8-pole filter
cwt_trace = cwt_filter(baseline_trace) #Use a continuous wavelet transform to remove noise, but keep time series info
#%%
@recipe function f(nt::NeuroTrace; plotby = :channel, display_stim = :subplot)
    grid := false
    if plotby == :channel
        if display_stim == :include
            n_sub_x = size(nt,3)-length(nt.stim_ch)
            n_sub_y = 1
        elseif display_stim == :subplot
            n_sub_x = size(nt,3)
            n_sub_y = 1
        end
        layout := (n_sub_x, n_sub_y)
        for (i,ch) in enumerate(eachchannel(nt))
            if display_stim == :include
                if i == nt.stim_ch
                    nothing
                else
                    println(i)
                    @series begin
                        subplot := i
                        x := nt.t
                        y := ch
                        yguide := "$(nt.chNames[i])($(nt.chUnits[i]))"
                        if i == n_sub_x
                            xguide := "Time ($(nt.tUnits))" 
                        end
                        ()
                    end
                    @series begin
                        subplot := i
                        series := :vline
                        x := [nt.t[findstimRng(nt)[end]]]
                        y := [nt.t[findstimRng(nt)[end]]]
                        ()
                    end
                end
            end
        end
    elseif plotby == :sweep
        layout := (size(nt,1), 1)
        for (i,swp) in enumerate(eachsweep(nt))
            @series begin
                subplot := i
                x := nt.t
                y := swp
                yguide := "Sweep $i"
                vline := [findstimRng(nt)[end]]
                if i == size(nt,1)
                    xguide := "Time ($(nt.tUnits))" 
                end
                ()
            end
        end
    elseif plotby == :both
        layout := (size(nt,3), size(nt,1))
        for i = 1:size(nt,1)
            for j = 1:size(nt,3)
                @series begin
                    subplot := i*j
                    x := nt.t
                    y := nt.data_array[i,:,j]
                    yguide := "$(nt.chNames[j])($(nt.chUnits[j])) Swp $i"
                    ()
                end
            end
        end
    end
end
#%%
plot(trace, plotby = :channel, display_stim = :include, c = :blue)
#%% Test the analysis of .abf files