"""
This function helps us to determine sweeps and channels in a layout for plotting
"""
function layout_helper(x, trace_size)
    if x == :sweeps
        return trace_size[1]
    elseif x == :channels
        return trace_size[3]
    elseif isa(x, Int64)
        return x
    end
end

function subplot_selector(x, trace_size)
    if x == :sweeps
        return 1:trace_size[1]
    elseif x == :channels
        return 1:trace_size[3]
    elseif isa(x, Int64)
        return [x]
    elseif isa(x, AbstractArray)
        return x
    end
end
"""
Plotting function. 

    to_plot: (sweep to plot, channel to plot) 
    This is a tuple of symbols, Int64 or arrays of Int64
        - using the keyword :sweeps will plot all sweeps
        - using the keyword :channels will plot all channels
        - using a number in index [1] will plot that specific sweep
        - using a number in index [2] will plot that specific channel
        - using an array in index [1] will plot all of the sweeps in the index
        - using an array in index [2] will plot all of the channels in the index
    layout: (rows, columns)
        - using the keyword :channels plots the number of channels in the respective column/region
        - using the keyword :sweeps plots the respective number of sweeps in the respective column/region
    plot_stim_mode:
    - :overlay_vline_start -> Beginning of stimuli are plotted as a line
    - :overlay_vline_end -> End of stimuli is plotted as a line
    - :overlay_vspan -> 
"""
@recipe function f(
        exp::Experiment{T}; 
        to_plot = (:sweeps, :channels), 
        layout = (:channels, 1),
        plot_stim_mode = :none, #We will set this as default none for now
        label = "", label_stim = false,
        kwargs...
    ) where T <: Real
    
    #Set the basic characteristics of each plot
    grid := false
    if haskey(kwargs, :c) 
        color := kwargs[:c]
    elseif haskey(kwargs, :color)
        color := kwargs[:color]
    else
        color := :black
    end

    layout := map(lay -> layout_helper(lay, size(exp)), layout)
    
    swp_rng, ch_rng = map(subp -> subplot_selector(subp, size(exp)), to_plot)
    #println(ch_rng)
    for swp in swp_rng, ch in ch_rng
        if size(exp,3) == 1
            xlabels = "Time ($(exp.tUnits))"
        else
            xlabels = reshape(repeat([""], size(exp,3)-1), (1, size(exp,3)-1))
            xlabels[end] = "Time ($(exp.tUnits))"
        end
        xguide := xlabels
        @series begin
            label := label
            subplot := ch
            x := exp.t #t_series
            y := exp[swp, :, ch]
            yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
            ()
        end
        if !isempty(exp.stim_protocol) && plot_stim_mode != :none
            @series begin
                subplot := ch
                seriescolor := :yellow
                
                if label_stim && swp == 1
                    label := "stimulus"
                else
                    label := ""
                end
                if plot_stim_mode == :overlay_vline_start
                    linewidth := 2.0
                    seriestype := :vline
                    y := [exp.stim_protocol[swp].timestamps[1]]
                    yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
                    ()
                elseif plot_stim_mode == :overlay_vline_end
                    linewidth := 2.0
                    seriestype := :vline
                    y := [exp.stim_protocol[swp].timestamps[2]]
                    yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
                    ()
                elseif plot_stim_mode == :overlay_vspan
                    linewidth := 1.0
                    linecolor := :yellow
                    seriesalpha := 0.5
                    seriestype := vspan
                    y := [exp.stim_protocol[swp].timestamps...]
                    yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
                    ()
                end
            end
        end
    end
end