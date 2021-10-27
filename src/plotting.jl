"""
This function helps us to determine sweeps and channels in a layout for plotting
"""
function layout_helper(x::Symbol, trace_size)
    if x == :channels
        return trace_size[1]
    elseif x == :sweeps
        return trace_size[2]
    end
end
layout_helper(x::Int64, trace_size) = x

#These are all just convienance functions to help select the subplot
subplot_selector(x::Int64, trace_size) = [x]
subplot_selector(x::AbstractArray, trace_size) = x
subplot_selector(x::UnitRange{T}, trace_size) where T<: Real = x

function subplot_selector(x::Symbol, trace_size)
    if x == :sweeps
        return 1:trace_size[1]
    elseif x == :channels
        return 1:trace_size[3]
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
        to_plot = (:channels, :sweeps), #(row, column)
        subplot = -1, 
        layout = (:channels, 1), 
        plot_stim_mode = :none, #We will set this as default none for now
        label = "", label_stim = false,
        xlabels = nothing, ylabels = nothing
    ) where T <: Real
    
    #Set the basic characteristics of each plot
    grid := false

    
    #swp_rng, ch_rng
    plt_rows, plt_cols = map(subp -> subplot_selector(subp, size(exp)), to_plot)
    lay = (map(lay -> layout_helper(lay, (plt_rows|>length, plt_cols|>length)), layout))
    #println(lay)
    layout := lay


    for (subp_row, row) in enumerate(plt_rows), (subp_col, col) in enumerate(plt_cols)
        #println("Row: $row")
        #println("Col: $col")
        if layout[1] == :channels || layout[2] == :sweeps       
            swp = col
            ch = row
        elseif layout[1] == :sweeps || layout[2] == :channels
            swp = col
            ch = row
        end
        
        subp_row = row > lay[1] ? lay[1] : subp_row
        subp_col = col > lay[2] ? lay[2] : subp_col
        subp = subp_col * subp_row

        #println(subp)
        if label != "" && swp == 1
            label := label
        end
        if size(exp,3) == 1 && isnothing(xlabels)
            xlabels = "Time (sec)"
        elseif isnothing(xlabels)
            xlabels = reshape(repeat([""], size(exp,3)-1), (1, size(exp,3)-1))
            xlabels[end] = "Time (sec)"
        else
            xlabels = xlabels
        end

        xguide := xlabels
        @series begin
            label := label
            if subplot != -1
                subplot := subplot
            else
                subplot := subp
            end
            x := exp.t #t_series
            y := exp[swp, :, ch]
            yguide := isnothing(ylabels) ? "$(exp.chNames[ch])($(exp.chUnits[ch]))" : ylabels
            ()
        end
        if !isempty(exp.stim_protocol) && plot_stim_mode != :none
            @series begin
                if subplot != -1
                    subplot := subplot
                else
                    subplot := subp
                end
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