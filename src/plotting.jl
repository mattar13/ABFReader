"""
Plotting function. 

    sweep_placement = 
    - :single -> All sweeps plotted on a single plot
    - :individual -> all sweeps plotted on a individual plot
    plot_stim_mode:
    - :overlay_vline_start -> Beginning of stimuli are plotted as a line
    - :overlay_vline_end -> End of stimuli is plotted as a line
    - :overlay_vspan -> 
"""
@recipe function f(exp::Experiment; sweep_placement = :single, plot_stim_mode = :overlay_vline_start)
    grid := false
    if sweep_placement == :single 
        #Place all the sweeps in a single plot
        layout := (size(exp,3), 1)
    elseif sweep_placement == :individual
        println("Not quite ready yet")
        layout := (size(exp,1), size(exp,3))
    end

    for swp in 1:size(exp,1), ch in 1:size(exp,3)
        xlabels = reshape(repeat([""], size(exp,3)-1), (1, size(exp,3)-1))
        xlabels[end] = "Time ($(exp.tUnits))"
        xguide := xlabels
        @series begin
            label := ""
            subplot := ch
            x := exp.t #t_series
            y := exp[swp, :, ch]
            yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
            ()
        end
        @series begin
            subplot := ch
            seriescolor := :yellow
            linewidth := 2.0
            label := ""
            if plot_stim_mode == :overlay_vline_start
                seriestype := :vline
                y := [exp.stim_protocol[swp].timestamps[1]]
                yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
                ()
            elseif plot_stim_mode == :overlay_vline_end
                seriestype := :vline
                y := [exp.stim_protocol[swp].timestamps[2]]
                yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
                ()
            elseif plot_stim_mode == :overlay_vspan
                alpha := 0.2
                seriestype := vspan
                y := [exp.stim_protocol[swp].timestamps...]
                yguide := "$(exp.chNames[ch])($(exp.chUnits[ch]))"
                ()
            end
        end
    end
end