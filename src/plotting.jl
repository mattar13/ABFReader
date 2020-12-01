"""
This function plots by channel. This is the most basic functionality of the trace plot
"""
@recipe function f(nt::NeuroTrace; stim_plot = :subplot)
    grid := false
    layout := (size(nt,3), 1)
    if stim_plot == :include
        layout := (size(nt,3)-1, 1)
        for (i,ch) in enumerate(eachchannel(nt; include_stim = false))
            @series begin
                subplot := i
                x := nt.t
                y := ch
                yguide := "$(nt.chNames[i])($(nt.chUnits[i]))"
                if i == size(nt,3)
                    xguide := "Time ($(nt.tUnits))" 
                end
                ()
            end
            @series begin
                t_stim_start, t_stim_end = findstimRng(nt)
                subplot := i
                seriestype := :vline
                label := "Stimulus"
                y := [t_stim_end]
            end
        end
    else
        layout := (size(nt,3), 1)
        for (i,ch) in enumerate(eachchannel(nt))
            @series begin
                subplot := i
                x := nt.t
                y := ch
                yguide := "$(nt.chNames[i])($(nt.chUnits[i]))"
                if i == size(nt,3)
                    xguide := "Time ($(nt.tUnits))" 
                end
                ()
            end
        end
    end
end