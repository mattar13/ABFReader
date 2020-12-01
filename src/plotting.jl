"""
This function plots by channel. This is the most basic functionality of the trace plot
"""
@recipe function f(nt::NeuroTrace; stim_plot = :subplot)
    grid := false    
    if stim_plot == :include
        layout := (size(nt,3)-1, 1)
        for (i,ch) in enumerate(eachchannel(nt; include_stim = false))
            xlabels = reshape(repeat([""], size(nt,3)-1), (1, size(nt,3)-1))
            xlabels[end] = "Time ($(nt.tUnits))"
            xguide := xlabels
            @series begin
                subplot := i
                x := nt.t
                y := ch
                yguide := "$(nt.chNames[i])($(nt.chUnits[i]))"
                ()
            end
            @series begin
                t_stim_start, t_stim_end = findstimRng(nt)
                subplot := i
                seriescolor := :yellow
                lw := 4.0
                seriestype := :vline
                label := "Stimulus"
                y := [nt.t[t_stim_end]]
                ()
            end
        end
    else
        layout := (size(nt,3), 1)
        for (i,ch) in enumerate(eachchannel(nt))
            xlabels = reshape(repeat([""], size(nt,3)), (1, size(nt,3)))
            xlabels[end] = "Time ($(nt.tUnits))"
            xguide := "Time ($(nt.tUnits))" 
            @series begin
                subplot := i
                x := nt.t
                y := ch
                yguide := "$(nt.chNames[i])($(nt.chUnits[i]))"
                ()
            end
        end
    end
end