"""
This function plots by channel. This is the most basic functionality of the trace plot
"""
@recipe function f(nt::NeuroTrace; stim_plot = :subplot)
    grid := false    
    if stim_plot == :include
        layout := (size(nt,3)-1, 1)
        for (i,ch) in enumerate(eachchannel(nt; include_stim = false))
            xlabels = repeat(["";], size(nt,3)-(1+1))
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
                seriescolor = :yellow
                seriestype := :vline
                label := "Stimulus"
                y := [nt.t[t_stim_end]]
                #yguide := "$(nt.chNames[i])($(nt.chUnits[i]))"
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