"""
This function plots by channel. This is the most basic functionality of the trace plot
"""
@recipe function f(nt::NeuroTrace; stim_plot = :include, time_adjusted = true)
    grid := false
    if nt.stim_ch == -1
        #If there is no stimulus, always use the subplot
        stim_plot = :subplot
    end

    if stim_plot == :include
        layout := (size(nt,3)-1, 1)
        for swp in 1:size(nt,1)
            for ch in 1:size(nt,3)-1
                xlabels = reshape(repeat([""], size(nt,3)-1), (1, size(nt,3)-1))
                xlabels[end] = "Time ($(nt.tUnits))"
                xguide := xlabels
                @series begin
                    subplot := ch
                    if time_adjusted 
                        x := nt.t .- nt.t[findstimRng(nt)[end]]
                    else
                        x := nt.t
                    end
                    y := nt[swp, :, ch]
                    yguide := "$(nt.chNames[ch])($(nt.chUnits[ch]))"
                    ()
                end
                @series begin
                    t_stim_start, t_stim_end = findstimRng(nt)
                    subplot := ch
                    seriescolor := :yellow
                    linewidth := 4.0
                    seriestype := :vline
                    label := "Stimulus"
                    y := [nt.t[t_stim_end]]
                    ()
                end
            end
        end
    else
        layout := (size(nt,3), 1)
        for swp in 1:size(nt,1)
            for ch in 1:size(nt,3)
                xlabels = reshape(repeat([""], size(nt,3)), (1, size(nt,3)))
                xlabels[end] = "Time ($(nt.tUnits))"
                xguide := "Time ($(nt.tUnits))" 
                @series begin
                    subplot := ch
                    if time_adjusted 
                        x := nt.t .- nt.t[findstimRng(nt)[end]]
                    else
                        x := nt.t
                    end
                    y := nt[swp, :, ch]
                    yguide := "$(nt.chNames[ch])($(nt.chUnits[ch]))"
                    ()
                end
            end
        end
    end
end