"""
This function plots by channel. This is the most basic functionality of the trace plot
"""
@recipe function f(nt::Experiment; stim_plot = :include, time_adjusted = true)
    grid := false
    if nt.stim_ch == -1
        #If there is no stimulus, always use the subplot
        stim_plot = :subplot
    end

    if stim_plot == :include
        layout := (size(nt,3)-1, 1)

        for swp in 1:size(nt,1)
            for ch in 1:size(nt,3)-1

                if time_adjusted 
                    t_series = nt.t .- nt.t[findstimRng(nt)[swp, 1]]
                else
                    t_series = nt.t
                end
                xlabels = reshape(repeat([""], size(nt,3)-1), (1, size(nt,3)-1))
                xlabels[end] = "Time ($(nt.tUnits))"
                xguide := xlabels
                @series begin
                    subplot := ch
                    x := nt.t #t_series
                    y := nt[swp, :, ch]
                    yguide := "$(nt.chNames[ch])($(nt.chUnits[ch]))"
                    ()
                end
                @series begin
                    subplot := ch
                    seriescolor := :yellow
                    linewidth := 2.0
                    seriestype := :vline
                    if swp == 1
                        label := "Stimulus"
                    end
                    y := [0.0]
                    yguide := "$(nt.chNames[ch])($(nt.chUnits[ch]))"
                    ()
                end
            end
        end
    else
        layout := (size(nt,3), 1)
        for swp in 1:size(nt,1)
            for ch in 1:size(nt,3)
                if time_adjusted 
                    t_series = nt.t .- nt.t[findstimRng(nt)[swp, 1]]
                else
                    t_series = nt.t
                end
                xlabels = reshape(repeat([""], size(nt,3)), (1, size(nt,3)))
                xlabels[end] = "Time ($(nt.tUnits))"
                xguide := "Time ($(nt.tUnits))" 
                @series begin
                    subplot := ch
                    x := nt.t #t_series
                    y := nt[swp, :, ch]
                    yguide := "$(nt.chNames[ch])($(nt.chUnits[ch]))"
                    ()
                end
            end
        end
    end
end