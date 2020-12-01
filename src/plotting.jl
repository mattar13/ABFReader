@recipe function plot_by_channel(nt::NeuroTrace; stim_plot = :sublot)
    println("Functional interface working")
    if stim_plot == :subplot
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
    else
        println("Not implemented yet")    
    end
end

#These are plotting recipes. The mose basic plots the 
@recipe function f(nt::NeuroTrace; plotby = :channel, stim_plot = :subplot)
    grid := false
    if plotby == :channel
        plot_by_channel(nt; stim_plot = stim_plot)
    elseif plotby == :sweep
        layout := (size(nt,1), 1)
        for (i,swp) in enumerate(eachsweep(nt))
            @series begin
                subplot := i
                x := nt.t
                y := swp
                yguide := "Sweep $i"
                if i == size(nt,1)
                    xguide := "Time ($(nt.tUnits))" 
                end
                ()
            end
        end
    elseif plotby == :both
        layout := (size(nt,1), size(nt,3))
        for (i,ch) in enumerate(eachchannel(nt))
            @series begin
                subplot := i
                x := nt.t
                y := nt.data_array[i,:,1]
                yguide := "$(nt.chNames[i])($(nt.chUnits[i]))"
                ()
            end
        end
    end
end