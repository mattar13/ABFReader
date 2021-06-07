root = "E:\\Data\\ERG\\Retinoschisis\\"
test_paths = root |> parse_abf
	#lets just summarize all files we have this far
	summary = DataFrame(
		:Path => test_paths, 
		:Year => 0, :Month => 0, :Date => 0,
		:Animal => 0, :Condition => "Nothing", :Wavelength => 525, 
		:Age => 0, :Genotype => "None", :Channels => 0
		#:Min => [0.0], :Mean => [0.0], :Max => [0.0]
	)
for (idx, path) in enumerate(test_paths)
    println(joinpath(splitpath(path)[1:6]...))
    println(path)
    println("file $idx")
    nt = formatted_split(path, format_bank)
    println(nt)
    summary[idx, :Year] = nt.Year
    summary[idx, :Month] = nt.Month
    summary[idx, :Date] = nt.Date
    summary[idx, :Animal] = nt.Animal
    summary[idx, :Condition] = nt.Condition
    summary[idx, :Wavelength] = nt.Wavelength
    summary[idx, :Age] = nt.Age
    #data = extract_abf(path)
    #summary[idx, :Channels] = length(data.chNames)
    summary[idx, :Genotype] = "$(nt.Genotype)"
end

q = summary |> 
    @unique({_.Year, _.Month, _.Date, _.Animal})|>
    @map({
        Root = joinpath(splitpath(_.Path)[1:6]...),
        _.Year, _.Month, _.Date, _.Animal, _.Age, _.Genotype
    }) |>
DataFrame


for row in eachrow(q)
    titlen = "$(row.Year)_$(row.Month)_$(row.Date)_$(row.Animal)_$(row.Age)_$(row.Genotype)"
    filepath = joinpath(root, "$(titlen).png")
    println(titlen)
    if !isfile(filepath)
        data_test = extract_abf(row.Root)
        truncate_data!(data_test, t_pre = 0.2, t_post = 0.5);
        baseline_cancel!(data_test, mode = :slope); 
        p_valid = plot(data_test, title = titlen)#, to_plot = (15, :channels))
        savefig(p_valid, joinpath(root, "$(titlen).png"))
    end
end
