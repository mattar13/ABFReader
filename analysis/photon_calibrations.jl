### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# ╔═╡ 41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
using NeuroPhys, DataFrames

# ╔═╡ 19b55080-1ee1-11eb-0b48-7f7133b81e6b
import NeuroPhys: stim_intensity, photons

# ╔═╡ 4572cda0-1ece-11eb-33da-212c1c12cb49
#In order to make this script work, this file is needed
target_folder = "D:\\Data\\Calibrations\\525Green_Curve"

# ╔═╡ 03afd290-1ecf-11eb-2738-b5c190ada7a4
stim_files = target_folder |> parse_abf

# ╔═╡ 26d1333e-1ecf-11eb-0400-2507bd73ad21
split(stim_files[1], "\\")[end] |> filename_extractor

# ╔═╡ 7cbada1e-1ee0-11eb-39bb-5d77235312ae
begin
	data_photon = DataFrame(
		ND = Int[], Percent = Int64[], Stim_time = Int64[], Photons = Float64[]
	)
	for file in stim_files
		info = split(file, "\\")[end] |> filename_extractor
		t, data_array, dt = extract_abf(file; chs = -1);
		stim_t, stim_y = stim_intensity(file; chs = -1);
		photon_d = stim_y .|> photons
		println(photon_d[1])
	end
end

# ╔═╡ Cell order:
# ╠═41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
# ╠═19b55080-1ee1-11eb-0b48-7f7133b81e6b
# ╠═4572cda0-1ece-11eb-33da-212c1c12cb49
# ╟─03afd290-1ecf-11eb-2738-b5c190ada7a4
# ╠═26d1333e-1ecf-11eb-0400-2507bd73ad21
# ╠═7cbada1e-1ee0-11eb-39bb-5d77235312ae
