### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ acc2d270-3427-11eb-37bf-33a3be171271
using PlutoUI

# ╔═╡ ceb9bce0-3427-11eb-3877-a779d14f63fc
using NeuroPhys

# ╔═╡ cebbb8b0-3427-11eb-04d5-7bf211ecd1b0
using DataFrames, Query

# ╔═╡ cebf8940-3427-11eb-0b10-210fe9e4ce7e
using Plots, StatsPlots, LsqFit

# ╔═╡ cec10fe0-3427-11eb-20a0-3b01c33d4284
using Statistics, Distributions

# ╔═╡ 1cd255ee-2ced-11eb-3256-b90aec063b67
md"
# Analysis of single concatenated files
"

# ╔═╡ d194c94e-3427-11eb-20d0-33898e117d26
target_folder = "D:\\Data\\ERG\\Data from paul\\"

# ╔═╡ d62d47be-3428-11eb-0b01-dbcfc2a257ac
paths = (target_folder |> parse_abf)

# ╔═╡ 0c91b030-3429-11eb-0eb1-7ffa6013aff4
begin
	data = extract_abf(paths[1]; swps = -1, chs = -1)
	truncate_data!(data)
end

# ╔═╡ 33290852-342a-11eb-1307-c992bb1a6df5
minimum(data.data_array, dims = 2)

# ╔═╡ Cell order:
# ╠═1cd255ee-2ced-11eb-3256-b90aec063b67
# ╠═acc2d270-3427-11eb-37bf-33a3be171271
# ╠═ceb9bce0-3427-11eb-3877-a779d14f63fc
# ╠═cebbb8b0-3427-11eb-04d5-7bf211ecd1b0
# ╠═cebf8940-3427-11eb-0b10-210fe9e4ce7e
# ╠═cec10fe0-3427-11eb-20a0-3b01c33d4284
# ╠═d194c94e-3427-11eb-20d0-33898e117d26
# ╠═d62d47be-3428-11eb-0b01-dbcfc2a257ac
# ╠═0c91b030-3429-11eb-0eb1-7ffa6013aff4
# ╠═33290852-342a-11eb-1307-c992bb1a6df5
