using ClusterDepth
using Random
using CairoMakie
using UnfoldSim
using Unfold
using UnfoldMakie
using Statistics

# ## How to use clusterDepth multiple comparison correction on multichannel data

# This tutorial is adapted from the first EEG example and uses the HArtMuT NYhead model (https://github.com/harmening/HArtMuT) to simulate multiple channels. 

# First set up the EEG simulation as before, with one subject and 40 trials:
design = SingleSubjectDesign(conditions=Dict(:condition => ["car", "face"])) |> x -> RepeatDesign(x, 40);
p1 = LinearModelComponent(;
    basis=p100(; sfreq=250),
    formula=@formula(0 ~ 1),
    β=[1.0]
);

n170 = LinearModelComponent(;
    basis=UnfoldSim.n170(; sfreq=250),
    formula=@formula(0 ~ 1 + condition),
    β=[1.0, 0.5], # condition effect - faces are more negative than cars
);
p300 = LinearModelComponent(;
    basis=UnfoldSim.p300(; sfreq=250),
    formula=@formula(0 ~ 1 + condition),
    β=[1.0, 0], # no p300 condition effect
);

# Now choose some source coordinates for each of the p100, n170, p300 that we want to simulate, and use the helper function closest_srcs to get the HArtMuT sources that are closest to these coordinates:
src_coords = [
    [20, -78, -10], #p100
    [-20, -78, -10], #p100
    [50, -40, -25], #n170
    [0, -50, 40], #p300
    [0, 5, 20], #p300
];

headmodel_HArtMuT = headmodel()
get_closest = coord -> UnfoldSim.closest_src(coord, headmodel_HArtMuT.cortical["pos"]) |> pi -> magnitude(headmodel_HArtMuT; type="perpendicular")[:, pi]

p1_l = p1 |> c -> MultichannelComponent(c, get_closest([-20, -78, -10]))
p1_r = p1 |> c -> MultichannelComponent(c, get_closest([20, -78, -10]))
n170_r = n170 |> c -> MultichannelComponent(c, get_closest([50, -40, -25]))
p300_do = p300 |> c -> MultichannelComponent(c, get_closest([0, -50, -40]))
p300_up = p300 |> c -> MultichannelComponent(c, get_closest([0, 5, 20]))

data, events = simulate(MersenneTwister(1), design, [p1_l, p1_r, n170_r, p300_do, p300_up],
    UniformOnset(; offset=0.5 * 250, width=100),
    RedNoise(noiselevel=1); return_epoched=true);


# ## Plotting
# This is what the data looks like, for one channel/trial respectively:
f = Figure()
Axis(f[1, 1], title="Single channel, all trials", xlabel="time", ylabel="y")
series!(data[1, :, :]', solid_color=:black)
lines!(mean(data[1, :, :], dims=2)[:, 1], color=:red)

hlines!([0], color=:gray)
Axis(f[2, 1], title="All channels, average over trials", xlabel="time", ylabel="y")
series!(mean(data, dims=3)[:, :, 1], solid_color=:black)
hlines!([0], color=:gray)
f

# And some topoplots: 
positions = [Point2f(p[1] + 0.5, p[2] + 0.5) for p in to_positions(headmodel_HArtMuT.electrodes["pos"]')]

df = UnfoldMakie.eeg_matrix_to_dataframe(mean(data, dims=3)[:, :, 1], string.(1:length(positions)));
Δbin = 20 # 20 samples / bin
plot_topoplotseries(df, Δbin; positions=positions, visual=(; enlarge=1, label_scatter=false))

# ## ClusterDepth
# Now that the simulation is done, let's try out ClusterDepth and plot our results
# Note that this is a simple test of "activity" vs. 0
pvals = clusterdepth(data; τ=1.6, nperm=100);
fig, ax, hm = heatmap(transpose(pvals))
ax.title = "pvals";
ax.xlabel = "time";
ax.ylabel = "channel";
Colorbar(fig[:, end+1], hm);
fig
