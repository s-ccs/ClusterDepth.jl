using ClusterDepth
using Random
using LinearAlgebra
using CairoMakie
using UnfoldSim
using StatsBase
using Distributions
using DataFrames
using Unfold
using UnfoldMakie
using MAT
using CSV
using CoordinateTransformations
using StaticArrays
using TopoPlots

# ## How to use clusterDepth multiple comparison correction on multichannel data

# This tutorial is adapted from the first EEG example and uses the HArtMuT NYhead model (https://github.com/harmening/HArtMuT) to simulate multiple channels. 

# First set up the EEG simulation as before, with one subject and 40 trials:
design = SingleSubjectDesign(conditions= Dict(:condition=>["car","face"]))|>x->RepeatDesign(x,40);	
p100 = LinearModelComponent(;
	basis=UnfoldSim.p100(;sfreq=250),
	formula=@formula(0 ~ 1),
	β=[1.]
);
n170 = LinearModelComponent(;
	basis=UnfoldSim.n170(;sfreq=250),
	formula=@formula(0 ~ 1+condition),
	β=[1.,-0.5], # condition effect - faces are more negative than cars
);
p300 = LinearModelComponent(;
	basis=UnfoldSim.p300(;sfreq=250),
	formula=@formula(0 ~ 1+condition),
	β=[1.,0], # no p300 condition effect
);

# Next, simulate data for P100 (two sources), N170, and P300 (two sources) - single channel, multiple trials.
dataList = [];
for sig = [p100,p100,n170,p300,p300]
	data, events = simulate(MersenneTwister(1), design, sig, UniformOnset(; offset=5, width=4), RedNoise(noiselevel=0.2); return_epoched=true);
	push!(dataList,data);
end

# Import the HArtMuT NYhead model
dlpath = download("https://github.com/harmening/HArtMuT/raw/main/HArtMuTmodels/HArtMuT_NYhead_small.mat");
file = matopen(dlpath);
hartmut = read(file, "HArtMuT");
close(file);
model = hartmut["cortexmodel"];

# and calculate the magnitude for the leadfield.
s = size(model["leadfield"]);
magnitude = zeros(s[1:2]);
for i=1:s[1] 
	for j=1:s[2]
		magnitude[i,j] = norm(model["leadfield"][i,j,:]);
	end
end

# Now choose some source coordinates for each of the p100, n170, p300 that we want to simulate, and use the helper function closest_srcs to get the HArtMuT sources that are closest to these coordinates:
src_coords = [
	[20 -78 -10], #p100
	[-20 -78 -10], #p100
	[50 -40 -25], #n170
	[0 -50 40], #p300
	[0 5 20], #p300
];
closest_src_indices = closest_srcs(src_coords, model["pos"]);

# Get multichannel data by multiplying the simulated one-channel data with the factors we get from the leadfield matrix - repeat this for each of the 5 sources
l =[];
for k = 1:length(dataList)
		push!(l,magnitude[:,closest_src_indices[k]] .* reshape(dataList[k],1,size(dataList[k])...));
end

# Add up the data from both sources for p1 and p3 respectively, then add along time (dimension 2) to get final data: channels x time x trials
p1_data = l[1] + l[2];
n1_data = l[3];
p3_data = l[4] + l[5];
finalData = zeros(size(l[1],1),maximum(size.(l,2)),size(l[1],3));
for l_one = l
	ix = 1:size(l_one,2);
	finalData[:,ix,:] .+= l_one;
end

# ## Plotting
# This is what the data looks like, for one channel/trial respectively:
f = Figure()
Axis(f[1, 1], title = "Single channel, all trials", xlabel="time", ylabel = "y")
for i=1:size(finalData,3)
	lines!(finalData[1,:,i])
end
hlines!([0],color=:gray)
Axis(f[2,1], title = "All channels, single trial", xlabel="time", ylabel = "y")
for i=1:size(finalData,1)
	lines!(finalData[i,:,1])
end
hlines!([0],color=:gray)
f

# And some topoplots: 
## first get the positions (channel locations from SEREEGA)
pos3d = CSV.read(download("https://raw.githubusercontent.com/lrkrol/SEREEGA/b5dec01e8b98260904906cc19c38853bfbb2ff19/leadfield/hartmut/chanlocs-hartmut-NYhead_small.xyz"),DataFrame,header=false);
rename!(pos3d,:Column1=>:ix,:Column2=>:x,:Column3=>:y,:Column4=>:z,:Column5=>:label)

## remove some whitespace in the labels
transform!(pos3d,:label=>(x->strip.(x))=>:label);

## remove some weird channels from position
subset!(pos3d,:label=>ByRow(l -> l ∉ ["Nk1","Nk2","Nk3","Nk4"]));

## getting index of these channels from imported hartmut model data, exclude them in the topoplot
remove_channels = findall(l -> l ∈ ["Nk1","Nk2","Nk3","Nk4"], hartmut["electrodes"]["label"]);
remove_indices = getindex.(remove_channels,1) 

out = UnfoldMakie.to_positions(pos3d.x,pos3d.y,pos3d.z);

# now plot
Δbin = 20
t = 1-Δbin;
f_topo = Figure()
for i=1:(size(finalData,2)÷Δbin)
	for j=1:3
		global t = t + Δbin
		if t>size(finalData,2) break;
		end
		ax_topo = Axis(f_topo[i,j], aspect=1, title=string("t=",t))
		eeg_topoplot!(f_topo[i,j],finalData[Not(remove_indices),t,1],pos3d.label;positions=out,enlarge=0.5,colorrange=(-0.0001,0.0001))	
	end
end
f_topo

# Or with UnfoldMakie plot_topoplotseries:
df = UnfoldMakie.eeg_matrix_to_dataframe(finalData[Not(remove_indices),:,1], string.(1:length(out)));
plot_topoplotseries(df, Δbin; positions = out)

# ## ClusterDepth
# Now that the simulation is done, let's try out ClusterDepth and plot our results
pvals = clusterdepth(finalData; τ=1.6,nperm=5);
fig, ax, hm = heatmap(transpose(pvals))
ax.title="pvals";
ax.xlabel="time";
ax.ylabel="channel";
Colorbar(fig[:, end+1], hm);
fig
