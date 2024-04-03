using ClusterDepth
using Random
using CairoMakie
using UnfoldSim
using StatsBase
using Distributions
using DataFrames
using Unfold
using UnfoldMakie

# ## How to use clusterDepth multiple comparison correction
# !!! info
#       This tutorial focuses on single-channel data. For multichannel data, see the tutorial "Further EEG Example".

# Let's setup an EEG simulation using UnfoldSim.jl. We simulate a simple 1x2 design with 20 subjects, each with 40 trials
n_subjects = 20
design = MultiSubjectDesign(n_subjects=n_subjects, n_items=40, items_between=Dict(:condition => ["car", "face"]))
first(generate_events(design), 5)

# next we define a ground-truth signal + relation to events/design with Wilkinson Formulas
# let's simulate a P100, a N170 and a P300 - but an effect only on the N170
p1 = MixedModelComponent(;
    basis=UnfoldSim.p100(; sfreq=250),
    formula=@formula(0 ~ 1 + (1 | subject)),
    β=[1.0],
    σs=Dict(:subject => [1]),
);
n170 = MixedModelComponent(;
    basis=UnfoldSim.n170(; sfreq=250),
    formula=@formula(0 ~ 1 + condition + (1 + condition | subject)),
    β=[1.0, -0.5], # condition effect - faces are more negative than cars
    σs=Dict(:subject => [1, 0.2]), # random slope yes please!
);

p300 = MixedModelComponent(;
    basis=UnfoldSim.p300(; sfreq=250),
    formula=@formula(0 ~ 1 + condition + (1 + condition | subject)),
    β=[1.0, 0], ## no p300 condition effect
    σs=Dict(:subject => [1, 1.0]), # but a random slope for condition
);

## Start the simulation
data, events = simulate(MersenneTwister(1), design, [p1, n170, p300], UniformOnset(; offset=500, width=100), RedNoise(noiselevel=1); return_epoched=true)
times = range(0, stop=size(data, 1) / 250, length=size(data, 1));

# let's fit an Unfold Model for each subject
# !!! note
#       In principle, we do not need Unfold here - we could simply calculate (subjectwise) means of the conditions, and their time-resolved difference. Using Unfold.jl here simply generalizes it to more complex design, e.g. with continuous predictors etc.

models = map((d, ev) -> (fit(UnfoldModel, @formula(0 ~ 1 + condition), DataFrame(ev), d, times), ev.subject[1]),
    eachslice(data; dims=3),
    groupby(events, :subject))


# now we can inspect the data easily, and extract the face-effect
function add_subject!(df, s)
    df[!, :subject] .= s
    return df
end
allEffects = map((x) -> (effects(Dict(:condition => ["car", "face"]), x[1]), x[2]) |> (x) -> add_subject!(x[1], x[2]), models) |> e -> reduce(vcat, e)

plot_erp(allEffects; mapping=(color=:condition, group=:subject))

# extract the face-coefficient from the linear model
allCoefs = map(m -> (coeftable(m[1]), m[2]) |> (x) -> add_subject!(x[1], x[2]), models) |> e -> reduce(vcat, e)
plot_erp(allCoefs; mapping=(group=:subject, col=:coefname))

# let's unstack the tidy-coef table into a matrix and put it to clusterdepth for clusterpermutation testing
faceCoefs = allCoefs |> x -> subset(x, :coefname => x -> x .== "condition: face")
erpMatrix = unstack(faceCoefs, :subject, :time, :estimate) |> x -> Matrix(x[:, 2:end])' |> collect
summary(erpMatrix)

# ## Clusterdepth
pvals = clusterdepth(erpMatrix; τ=quantile(TDist(n_subjects - 1), 0.95), nperm=5000);
# well - that was fast, less than a second for a cluster permutation test. not bad at all!

# ## Plotting
# Some plotting, and we add the identified cluster

# first calculate the ERP
faceERP = groupby(faceCoefs, [:time, :coefname]) |>
          x -> combine(x, :estimate => mean => :estimate,
    :estimate => std => :stderror);

# put the pvalues into a nicer format
pvalDF = ClusterDepth.cluster(pvals .<= 0.05) |> x -> DataFrame(:from => x[1] ./ 250, :to => (x[1] .+ x[2]) ./ 250, :coefname => "condition: face")
plot_erp(faceERP; stderror=true, pvalue=pvalDF)

# Looks good to me! We identified the cluster :-)

# old unused code to use extra=(;pvalue=pvalDF) in the plotting function, but didnt work.

