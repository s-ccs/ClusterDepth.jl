using ClusterDepth
using Random
using CairoMakie
using UnfoldSim
using StatsBase
using Distributions
using DataFrames
using Unfold
using UnfoldMakie

# # How to use the ClusterDepth multiple comparison correction

# !!! info
#       This tutorial focuses on single-channel data. For multichannel data, see the tutorial "Further EEG Example".
# ## Simulating test-data
# Let's setup an EEG simulation using UnfoldSim.jl. We simulate a one factor 1x2 design with 20 subjects, each with 40 trials
n_subjects = 20
design = MultiSubjectDesign(
    n_subjects=n_subjects,
    n_items=40,
    items_between=Dict(:condition => ["car", "face"]),
)
first(generate_events(design), 3)

# Next we define a ground-truth signal based on a Linear Mixed Model.
#
# We simulate a **P100**, a **N170** and a **P300** - but an effect only on the **N170**
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

data, events = simulate(
    MersenneTwister(1),
    design,
    [p1, n170, p300],
    UniformOnset(; offset=500, width=100),
    RedNoise(noiselevel=1);
    return_epoched=true,
)
times = range(0, stop=size(data, 1) / 250, length=size(data, 1));

# ## Fit a model to each subject
# We now have data, but no multiple comparison problem yet - we have to fit one **analysis** regression model to each time-point of each subject. Thereby, we perform 113 tests and would (without multiple testing correction) expect 5-6 samples to be significant, even if there is no true effect (but there is one ;)).
# !!! note
#       In principle, we do not need Unfold here - we could simply calculate (subjectwise) means of the conditions, and their time-resolved difference. Using Unfold.jl here simply generalizes it to more complex design, e.g. with continuous predictors etc.

models = map(
    (d, ev) -> (
        fit(UnfoldModel, @formula(0 ~ 1 + condition), DataFrame(ev), d, times),
        ev.subject[1],
    ),
    eachslice(data; dims=3),
    groupby(events, :subject),
);


# now we can inspect the data easily, and extract the face-effect
function add_subject!(df, s)
    df[!, :subject] .= s
    return df
end
allEffects =
    map(
        (x) ->
            (effects(Dict(:condition => ["car", "face"]), x[1]), x[2]) |>
            (x) -> add_subject!(x[1], x[2]),
        models,
    ) |> e -> reduce(vcat, e)

plot_erp(allEffects; mapping=(color=:condition, group=:subject))
# Every line is from one subject, the color indicates our two conditions. 
#
# It is easier to see potential differences if we would plot the difference: 
#
# First we extract all coefficients in a nice, tidy DataFrame
allCoefs =
    map(m -> (coeftable(m[1]), m[2]) |>
             (x) -> add_subject!(x[1], x[2]), models) |>
    e -> reduce(vcat, e)
plot_erp(allCoefs; mapping=(group=:subject, col=:coefname))
# This plot now shows the intercept (in our contrast-case the `condition:"car"` ERP), and the difference curve.
#
# Next we unstack the tidy-coef table into a matrix and put it to clusterdepth for clusterpermutation testing
faceCoefs = allCoefs |> x -> subset(x, :coefname => x -> x .== "condition: face")
erpMatrix =
    unstack(faceCoefs, :subject, :time, :estimate) |> x -> Matrix(x[:, 2:end])' |> collect
summary(erpMatrix)

# ## Clusterdepth
pvals = clusterdepth(erpMatrix; τ=quantile(TDist(n_subjects - 1), 0.95), nperm=5000);
# well - that was fast, less than a second for a cluster permutation test. not bad at all!

# ## Plotting
# Some plotting, and we add the identified cluster

# first calculate the ERP
faceERP =
    groupby(faceCoefs, [:time, :coefname]) |>
    x -> combine(x, :estimate => mean => :estimate, :estimate => x -> std(x) / sqrt(length(x)) => :stderror);

# put the pvalues into a nicer format
pvalDF =
    ClusterDepth.cluster(pvals .<= 0.05) |>
    x -> DataFrame(
        :from => x[1] ./ 250,
        :to => (x[1] .+ x[2]) ./ 250,
        :coefname => "condition: face",
    )
plot_erp(faceERP; stderror=true, pvalue=pvalDF)

# Looks good! We identified the cluster :-)