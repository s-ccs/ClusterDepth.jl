using ClusterDepth
using Random
using CairoMakie
using UnfoldSim
using StatsBase
using Distributions
using DataFrames

# # Family Wise Error of ClusterDepth
# Here we calculate the Family Wise Error the ClusterDepth Correct
# That is, we want to check that pvalues we get, indeed return us a type-1 of 5% for all time-points
# The point being, that if you do tests on 113 timepoints, the chance that one is significant is not 5% but 
(1-(1-0.05)^113)*100 ##%

# ## Setup Simulation
# Let's setup a simulation using UnfoldSim.jl. We simulate a simple 1x2 design with 20 subjects
n_subjects=20
design = MultiSubjectDesign(n_subjects=n_subjects, n_items = 2, items_between= Dict(:condition=>["small","large"]))
first(generate(design),5)

# 
# Next we define a ground-truth signal + relation to events/design with Wilkinson Formulas.
# we want no condition effect, therefore β for the condition should be 0. We further add some inter-subject variability with the mixed models.
# We will use a simulated P300 signal, which at 250Hz has 113 samples.
signal = MixedModelComponent(;
    basis=UnfoldSim.p300(;sfreq=250),
    formula=@formula(dv ~ 1+condition+(1|subject)),
    β=[1.,0.0],
    σs = Dict(:subject=>[1]),
);

# 
# Let's move the actual simulation into a function, so we can call it many times.
# Note that we use (`RedNoise`)[https://unfoldtoolbox.github.io/UnfoldSim.jl/dev/literate/reference/noisetypes/] which has lot's of Autocorrelation between timepoints. nice!
function run_fun(r)
    data, events = simulate(MersenneTwister(r), design, signal, UniformOnset(; offset=5, width=4), RedNoise(noiselevel=1); return_epoched=true)
    data = data[:,events.condition .== "small"] .- data[:,events.condition .== "large"]
    
    return data,clusterdepth(data'; τ=quantile(TDist(n_subjects - 1), 0.95),nperm=1000)
end;

# ## Understanding the simulation
# let's have a look at the actual data by running it once, plotting condition wise trials, the ERP and histograms of uncorrected and corrected p-values
data,pval = run_fun(5)
conditionSmall = data[:,1:2:end]
conditionLarge = data[:,2:2:end]
pval_uncorrected = 1 .- cdf.(TDist(n_subjects-1),abs.(ClusterDepth.studentt(conditionSmall.-conditionLarge)))
sig = pval_uncorrected.<=0.025;

# For the uncorrected p-values based on the t-distribution, we get a type1 error over "time":
mean(sig)

# this is the type 1 error of 5% we expected.

# !!! note
#       Type-I error is not the FWER (family wise error rate). FWER is the property of a set of tests (in this case tests per time-point), we can calculate it by repeating such tests,
#       and checking for each repetition whether any sample of a repetition is significant (e.g. `any(sig)` followed by a `mean(repetitions_anysig)`).


f = Figure();
series!(Axis(f[1,1],title="condition==small"),conditionSmall',solid_color   =:red)
series!(Axis(f[1,2],title="condition==large"),conditionLarge',solid_color   =:blue)
ax = Axis(f[2,1:2],title="ERP (mean over trials)")

sig = allowmissing(sig)
sig[sig.==0] .= missing
@show sum(skipmissing(sig))
lines!(sig,color=:gray,linewidth=4)
lines!(ax,mean(conditionSmall,dims=2)[:,1],solid_color   =:red)
lines!(ax,mean(conditionLarge,dims=2)[:,1],solid_color   =:blue)

hist!(Axis(f[3,1],title="uncorrected pvalues"),pval_uncorrected,bins=0:0.01:1.1)
hist!(Axis(f[3,2],title="clusterdepth corrected pvalues"),pval,bins=0:0.01:1.1)
f
#----
# ## Run simulations
# This takes some seconds (depending on your infrastructure)
reps = 500
res = fill(NaN, reps, 2)
Threads.@threads for r = 1:reps
    data,pvals = run_fun(r)
    res[r, 1] = mean(pvals .<= 0.05)
    res[r, 2] = mean(abs.(ClusterDepth.studentt(data)) .>= quantile(TDist(n_subjects - 1), 0.975))
end;
# Finally, let's calculate the percentage of simulations where we find a significant effect somewhere
mean(res.>0,dims=1) |> x-> (:clusterdepth=>x[1],:uncorrected=>x[2])

# Nice, correction seems to work in principle :) Clusterdepth is not exactly 5%, but with more repetitions we should get there (e.g. with 5000 repetitions, we get 0.051%).

# !!! info
#       if you look closely, the `:uncorrected` value (around 60%) is not as bad as the 99% promised in the introduction. This is due to the correlation between the tests introduced by the noise. Indeed, a good exercise is to repeat everything, but put `RedNoise` to `WhiteNoise`