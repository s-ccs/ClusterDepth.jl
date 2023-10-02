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
#       Currently ClusterDepth.jl supports only a single channel. Multi-Channel support is described in the paper, but I dont know how well that algorithm works.

# Let's setup an EEG simulation using UnfoldSim.jl. We simulate a simple 1x2 design with 20 subjects, each with 40 trials
n_subjects=20
design = MultiSubjectDesign(n_subjects=n_subjects, n_items = 40, items_between= Dict(:condition=>["car","face"]))
first(generate(design),5)

# next we define a ground-truth signal + relation to events/design with Wilkinson Formulas
# let's simulate a P100, a N170 and a P300 - but an effect only on the N170
p100 = MixedModelComponent(;
    basis=UnfoldSim.p100(;sfreq=250),
    formula=@formula(0 ~ 1+(1|subject)),
    β=[1.],
    σs = Dict(:subject=>[1]),
);
n170 = MixedModelComponent(;
    basis=UnfoldSim.n170(;sfreq=250),
    formula=@formula(0 ~ 1+condition+(1+condition|subject)),
    β=[1.,-0.5], # condition effect - faces are more negative than cars
    σs = Dict(:subject=>[1,0.2]), # random slope yes please!
);

p300 = MixedModelComponent(;
    basis=UnfoldSim.p300(;sfreq=250),
    formula=@formula(0 ~ 1+condition+(1+condition|subject)),
    β=[1.,0], ## no p300 condition effect
    σs = Dict(:subject=>[1,1.]), # but a random slope for condition
);

## Start the simulation
data, events = simulate(MersenneTwister(1), design, [p100,n170,p300], UniformOnset(; offset=5, width=4), RedNoise(noiselevel=1); return_epoched=true)
times = range(0,stop=size(data,1)/250,length=size(data,1));

# because the writer of this documentation is also  the author of Unfold, we are using Unfold on each subject to get one difference curve per subject
allData = hcat(events,DataFrame(:eeg=>[[col...] for col in eachcol(data)]));

# let's fit an Unfold Model for each subject
# !!! note
#       In principle, we do not need Unfold here - we could simply calculate (subjectwise) means of the conditions, and their time-resolved difference. Using Unfold.jl here simply generalizes it to more complex design, e.g. with continuous predictors etc.
models = groupby(allData,:subject) |> 
    x->combine(x,AsTable(:)=>(oneSubject->fit(UnfoldModel,
                                            @formula(0~1+condition),
                                            DataFrame(oneSubject),##events
                                            reshape(hcat(oneSubject.eeg...),1,length(oneSubject.eeg[1]),length(oneSubject.eeg)),
                                            times)) =>"unfoldModel");


# now we can inspect the data easily, and extract the face-effect

allEffects = combine(groupby(models,:subject),
    :unfoldModel=>(x->effects(Dict(:condition=>["car","face"]),x[1]))=>AsTable)
plot_erp(allEffects;mapping=(:color=>:condition,:group=>:subject))

# extract the face-coefficient from the linear model
allCoefs = combine(groupby(models,:subject),:unfoldModel=>(x->coeftable(x[1]))=>AsTable)
plot_erp(allCoefs;mapping=(:group=>:subject,col=:coefname))

# let's unstack the tidy-coef table into a matrix and put it to clusterdepth for clusterpermutation testing
faceCoefs = allCoefs|>x->subset(x,:coefname=>x->x.=="condition: face")
erpMatrix= unstack(faceCoefs,:subject,:time,:estimate)|>x->Matrix(x[:,2:end])'|>collect
summary(erpMatrix)

# ## Clusterdepth
pvals = clusterdepth(erpMatrix; τ=quantile(TDist(n_subjects - 1), 0.95),nperm=5000);
# well - that was fast, less than a second for a cluster permutation test. not bad at all!

# ## Plotting
# Some plotting, and we add the identified cluster

# first calculate the ERP
faceERP = groupby(faceCoefs,[:time,:coefname])|>
    x->combine(x,:estimate=>mean=>:estimate,
                 :estimate=>std=>:stderror);

# put the pvalues into a nicer format (in the future UnfoldMakie should do)
plot_erp(faceERP; extra=(;stderror=true))
significant = (pvals.<=0.05) |>allowmissing
significant[significant.==0] .= missing
lines!(times,significant.*-1,color=:darkblue,linewidth=4)
hlines!([0],color=:gray)
current_figure()

# Looks good to me! We identified the cluster :-)

# old unused code to use extra=(;pvalue=pvalDF) in the plotting function, but didnt work.
#pvalDF = ClusterDepth.cluster(pvals.<=0.05)|>x->DataFrame(:from=>x[1],:to=>x[1].+x[2],:coefname=>"condition: face")
