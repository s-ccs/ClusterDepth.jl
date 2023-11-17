var documenterSearchIndex = {"docs":
[{"location":"tutorials/demo/","page":"-","title":"-","text":"EditURL = \"https://github.com/s-ccs/ClusterDepth.jl/blob/main/docs/src/tutorials/demo.jl\"","category":"page"},{"location":"tutorials/demo/","page":"-","title":"-","text":"using ClusterDepth\nusing Random\nusing CairoMakie\n\n\nn_t =40 # timepoints\nn_sub = 50\nn_perm = 5000\n\nsnr = 0.5 # signal to nois\n\n# add a signal to the middle\nsignal = vcat(zeros(n_t÷4), sin.(range(0,π,length=n_t÷2)), zeros(n_t÷4))\n\n# same signal for all subs\nsignal = repeat(signal,1,n_sub)\n\n\n# add noise\ndata = randn(MersenneTwister(123),n_t,n_sub).+ snr .* signal\n\n# by default assumes τ=2.3 (~alpha=0.05), and one-sample ttest\n@time pvals = clusterdepth(data);\n\nf = Figure()\nax = f[1,1] = Axis(f)\n\n\nlines!(abs.(ClusterDepth.studentt(data)))\nh1 = scatter(f[1,2],pvals;axis=(;yscale=log10),label=\"troendle\")\n\npvals2 = clusterdepth(data;pval_type=:naive)\nh2 = scatter!(1.2:40.2,pvals2,color=\"red\",label=\"naive\")\nhlines!([0.05,0.01])\naxislegend()\nf","category":"page"},{"location":"tutorials/demo/","page":"-","title":"-","text":"","category":"page"},{"location":"tutorials/demo/","page":"-","title":"-","text":"This page was generated using Literate.jl.","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"EditURL = \"https://github.com/s-ccs/ClusterDepth.jl/blob/main/docs/src/reference/type1_troendle.jl\"","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"using ClusterDepth\nusing Random\nusing CairoMakie\nusing UnfoldSim\nusing StatsBase\nusing ProgressMeter\nusing Distributions","category":"page"},{"location":"reference/type1_troendle/#Family-Wise-Error-of-Troendle","page":"Troendle FWER","title":"Family Wise Error of Troendle","text":"","category":"section"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Here we calculate the Family Wise Error of doing ntests at the same time. That is, we want to check that Troendle indeed returns us a type-1 of 5% for a set of tests.","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"The point being, that if you do 30 tests, the chance that one is significant is not 5% but actually","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"(1-(1-0.05)^30)*100 ##%","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Let's setup some simulation parameters","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"reps = 1000\nperms = 1000\nntests = 30;\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"we will use the student-t in it's 2-sided variant (abs of it)","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"fun = x->abs.(ClusterDepth.studentt(x));\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"this function simulates data without any effect (H0), then the permutations, and finally calls troendle","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"function run_fun(r,perms,fun,ntests)\n    rng = MersenneTwister(r)\n    data = randn(rng,ntests,50)\n    perm = Matrix{Float64}(undef,size(data,1),perms)\n    for p = 1:perms\n        perm[:,p] = ClusterDepth.sign_permute(rng,data,fun)\n    end\n    stat = fun(data)\n    return data,stat, troendle(perm,stat)\nend;\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"let's test it once","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"data,stats_t,pvals = run_fun(1,perms,fun,ntests);\nprintln(\"data:\", size(data),\" t-stats:\",size(stats_t),\" pvals:\",size(pvals))","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"run the above function reps=1000` times - we also save the uncorrected t-based pvalue","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"pvals_all = fill(NaN,reps,2,ntests)\n@Threads.threads  for r = 1:reps\n    data,stat,pvals = run_fun(r,perms,fun,ntests)\n    pvals_all[r,1,:] = pvals\n    pvals_all[r,2,:] = (1 .- cdf.(TDist(size(data,2)),abs.(stat))).*2 # * 2 becaue of twosided. Troendle takes this into account already\nend;\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Let's check in how many of our simlations we have a significant p-value =<0.05","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"res = any(pvals_all[:,:,:] .<= 0.05,dims=3)[:,:,1]\nmean(res.>0,dims=1) |> x-> (:troendle=>x[1],:uncorrected=>x[2])","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Nice. Troendle fits perfectly and the uncorrected is pretty close to what we calculated above!","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Finally we end this with a short figure to get a better idea of how this data looks like and a histogram of the p-values","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"f = Figure()\nax = f[1,1] = Axis(f)\n\n\n\nlines!(ax,abs.(ClusterDepth.studentt(data)))\nheatmap!(Axis(f[2,1]),data)\nseries!(Axis(f[2,2]),data[:,1:7]')\nh1 = scatter!(Axis(f[1,2];yscale=log10),pvals,label=\"troendle\")\n\nhlines!([0.05,0.01])\n\nhist!(Axis(f[3,1]),pvals_all[:, 1,:][:],bins=0:0.01:1.1)\nhist!(Axis(f[3,2]),pvals_all[:, 2,:][:],bins=0:0.01:1.1)\nf","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"This page was generated using Literate.jl.","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"EditURL = \"https://github.com/s-ccs/ClusterDepth.jl/blob/main/docs/src/reference/type1.jl\"","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"using ClusterDepth\nusing Random\nusing CairoMakie\nusing UnfoldSim\nusing StatsBase\nusing Distributions\nusing DataFrames","category":"page"},{"location":"reference/type1/#Family-Wise-Error-of-ClusterDepth","page":"Clusterdepth FWER","title":"Family Wise Error of ClusterDepth","text":"","category":"section"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Here we calculate the Family Wise Error the ClusterDepth Correct That is, we want to check that pvalues we get, indeed return us a type-1 of 5% for all time-points The point being, that if you do tests on 113 timepoints, the chance that one is significant is not 5% but","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"(1-(1-0.05)^113)*100 ##%","category":"page"},{"location":"reference/type1/#Setup-Simulation","page":"Clusterdepth FWER","title":"Setup Simulation","text":"","category":"section"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Let's setup a simulation using UnfoldSim.jl. We simulate a simple 1x2 design with 20 subjects","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"n_subjects=20\ndesign = MultiSubjectDesign(n_subjects=n_subjects, n_items = 2, items_between= Dict(:condition=>[\"small\",\"large\"]))\nfirst(generate(design),5)","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Next we define a ground-truth signal + relation to events/design with Wilkinson Formulas. we want no condition effect, therefore β for the condition should be 0. We further add some inter-subject variability with the mixed models. We will use a simulated P300 signal, which at 250Hz has 113 samples.","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"signal = MixedModelComponent(;\n    basis=UnfoldSim.p300(;sfreq=250),\n    formula=@formula(0 ~ 1+condition+(1|subject)),\n    β=[1.,0.0],\n    σs = Dict(:subject=>[1]),\n);\nnothing #hide","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Let's move the actual simulation into a function, so we can call it many times. Note that we use (RedNoise)[https://unfoldtoolbox.github.io/UnfoldSim.jl/dev/literate/reference/noisetypes/] which has lot's of Autocorrelation between timepoints. nice!","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"function run_fun(r)\n    data, events = simulate(MersenneTwister(r), design, signal, UniformOnset(; offset=5, width=4), RedNoise(noiselevel=1); return_epoched=true)\n    data = data[:,events.condition .== \"small\"] .- data[:,events.condition .== \"large\"]\n    return data,clusterdepth(data; τ=quantile(TDist(n_subjects - 1), 0.95),nperm=1000)\nend;\nnothing #hide","category":"page"},{"location":"reference/type1/#Understanding-the-simulation","page":"Clusterdepth FWER","title":"Understanding the simulation","text":"","category":"section"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"let's have a look at the actual data by running it once, plotting condition wise trials, the ERP and histograms of uncorrected and corrected p-values","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"data,pval = run_fun(5)\nconditionSmall = data[:,1:2:end]\nconditionLarge = data[:,2:2:end]\npval_uncorrected = 1 .- cdf.(TDist(n_subjects-1),abs.(ClusterDepth.studentt(conditionSmall.-conditionLarge)))\nsig = pval_uncorrected.<=0.025;\nnothing #hide","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"For the uncorrected p-values based on the t-distribution, we get a type1 error over \"time\":","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"mean(sig)","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"this is the type 1 error of 5% we expected.","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"note: Note\nType-I error is not the FWER (family wise error rate). FWER is the property of a set of tests (in this case tests per time-point), we can calculate it by repeating such tests,   and checking for each repetition whether any sample of a repetition is significant (e.g. any(sig) followed by a mean(repetitions_anysig)).","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"f = Figure();\nseries!(Axis(f[1,1],title=\"condition==small\"),conditionSmall',solid_color   =:red)\nseries!(Axis(f[1,2],title=\"condition==large\"),conditionLarge',solid_color   =:blue)\nax = Axis(f[2,1:2],title=\"ERP (mean over trials)\")\n\nsig = allowmissing(sig)\nsig[sig.==0] .= missing\n@show sum(skipmissing(sig))\nlines!(sig,color=:gray,linewidth=4)\nlines!(ax,mean(conditionSmall,dims=2)[:,1],solid_color   =:red)\nlines!(ax,mean(conditionLarge,dims=2)[:,1],solid_color   =:blue)\n\nhist!(Axis(f[3,1],title=\"uncorrected pvalues\"),pval_uncorrected,bins=0:0.01:1.1)\nhist!(Axis(f[3,2],title=\"clusterdepth corrected pvalues\"),pval,bins=0:0.01:1.1)\nf","category":"page"},{"location":"reference/type1/#Run-simulations","page":"Clusterdepth FWER","title":"Run simulations","text":"","category":"section"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"This takes some seconds (depending on your infrastructure)","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"reps = 500\nres = fill(NaN, reps, 2)\nThreads.@threads for r = 1:reps\n    data,pvals = run_fun(r)\n    res[r, 1] = mean(pvals .<= 0.05)\n    res[r, 2] = mean(abs.(ClusterDepth.studentt(data)) .>= quantile(TDist(n_subjects - 1), 0.975))\nend;\nnothing #hide","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Finally, let's calculate the percentage of simulations where we find a significant effect somewhere","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"mean(res.>0,dims=1) |> x-> (:clusterdepth=>x[1],:uncorrected=>x[2])","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Nice, correction seems to work in principle :) Clusterdepth is not exactly 5%, but with more repetitions we should get there (e.g. with 5000 repetitions, we get 0.051%).","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"info: Info\nif you look closely, the :uncorrected value (around 60%) is not as bad as the 99% promised in the introduction. This is due to the correlation between the tests introduced by the noise. Indeed, a good exercise is to repeat everything, but put RedNoise to WhiteNoise","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"This page was generated using Literate.jl.","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"EditURL = \"https://github.com/s-ccs/ClusterDepth.jl/blob/main/docs/src/tutorials/eeg.jl\"","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"using ClusterDepth\nusing Random\nusing CairoMakie\nusing UnfoldSim\nusing StatsBase\nusing Distributions\nusing DataFrames\nusing Unfold\nusing UnfoldMakie","category":"page"},{"location":"tutorials/eeg/#How-to-use-clusterDepth-multiple-comparison-correction","page":"An EEG Example","title":"How to use clusterDepth multiple comparison correction","text":"","category":"section"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"info: Info\nCurrently ClusterDepth.jl supports only a single channel. Multi-Channel support is described in the paper, but I dont know how well that algorithm works.","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"Let's setup an EEG simulation using UnfoldSim.jl. We simulate a simple 1x2 design with 20 subjects, each with 40 trials","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"n_subjects=20\ndesign = MultiSubjectDesign(n_subjects=n_subjects, n_items = 40, items_between= Dict(:condition=>[\"car\",\"face\"]))\nfirst(generate(design),5)","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"next we define a ground-truth signal + relation to events/design with Wilkinson Formulas let's simulate a P100, a N170 and a P300 - but an effect only on the N170","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"p100 = MixedModelComponent(;\n    basis=UnfoldSim.p100(;sfreq=250),\n    formula=@formula(0 ~ 1+(1|subject)),\n    β=[1.],\n    σs = Dict(:subject=>[1]),\n);\nn170 = MixedModelComponent(;\n    basis=UnfoldSim.n170(;sfreq=250),\n    formula=@formula(0 ~ 1+condition+(1+condition|subject)),\n    β=[1.,-0.5], # condition effect - faces are more negative than cars\n    σs = Dict(:subject=>[1,0.2]), # random slope yes please!\n);\n\np300 = MixedModelComponent(;\n    basis=UnfoldSim.p300(;sfreq=250),\n    formula=@formula(0 ~ 1+condition+(1+condition|subject)),\n    β=[1.,0], ## no p300 condition effect\n    σs = Dict(:subject=>[1,1.]), # but a random slope for condition\n);\n\n# Start the simulation\ndata, events = simulate(MersenneTwister(1), design, [p100,n170,p300], UniformOnset(; offset=5, width=4), RedNoise(noiselevel=1); return_epoched=true)\ntimes = range(0,stop=size(data,1)/250,length=size(data,1));\nnothing #hide","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"because the writer of this documentation is also  the author of Unfold, we are using Unfold on each subject to get one difference curve per subject","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"allData = hcat(events,DataFrame(:eeg=>[[col...] for col in eachcol(data)]));\nnothing #hide","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"let's fit an Unfold Model for each subject","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"note: Note\nIn principle, we do not need Unfold here - we could simply calculate (subjectwise) means of the conditions, and their time-resolved difference. Using Unfold.jl here simply generalizes it to more complex design, e.g. with continuous predictors etc.","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"models = groupby(allData,:subject) |>\n    x->combine(x,AsTable(:)=>(oneSubject->fit(UnfoldModel,\n                                            @formula(0~1+condition),\n                                            DataFrame(oneSubject),##events\n                                            reshape(hcat(oneSubject.eeg...),1,length(oneSubject.eeg[1]),length(oneSubject.eeg)),\n                                            times)) =>\"unfoldModel\");\nnothing #hide","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"now we can inspect the data easily, and extract the face-effect","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"allEffects = combine(groupby(models,:subject),\n    :unfoldModel=>(x->effects(Dict(:condition=>[\"car\",\"face\"]),x[1]))=>AsTable)\nplot_erp(allEffects;mapping=(:color=>:condition,:group=>:subject))","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"extract the face-coefficient from the linear model","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"allCoefs = combine(groupby(models,:subject),:unfoldModel=>(x->coeftable(x[1]))=>AsTable)\nplot_erp(allCoefs;mapping=(:group=>:subject,col=:coefname))","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"let's unstack the tidy-coef table into a matrix and put it to clusterdepth for clusterpermutation testing","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"faceCoefs = allCoefs|>x->subset(x,:coefname=>x->x.==\"condition: face\")\nerpMatrix= unstack(faceCoefs,:subject,:time,:estimate)|>x->Matrix(x[:,2:end])'|>collect\nsummary(erpMatrix)","category":"page"},{"location":"tutorials/eeg/#Clusterdepth","page":"An EEG Example","title":"Clusterdepth","text":"","category":"section"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"pvals = clusterdepth(erpMatrix; τ=quantile(TDist(n_subjects - 1), 0.95),nperm=5000);\nnothing #hide","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"well - that was fast, less than a second for a cluster permutation test. not bad at all!","category":"page"},{"location":"tutorials/eeg/#Plotting","page":"An EEG Example","title":"Plotting","text":"","category":"section"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"Some plotting, and we add the identified cluster","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"first calculate the ERP","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"faceERP = groupby(faceCoefs,[:time,:coefname])|>\n    x->combine(x,:estimate=>mean=>:estimate,\n                 :estimate=>std=>:stderror);\nnothing #hide","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"put the pvalues into a nicer format (in the future UnfoldMakie should do)","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"plot_erp(faceERP; extra=(;stderror=true))\nsignificant = (pvals.<=0.05) |>allowmissing\nsignificant[significant.==0] .= missing\nlines!(times,significant.*-1,color=:darkblue,linewidth=4)\nhlines!([0],color=:gray)\ncurrent_figure()","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"Looks good to me! We identified the cluster :-)","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"old unused code to use extra=(;pvalue=pvalDF) in the plotting function, but didnt work.","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"#pvalDF = ClusterDepth.cluster(pvals.<=0.05)|>x->DataFrame(:from=>x[1],:to=>x[1].+x[2],:coefname=>\"condition: face\")","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"","category":"page"},{"location":"tutorials/eeg/","page":"An EEG Example","title":"An EEG Example","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = ClusterDepth","category":"page"},{"location":"#ClusterDepth","page":"Home","title":"ClusterDepth","text":"","category":"section"},{"location":"#Comparison-to-permuco-R-implementation","page":"Home","title":"Comparison to permuco R implementation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The implementation to Permuco is similar, but ClusterDepth.jl is more barebone - that is, we dont offer many permutation schemes, focus on the ClusterDepth Algorithm, and don't provide the nice wrappers like clusterLM.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Timing wise, a simple test on 50 subjects, 100 repetitions, 5000 permutations shows the following results:","category":"page"},{"location":"","page":"Home","title":"Home","text":"timepoints ClusterDepth.jl permuco julia-speedup\n40 0.1s 2.9s 29x\n400 0.6s 22s 36x\n4000 7s 240s 34x","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ClusterDepth]","category":"page"},{"location":"#ClusterDepth.calc_clusterdepth-Tuple{AbstractMatrix{<:Real}, Any}","page":"Home","title":"ClusterDepth.calc_clusterdepth","text":"calc_clusterdepth(data,τ)\n\nreturns tuple with three entries: 1:maxLength, maximal clustervalue per clusterdepth head, same for tail\n\nWe assume data and τ have already been transformed for one/two sided testing, so that we can do d0.>τ for finding clusters\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.cluster-Tuple{Any}","page":"Home","title":"ClusterDepth.cluster","text":"finds neighbouring clusters in the vector and returns start + length vectors\n\nif the first and last cluster start on the first/last sample, we dont know their real depth\n\nInput is assumed to be a thresholded Array with only 0/1\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.clusterdepth-Tuple{AbstractArray, Vararg{Any}}","page":"Home","title":"ClusterDepth.clusterdepth","text":"clusterdepth(rng,data::AbstractArray;τ=2.3, statFun=x->abs.(studentt(x)),permFun=signpermute!,nperm=5000,pvaltype=:troendle)\n\ncalculate clusterdepth of given datamatrix. \n\ndata: statFun will be applied on second dimension of data (typically this will be subjects)\n\nOptional\n\nτ: Cluster-forming threshold \nstatFun: default  the one-sample studenttest, can be any custom function on a Matrix returning a Vector\npermFun: default to sign-flip (for one-sample case)\nnperm: number of permutations, default 5000\npval_type: how to calculate pvalues within each cluster, default :troendle, see ?pvals\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.ix_sortUnique-Tuple{Any}","page":"Home","title":"ClusterDepth.ix_sortUnique","text":"in some sense: argsort(argunique(x)), returns the indices to get a sorted unique of x\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.multicol_minimum-Tuple{AbstractMatrix, AbstractVector}","page":"Home","title":"ClusterDepth.multicol_minimum","text":"calculates the minimum in `X` along `dims=2` in the columns specified by àrrayOfIndicearrays` which could be e.g. `[[1,2],[5,6],[3,4,7]]`\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.pvals-Tuple{Any}","page":"Home","title":"ClusterDepth.pvals","text":"pvals(data;kwargs...) = pvals(data[2:end],data[1];kwargs...)\n\npvals(data::AbstractVector,stat::Real;type=:twosided)\n\ncalculates pvalues based on permutation results\n\nif called with stat, first entry is assumed to be the observation \n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.pvals-Tuple{Matrix, Vararg{Any}}","page":"Home","title":"ClusterDepth.pvals","text":"Calculate pvals from cluster-depth permutation matrices\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.sign_permute!-Tuple{Any, AbstractArray, Any}","page":"Home","title":"ClusterDepth.sign_permute!","text":"Permutation via random sign-flip Flips signs along the last dimension\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.troendle-Tuple{AbstractMatrix, AbstractVector}","page":"Home","title":"ClusterDepth.troendle","text":"function troendle(perm::AbstractMatrix,stat::AbstractVector;type=:twosided)\n\nMultiple Comparison Correction as in Troendle 1995\n\nperm with size  ntests x nperms\n\nstat with size ntests\n\ntype can be :twosided (default), :lesser, :greater\n\nHeavily inspired by the R implementation in permuco from Jaromil Frossard\n\nNote: While permuco is released under BSD, the author Jaromil Frossard gave us an MIT license for the troendle and the clusterdepth R-functions.\n\n\n\n\n\n","category":"method"}]
}
