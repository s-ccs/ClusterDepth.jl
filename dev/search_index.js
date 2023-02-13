var documenterSearchIndex = {"docs":
[{"location":"tutorials/demo/","page":"Demo","title":"Demo","text":"EditURL = \"https://github.com/s-ccs/ClusterDepth.jl/blob/main/docs/src/tutorials/demo.jl\"","category":"page"},{"location":"tutorials/demo/","page":"Demo","title":"Demo","text":"using ClusterDepth\nusing Random\nusing CairoMakie\n\n\nn_t =40 # timepoints\nn_sub = 50\nn_perm = 5000\n\nsnr = 0.5 # signal to nois\n\n# add a signal to the middle\nsignal = vcat(zeros(n_t÷4), sin.(range(0,π,length=n_t÷2)), zeros(n_t÷4))\n\n# same signal for all subs\nsignal = repeat(signal,1,n_sub)\n\n\n# add noise\ndata = randn(MersenneTwister(123),n_t,n_sub).+ snr .* signal\n\n# by default assumes τ=2.3 (~alpha=0.05), and one-sample ttest\n@time pvals = clusterdepth(data);\n\nf = Figure()\nax = f[1,1] = Axis(f)\n\n\nlines!(abs.(ClusterDepth.studentt(data)))\nh1 = scatter(f[1,2],pvals;axis=(;yscale=log10),label=\"troendle\")\n\npvals2 = clusterdepth(data;pval_type=:naive)\nh2 = scatter!(1.2:40.2,pvals2,color=\"red\",label=\"naive\")\nhlines!([0.05,0.01])\naxislegend()\nf","category":"page"},{"location":"tutorials/demo/","page":"Demo","title":"Demo","text":"","category":"page"},{"location":"tutorials/demo/","page":"Demo","title":"Demo","text":"This page was generated using Literate.jl.","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"EditURL = \"https://github.com/s-ccs/ClusterDepth.jl/blob/main/docs/src/reference/type1_troendle.jl\"","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"using ClusterDepth\nusing Random\nusing CairoMakie\nusing UnfoldSim\nusing StatsBase\nusing ProgressMeter\nusing Distributions","category":"page"},{"location":"reference/type1_troendle/#Family-Wise-Error-of-Troendle","page":"Troendle FWER","title":"Family Wise Error of Troendle","text":"","category":"section"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Here we calculate the Family Wise Error of doing ntests at the same time. That is, we want to check that Troendle indeed returns us a type-1 of 5% for a set of tests.","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"The point being, that if you do 30 tests, the chance that one is significant is not 5% but actually","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"(1-(1-0.05)^30)*100 ##%","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Let's setup some simulation parameters","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"reps = 1000\nperms = 1000\nntests = 30;\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"we will use the student-t in it's 2-sided variant (abs of it)","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"fun = x->abs.(ClusterDepth.studentt(x));\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"this function simulates data without any effect (H0), then the permutations, and finally calls troendle","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"function run_fun(r,perms,fun,ntests)\n    rng = MersenneTwister(r)\n    data = randn(rng,ntests,50)\n    perm = Matrix{Float64}(undef,size(data,1),perms)\n    for p = 1:perms\n        perm[:,p] = ClusterDepth.sign_permute(rng,data,fun)\n    end\n    stat = fun(data)\n    return data,stat, troendle(perm,stat)\nend;\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"let's test it once","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"data,stats_t,pvals = run_fun(1,perms,fun,ntests);\nprintln(\"data:\", size(data),\" t-stats:\",size(stats_t),\" pvals:\",size(pvals))","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"run the above function reps=1000` times - we also save the uncorrected t-based pvalue","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"pvals_all = fill(NaN,reps,2,ntests)\n@Threads.threads  for r = 1:reps\n    data,stat,pvals = run_fun(r,perms,fun,ntests)\n    pvals_all[r,1,:] = pvals\n    pvals_all[r,2,:] = (1 .- cdf.(TDist(size(data,2)),abs.(stat))).*2 # * 2 becaue of twosided. Troendle takes this into account already\nend;\nnothing #hide","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Let's check in how many of our simlations we have a significant p-value =<0.05","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"res = any(pvals_all[:,:,:] .<= 0.05,dims=3)[:,:,1]\nmean(res.>0,dims=1) |> x-> (:troendle=>x[1],:uncorrected=>x[2])","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Nice. Troendle fits perfectly and the uncorrected is pretty close to what we calculated above!","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"Finally we end this with a short figure to get a better idea of how this data looks like and a histogram of the p-values","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"f = Figure()\nax = f[1,1] = Axis(f)\n\n\n\nlines!(ax,abs.(ClusterDepth.studentt(data)))\nheatmap!(Axis(f[2,1]),data)\nseries!(Axis(f[2,2]),data[:,1:7]')\nh1 = scatter!(Axis(f[1,2];yscale=log10),pvals,label=\"troendle\")\n\nhlines!([0.05,0.01])\n\nhist!(Axis(f[3,1]),pvals_all[:, 1,:][:],bins=0:0.01:1.1)\nhist!(Axis(f[3,2]),pvals_all[:, 2,:][:],bins=0:0.01:1.1)\nf","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"","category":"page"},{"location":"reference/type1_troendle/","page":"Troendle FWER","title":"Troendle FWER","text":"This page was generated using Literate.jl.","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"EditURL = \"https://github.com/s-ccs/ClusterDepth.jl/blob/main/docs/src/reference/type1.jl\"","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"using ClusterDepth\nusing Random\nusing CairoMakie\nusing UnfoldSim\nusing StatsBase\nusing Distributions\nusing DataFrames","category":"page"},{"location":"reference/type1/#Family-Wise-Error-of-ClusterDepth","page":"Clusterdepth FWER","title":"Family Wise Error of ClusterDepth","text":"","category":"section"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Here we calculate the Family Wise Error the ClusterDepth Correct That is, we want to check that pvalues we get, indeed return us a type-1 of 5% for all time-points The point being, that if you do tests on 113 timepoints, the chance that one is significant is not 5% but","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"(1-(1-0.05)^113)*100 ##%","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Let's setup a simulation using UnfoldSim.jl. We simulate a simple 1x2 design with 20 subjects","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"n_subjects=20\ndesign = MultiSubjectDesign(n_subjects=n_subjects, n_items = 2, items_between= Dict(:condition=>[\"small\",\"large\"]))\nfirst(generate(design),5)","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"next define a ground-truth signal + relation to events/design with Wilkinson Formulas we want no condition effect, so that β should be 0. We further add some inter-subject variability with the mixed models We will use a simulated P300 signal, which has 113 samples","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"signal = MixedModelComponent(;\n    basis=UnfoldSim.p300(;sfreq=250),\n    formula=@formula(dv ~ 1+condition+(1|subject)),\n    β=[1.,0.0],\n    σs = Dict(:subject=>[1]),\n);\nnothing #hide","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"setup simulation function Note that we use RedNoise which has lot's of Autocorrelation between timepoints. nice!","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"function run_fun(r)\n    data, events = simulate(MersenneTwister(r), design, signal, UniformOnset(; offset=5, width=4), RedNoise(noiselevel=1); return_epoched=true)\n    data = data[:,events.condition .== \"small\"] .- data[:,events.condition .== \"large\"]\n    return data,clusterdepth(data; τ=quantile(TDist(n_subjects - 1), 0.95),nperm=1000)\nend;\nnothing #hide","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"let's have a look at the actual data by running it once, plotting condition wise trials, the ERP and histograms of uncorrected and corrected p-values","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"data,pval = run_fun(5)\nconditionSmall = data[:,1:2:end]\nconditionLarge = data[:,2:2:end]\npval_uncorrected = 1 .- cdf.(TDist(n_subjects-1),abs.(ClusterDepth.studentt(conditionSmall.-conditionLarge)))\n\nf = Figure();\nseries!(Axis(f[1,1],title=\"condition==small\"),conditionSmall',solid_color   =:red)\nseries!(Axis(f[1,2],title=\"condition==large\"),conditionLarge',solid_color   =:blue)\nax = Axis(f[2,1:2],title=\"ERP (mean over trials)\")\nsig = pval_uncorrected.<=0.025\nsig = allowmissing(sig)\nsig[sig.==0] .= missing\n@show sum(skipmissing(sig))\nlines!(sig,color=:gray,linewidth=4)\nlines!(ax,mean(conditionSmall,dims=2)[:,1],solid_color   =:red)\nlines!(ax,mean(conditionLarge,dims=2)[:,1],solid_color   =:blue)\n\nhist!(Axis(f[3,1],title=\"uncorrected pvalues\"),pval_uncorrected,bins=0:0.01:1.1)\nhist!(Axis(f[3,2],title=\"clusterdepth corrected pvalues\"),pval,bins=0:0.01:1.1)\nf","category":"page"},{"location":"reference/type1/#run-simulation","page":"Clusterdepth FWER","title":"run simulation","text":"","category":"section"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"This takes some seconds (depending on your infrastructure)","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"reps = 500\nres = fill(NaN, reps, 2)\nThreads.@threads for r = 1:reps\n    data,pvals = run_fun(r)\n    res[r, 1] = mean(pvals .<= 0.05)\n    res[r, 2] = mean(abs.(ClusterDepth.studentt(data)) .>= quantile(TDist(n_subjects - 1), 0.975))\nend;\nnothing #hide","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Finally, let's calculate the percentage of simulations where we find a significant effect somewhere","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"mean(res.>0,dims=1) |> x-> (:clusterdepth=>x[1],:uncorrected=>x[2])","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"Nice, correction seems to work in principle :) Clusterdepth is not exactly 5%, but with more repetitions we should get there (e.g. with 5000 repetitions, we get 0.051%).","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"info: Info\nif you look closely, the :uncorrected value (around 60%) is not as bad as the 99% promised in the introduction. This is due to the correlation between the tests introduced by the noise. Indeed, a good exercise is to repeat everything, but put RedNoise to WhiteNoise","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"","category":"page"},{"location":"reference/type1/","page":"Clusterdepth FWER","title":"Clusterdepth FWER","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = ClusterDepth","category":"page"},{"location":"#ClusterDepth","page":"Home","title":"ClusterDepth","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ClusterDepth.","category":"page"},{"location":"#Comparison-to-permuco-R-implementation","page":"Home","title":"Comparison to permuco R implementation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The implementation to Permuco is similar, but ClusterDepth.jl is more barebone - that is, we dont offer many permutation schemes, focus on the ClusterDepth Algorithm, and don't provide the nice wrappers like clusterLM.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Timing wise, a simple test on 50 subjects, 100 repetitions, 5000 permutations shows the following results:","category":"page"},{"location":"","page":"Home","title":"Home","text":"timepoints ClusterDepth.jl permuco julia-speedup\n40 0.1s 2.9s 29x\n400 0.6s 22s 36x\n4000 7s 240s 34x","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ClusterDepth]","category":"page"},{"location":"#ClusterDepth.calc_clusterdepth-Tuple{Any, Any}","page":"Home","title":"ClusterDepth.calc_clusterdepth","text":"calc_clusterdepth(data,τ)\n\nreturns tuple with three entries: 1:maxLength, maximal clustervalue per clusterdepth head, same for tail\n\nWe assume data and τ have already been transformed for one/two sided testing, so that we can do d0.>τ for finding clusters\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.cluster-Tuple{Any}","page":"Home","title":"ClusterDepth.cluster","text":"finds neighbouring clusters in the vector and returns start + length vectors\n\nif the first and last cluster start on the first/last sample, we dont know their real depth\n\nInputis assumed to be a thresholded Array with only 0/1\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.clusterdepth-Tuple{AbstractMatrix, Vararg{Any}}","page":"Home","title":"ClusterDepth.clusterdepth","text":"calculate clusterdepth of given datamatrix. \n\nclusterdepth([rng],data::AbstractMatrix,τ=2.3, statFun=twosidedstudentt,nperm=5000;pvaltype=:troendle)\n\n- `data`: `statFun` will be applied on second dimension of data (typically this will be subjects)\n\nOptional \t- τ: Cluster-forming threshold  \t- statFun: default  studenttest, can be any custom function on a Matrix returning a Vector \t- nperm: number of permutations, default 5000 \t- pval_type: how to calculate pvalues within each cluster, default :troendle, see ?pvals\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.ix_sortUnique-Tuple{Any}","page":"Home","title":"ClusterDepth.ix_sortUnique","text":"in some sense: argsort(argunique(x)), returns the indices to get a sorted unique of x\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.multicol_minimum-Tuple{AbstractMatrix, AbstractVector}","page":"Home","title":"ClusterDepth.multicol_minimum","text":"calculates the minimum in `X` along `dims=2` in the columns specified by àrrayOfIndicearrays` which could be e.g. `[[1,2],[5,6],[3,4,7]]`\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.pvals-Tuple{Any, ClusterDepth.ClusterDepthMatrix, Vararg{Any}}","page":"Home","title":"ClusterDepth.pvals","text":"Calculate pvals from cluster-depth permutation matrices\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.pvals-Tuple{Any}","page":"Home","title":"ClusterDepth.pvals","text":"pvals(data;kwargs...) = pvals(data[2:end],data[1];kwargs...)\n\npvals(data::AbstractVector,stat::Real;type=:twosided)\n\ncalculates pvalues based on permutation results\n\nif called with stat, first entry is assumed to be the observation \n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.sign_permute-Tuple{Any, AbstractMatrix, Any}","page":"Home","title":"ClusterDepth.sign_permute","text":"Permutation via random sign-flip Flips signs along the second dimension\n\n\n\n\n\n","category":"method"},{"location":"#ClusterDepth.troendle-Tuple{AbstractMatrix, AbstractVector}","page":"Home","title":"ClusterDepth.troendle","text":"function troendle(perm::AbstractMatrix,stat::AbstractVector;type=:twosided)\n\nMultiple Comparison Correction as in Troendle 1995\n\nperm with size  ntests x nperms\n\nstat with size ntests\n\ntype can be :twosided (default), :lesser, :greater\n\nHeavily inspired by the R implementation in permuco from Jaromil Frossard\n\nNote: While permuco is released under BSD, the author Jaromil Frossard gave us an MIT license for the troendle and the clusterdepth R-functions.\n\n\n\n\n\n","category":"method"}]
}
