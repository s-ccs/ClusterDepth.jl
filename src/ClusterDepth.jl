module ClusterDepth

	using Random
	using Images
	using SparseArrays
	using StatsBase
	
    include("cluster.jl")
    include("pvals.jl")
    include("troendle.jl")
    
    export troendle
    clusterdepth(data::AbstractMatrix,args...;kwargs...) = clusterdepth(MersenneTwister(1),data,args...;kwargs...)
    function clusterdepth(rng,data::AbstractMatrix,τ=2.3, statFun=twosided_studentt,nperm=5000;pval_type=:troendle)
        (h,t) = perm_clusterdepths_both(rng,d,statFun,τ;nₚ=nperm)
        return pvals(statFun(d),(h,t);type=pval_type)
    end

end
