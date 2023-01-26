module ClusterDepth

	using Random
	using Images
	using SparseArrays
	using StatsBase
	import Base.show

    struct ClusterDepthMatrix{T} <: AbstractMatrix{T}
        J
        ClusterDepthMatrix(x) = new{SparseMatrixCSC}(sparse(x))

    end
    Base.show(io::IO, x::ClusterDepthMatrix) = show(io, x.J)

    Base.show(io::IO, m::MIME"text/plain", x::ClusterDepthMatrix) = show(io, m, x.J)

    struct test{T} <: AbstractArray{T,2}
        J
    end
    include("cluster.jl")
    include("pvals.jl")
    include("utils.jl")
    include("troendle.jl")
    
    export troendle
    export clusterdepth

end
