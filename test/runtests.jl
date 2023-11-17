using ClusterDepth

include("setup.jl")
@testset "ClusterDepth.jl" begin

    include("troendle.jl")
    include("cluster.jl")
    include("utils.jl")
    include("pvals.jl")
end
