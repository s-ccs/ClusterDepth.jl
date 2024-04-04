@testset "cluster" begin
    s, l = ClusterDepth.cluster(
        [4.0, 0.0, 10.0, 0.0, 3.0, 4.0, 0, 4.0, 4.0, 0.0, 0.0, 5.0] .> 0.9,
    )
    @test s == [3, 5, 8]
    @test l == [0, 1, 1]

    s, l = ClusterDepth.cluster([0.0, 0.0, 0.0, 0.0] .> 0.9)
    @test s == []
    @test l == []
end

@testset "Tests for 2D data" begin
    data = randn(StableRNG(1), 4, 5)
    @show ClusterDepth.calc_clusterdepth(data, 0)
end

@testset "Tests for 3D data" begin
    data = randn(StableRNG(1), 3, 20, 5)
    @show ClusterDepth.clusterdepth(data; τ = 0.4, nperm = 5)
end
