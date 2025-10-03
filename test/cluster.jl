

@testset "cluster" begin
    s, l = ClusterDepth.cluster(
        [4.0, 0.0, 10.0, 0.0, 3.0, 4.0, 0, 4.0, 4.0, 0.0, 0.0, 5.0] .> 0.9,
    )

    @test s == [3, 5, 8]
    @test l == [0, 1, 1]

    s, l = ClusterDepth.cluster([0.0, 0.0, 0.0, 0.0] .> 0.9)
    @test s == []
    @test l == []


    s, l = ClusterDepth.cluster(
        [4.0, 0.0, 10.0, 0.0, 3.0, 4.0, 0, 4.0, 4.0] .> 0.9,
    )
    @test s == [3, 5]
    @test l == [0, 1]

    s, l = ClusterDepth.cluster(
        [4.0, 0.0, 10.0, 0.0, 3.0, 4.0, 0, 4.0, 4.0] .> 0.9,
    )
end

@testset "Tests for 2D data" begin
    data = randn(StableRNG(1), 4, 5)
    res = ClusterDepth.clusterdepth(data; τ=0.4)
    @test size(res) == (4,)
end

@testset "Tests for 3D data" begin
    data = randn(StableRNG(1), 3, 20, 5)
    res = ClusterDepth.clusterdepth(data; τ=0.4, nperm=5)
    @test size(res) == (3, 20)
end
@testset "Test sidefun" begin
    data = randn(StableRNG(1), 23, 20)
    data[3:8, :] .+= 3
    data[12:17, :] .-= 3
    res = ClusterDepth.clusterdepth(data; τ=0.4, nperm=5)
    res_negated = ClusterDepth.clusterdepth(.-data; τ=0.4, nperm=5)
    @test res ≈ res_negated # should be same if side_type=:abs


    # testing the default is abs
    res_abs = ClusterDepth.clusterdepth(data; τ=0.4, nperm=5, side_type=:abs)
    @test res ≈ res_abs


    res_pos = ClusterDepth.clusterdepth(data; τ=0.4, nperm=5, side_type=:positive)
    @test all(res_pos[3:8] .< 0.8)
    res_neg = ClusterDepth.clusterdepth(data; τ=0.4, nperm=5, side_type=:negative)
    @test all(res_neg[12:17] .< 0.8)


end
@testset "Test warning clusterbegin/end" begin

    # test the warning that clusters must not begin/end with a potentially significant cluster
    data = randn(StableRNG(1), 23, 20)
    data[1:5, :] .+= 3
    @test_warn x -> occursin("Your data shows a cluster", x) ClusterDepth.clusterdepth(data; τ=0.4, nperm=5)
    data = randn(StableRNG(1), 23, 20)
    data[23, :] .-= 3
    @test_warn x -> occursin("Your data shows a cluster", x) ClusterDepth.clusterdepth(data; τ=0.4, nperm=5)

end

@testset "sparse matrix size" begin
    data = rand(StableRNG(42), 5, 10) .- 0.5
    data[2, :] .+= 1
    data[3, :] .+= 3
    cdmTuple = ClusterDepth.perm_clusterdepths_both(StableRNG(1), data, (ClusterDepth.sign_permute!), 1.5; nₚ=1000, (statfun!)=(ClusterDepth.studentt!), sidefun=abs)
    @test size(cdmTuple[1].J, 2) == 1000 # sparse matrix has to have the size of n-perms
end