using HypothesisTests
@testset "sign_permute" begin

    m = [1 1 1; 2 2 2; 3 3 3; 4 4 4]
    p = ClusterDepth.sign_permute!(StableRNG(2), deepcopy(m))

    @test p[1, :] == [1, -1, 1]

    # different seeds are different
    @test p != ClusterDepth.sign_permute!(StableRNG(3), deepcopy(m))
    # same seeds are the same
    @test p == ClusterDepth.sign_permute!(StableRNG(2), deepcopy(m))

    m = ones(1, 1000000)
    @test abs(mean(ClusterDepth.sign_permute!(StableRNG(1), deepcopy(m)))) < 0.001

    m = ones(1, 2, 3, 4, 5, 6, 7, 100)
    o = ClusterDepth.sign_permute!(StableRNG(1), deepcopy(m))
    @test sort(unique(mean(o, dims=1:ndims(o)-1))) == [-1.0, 1.0]


end

@testset "studentt" begin
    x = randn(StableRNG(1), 10000, 50)
    t = ClusterDepth.studentt(x)

    t_true = [HypothesisTests.OneSampleTTest(r).t for r in eachrow(x)]
    @test all(t .≈ t_true)
    @test length(t) == 10000
    @test maximum(abs.(t)) < 10 # we'd need to be super lucky ;)
    @test mean(abs.(t) .> 2) < 0.06

    #2D input data
    data = randn(StableRNG(1), 4, 5)
    @test size(ClusterDepth.studentt(data)) == (4,)

    #3D input data
    data = randn(StableRNG(1), 3, 4, 5)
    @test size(ClusterDepth.studentt(data)) == (3, 4)

    #
    t = rand(10000)
    ClusterDepth.studentt!(t, x)
    @test t ≈ ClusterDepth.studentt(x)

    @test length(t) == 10000
    @test maximum(abs.(t)) < 10 # we'd need to be super lucky ;)
    @test mean(abs.(t) .> 2) < 0.06

    #2D input data
    data = randn(StableRNG(1), 4, 5)
    t = rand(4)
    ClusterDepth.studentt!(t, data)
    @test size(t) == (4,)

    #3D input data
    data = randn(StableRNG(1), 3, 4, 5)
    @test size(ClusterDepth.studentt(data)) == (3, 4)
end


@testset begin
    "ttest unpaired"
    _x = randn(StableRNG(1), 10000, 50, 2)
    x1 = _x[:, :, 1]
    x2 = _x[:, :, 2]
    x = hcat(x1, x2)
    group = repeat([false, true], inner=size(x1, 2))
    @benchmark t = ClusterDepth.studentt_unpaired(x, group)

    @benchmark t_true = [HypothesisTests.UnequalVarianceTTest(r[group], r[.!group]).t for r in eachrow(x)]
    @test all(t .≈ t_true)
    @test length(t) == 10000


    _x = randn(StableRNG(1), 10000, 50, 30, 2)
    x1 = _x[:, :, :, 1]
    x2 = _x[:, :, :, 2]
    x = cat(x1, x2, dims=3)
    group = repeat([false, true], inner=size(x1, 3))
    t = ClusterDepth.studentt_unpaired(x, group .== 1)
    @test size(t) == (10000, 50)
end