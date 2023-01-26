@testset "sign_permute" begin

    m = [1 1 1;2 2 2;3 3 3;4 4 4]
    p = ClusterDepth.sign_permute(StableRNG(1),m,x->x)

    @test p[1,:] == [1, -1, 1]

    # different seeds are different
    @test p!= ClusterDepth.sign_permute(StableRNG(3),m,x->x)
    # same seeds are the same
    @test p == ClusterDepth.sign_permute(StableRNG(1),m,x->x)

    m = ones(1,10000)
    @test abs.(ClusterDepth.sign_permute(StableRNG(1),m,mean))<0.001
    
end

@testset "studentt" begin
    x = randn(StableRNG(1),10000,50)
    t = ClusterDepth.studentt(x)
    @test length(t) == 10000
    @test maximum(abs.(t))<10 # we'd need to be super lucky ;)
    @test mean(abs.(t).>2) < 0.06

    
end