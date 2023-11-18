@testset "sign_permute" begin

    m = [1 1 1;2 2 2;3 3 3;4 4 4]
    p = ClusterDepth.sign_permute!(StableRNG(1),m,x->x)
    
    @test p[1,:] == [1, -1, -1]

    # different seeds are different
    @test p!= ClusterDepth.sign_permute!(StableRNG(3),deepcopy(m),x->x)
    # same seeds are the same
    @test p == ClusterDepth.sign_permute!(StableRNG(1),deepcopy(m),x->x)

    m = ones(1,1000000)
    @test abs.(ClusterDepth.sign_permute!(StableRNG(1),deepcopy(m),mean))<0.001

    m = ones(1,2,3,4,5,6,7,100)
    o = ClusterDepth.sign_permute!(StableRNG(1),deepcopy(m),x->x)
    @test sort(unique(mean(o,dims=1:ndims(o)-1))) == [-1.,1.]
    
    #2D input data
    data = randn(StableRNG(1),4,5)
    @test size(ClusterDepth.sign_permute!(StableRNG(1),data,ClusterDepth.studentt)) == (4,)

    #3D input data
    data = randn(StableRNG(1),3,4,5);
    @test size(ClusterDepth.sign_permute!(StableRNG(1),data,ClusterDepth.studentt)) == (3,4) 
    
end

@testset "studentt" begin
    x = randn(StableRNG(1),10000,50)
    t = ClusterDepth.studentt(x)
    @test length(t) == 10000
    @test maximum(abs.(t))<10 # we'd need to be super lucky ;)
    @test mean(abs.(t).>2) < 0.06

    #2D input data
    data = randn(StableRNG(1),4,5)
    @test size(ClusterDepth.studentt(data)) == (4,)

    #3D input data
    data = randn(StableRNG(1),3,4,5);
    @test size(ClusterDepth.studentt(data)) == (3,4) 
end
