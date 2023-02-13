@testset "cluster" begin
    s,l = ClusterDepth.cluster([4., 0., 10., 0., 3., 4., 0, 4., 4.,0., 0., 5.].>0.9)
    @test s == [3,5,8]
    @test l == [0,1,1]

    s,l = ClusterDepth.cluster([0.,0.,0.,0.].>0.9)
    @test s == []
    @test l == []

    
end