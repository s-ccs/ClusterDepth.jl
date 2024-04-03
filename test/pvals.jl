@testset "pvals" begin
    @testset "pvals-onePerm" begin
        @test ClusterDepth.pvals(zeros(100), 1) == 1 / 101
        @test ClusterDepth.pvals(vcat(1, zeros(100))) == 1 / 101
        @test ClusterDepth.pvals(vcat(0, ones(100))) == 101 / 101
        @test ClusterDepth.pvals(vcat(80, 0:99)) == 21 / 101

    end
    @testset "pvals-ClusterDepthMatrix" begin
        cdm = ClusterDepth.ClusterDepthMatrix(sparse(ones(10, 1000)))
        p = ClusterDepth.pvals([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0], cdm, 0.1)

        @test p[7] ≈ 1 / 1001


        p = ClusterDepth.pvals([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0], cdm, 2.1)
        p .≈ 1000 / 1000



        p = ClusterDepth.pvals(
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0],
            cdm,
            0.1,
            type = :naive,
        )
        @test p[7] ≈ 1 / 1001


        J = zeros(10, 1000)
        J[1, :] .= 5
        J[2, :] .= 3
        J[3, :] .= 1
        cdm = ClusterDepth.ClusterDepthMatrix(sparse(J))

        p = ClusterDepth.pvals(
            [4.0, 0.0, 10.0, 0.0, 3.0, 4.0, 0.0, 4.0, 4.0, 0.0, 0.0],
            cdm,
            0.9,
        )
        @test all((p .> 0.05) .== [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1])

        #two tailed
        p = ClusterDepth.pvals(
            [4.0, 0.0, 10.0, 0.0, 3.0, 4.0, 0.0, 4.0, 4.0, 0.0, 0.0],
            (cdm, cdm),
            0.9,
        )
        @test all((p .> 0.05) .== [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1])


        p = ClusterDepth.pvals([0.0, 1.0, 2.0, 3.0, 4.0, 5.1, 2.0, 6.0], (cdm, cdm), 0.9)
        @test all((p .> 0.05) .== [1, 1, 1, 1, 1, 1, 1, 1])

        p = ClusterDepth.pvals([0.0, 1.0, 2.0, 3.0, 4.0, 5.1, 2.0, 0], (cdm, cdm), 0.9)
        @test all((p .> 0.05) .== [1, 1, 1, 0, 0, 0, 1, 1])


        p = ClusterDepth.pvals([0.0, 1.0, 2.0, 3.0, 4.0, 2.0, 1.0, 0.0, 6], (cdm, cdm), 0.9)
        @test all((p .> 0.05) .== [1, 1, 1, 0, 0, 1, 1, 1, 1])
    end


end
