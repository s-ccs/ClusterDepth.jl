@testset "troendle" begin
    @testset "troendle()" begin
        nperm = 900
        ntests = 9
        rng = StableRNG(1)
        perm = randn(rng, ntests, nperm)
        stat = [0.0, 0.0, -1, 1, -3.5, 3, 10, 10, -10]

        p_cor = troendle(perm, stat)

        @test length(p_cor) == ntests
        @test all(p_cor[1:2] .== 1.0)
        @test all(p_cor[7:9] .< 0.008)
        #
        statSig = [100.0, 100, 100, 100, 100, 100, 100, 100, 100]
        p_cor = troendle(perm, statSig)
        @test length(unique(p_cor)) == 1
        @test all(p_cor .== (1 / (1 + nperm)))



        # sidedness
        h_cor = troendle(perm, stat; type = :greater) .< 0.05
        @test h_cor == [0, 0, 0, 0, 0, 1, 1, 1, 0]
        h_cor = troendle(perm, stat; type = :lesser) .< 0.05
        @test h_cor == [0, 0, 0, 0, 1, 0, 0, 0, 1]
    end

    @testset "pvals_rankbased" begin
        nperm = 9
        ntests = 5
        perm = randn(StableRNG(1), ntests, nperm)
        perm[4:5, :] .= 0
        J = sparse(perm)

        stat_short = sparsevec(1:2, [3, 3])
        stat_longer = sparsevec(1:2, [3, 3], ntests + 2)
        stat = sparsevec(1:2, [3, 3], ntests)
        cdm = ClusterDepth.ClusterDepthMatrix(J)

        p = ClusterDepth.pvals_rankbased(perm, stat)
        @test ClusterDepth.pvals_rankbased(cdm, stat) == p

        # trimming shorter vec
        p_short = ClusterDepth.pvals_rankbased(cdm, stat_short)
        @test size(p_short) == (2, nperm + 1)
        @test p_short == p[1:2, :]

        # extending longer vec
        p_long = ClusterDepth.pvals_rankbased(cdm, stat_longer)
        @test size(p_long) == (ntests + 2, nperm + 1)

        p = ClusterDepth.pvals_rankbased([1 2 3 4 5; 1 2 3 4 5], [3.5, 3])
        @test p[:, 1] ≈ [0.5, 2.0 / 3]
        @test p[1, 2:end] ≈ p[2, 2:end]

        p = ClusterDepth.pvals_rankbased([1 2 3 4 5; 1 2 3 4 5], [0, 6])
        @test p[:, 1] ≈ [1.0, 1.0 / 6]


        p = ClusterDepth.pvals_rankbased(
            ClusterDepth.ClusterDepthMatrix([1 2 3 4 5; 1 2 3 4 5]),
            sparse([-6, -6, -6, -6, -6]),
        )

        @test all(p[:, 1] .== [1 / 6.0])

        # test sidedness
        p = ClusterDepth.pvals_rankbased([1 2 3 4 5; 1 2 3 4 5], [-6, 6])
        @test p[:, 1] ≈ [1.0 / 6, 1.0 / 6]

        p = ClusterDepth.pvals_rankbased([1 2 3 4 5; 1 2 3 4 5], [-6, 6]; type = :lesser)
        @test p[:, 1] ≈ [1.0 / 6, 1.0]

        p = ClusterDepth.pvals_rankbased([1 2 3 4 5; 1 2 3 4 5], [-6, 6]; type = :greater)
        @test p[:, 1] ≈ [1.0, 1.0 / 6]

    end



end
