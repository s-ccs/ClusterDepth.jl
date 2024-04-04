using ClusterDepth
using Random
using CairoMakie
using UnfoldSim
using StatsBase
using ProgressMeter
using Distributions

# ### Family Wise Error of Troendle
# Here we calculate the Family Wise Error of doing `ntests` at the same time.
# That is, we want to check that Troendle indeed returns us a type-1 of 5% for a _set_ of tests.
#
# The point being, that if you do 30 tests, the chance that one is significant is not 5% but actually 
(1 - (1 - 0.05)^30) * 100 ##%

# Let's setup some simulation parameters
reps = 1000
perms = 1000
ntests = 30;

# we will use the student-t in it's 2-sided variant (abs of it)
fun = x -> abs.(ClusterDepth.studentt(x));

# this function simulates data without any effect (H0), then the permutations, and finally calls troendle
function run_fun(r, perms, fun, ntests)
    rng = MersenneTwister(r)
    data = randn(rng, ntests, 50)
    perm = Matrix{Float64}(undef, size(data, 1), perms)
    stat = fun(data)
    for p = 1:perms
        ClusterDepth.sign_permute!(rng, data)
        perm[:, p] = fun(data)
    end
    return data, stat, troendle(perm, stat)
end;
# let's test it once
data, stats_t, pvals = run_fun(1, perms, fun, ntests);
println("data:", size(data), " t-stats:", size(stats_t), " pvals:", size(pvals))


# run the above function `reps=1000`` times - we also save the uncorrected t-based pvalue
pvals_all = fill(NaN, reps, 2, ntests)
Threads.@threads for r = 1:reps
    data, stat, pvals = run_fun(r, perms, fun, ntests)
    pvals_all[r, 1, :] = pvals
    pvals_all[r, 2, :] = (1 .- cdf.(TDist(size(data, 2)), abs.(stat))) .* 2 # * 2 becaue of twosided. Troendle takes this into account already
end;

# Let's check in how many of our simlations we have a significant p-value =<0.05 
res = any(pvals_all[:, :, :] .<= 0.05, dims=3)[:, :, 1]
mean(res .> 0, dims=1) |> x -> (:troendle => x[1], :uncorrected => x[2])

# Nice. Troendle fits perfectly and the uncorrected is pretty close to what we calculated above!

# Finally we end this with a short figure to get a better idea of how this data looks like and a histogram of the p-values

f = Figure()
ax = f[1, 1] = Axis(f)



lines!(ax, abs.(ClusterDepth.studentt(data)))
heatmap!(Axis(f[2, 1]), data)
series!(Axis(f[2, 2]), data[:, 1:7]')
h1 = scatter!(Axis(f[1, 2]; yscale=log10), pvals, label="troendle")

hlines!([0.05, 0.01])

hist!(Axis(f[3, 1]), pvals_all[:, 1, :][:], bins=0:0.01:1.1)
hist!(Axis(f[3, 2]), pvals_all[:, 2, :][:], bins=0:0.01:1.1)
f

