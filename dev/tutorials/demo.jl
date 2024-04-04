using ClusterDepth
using Random
using CairoMakie


n_t =40 # timepoints
n_sub = 50
n_perm = 5000

snr = 0.5 # signal to nois

## add a signal to the middle
signal = vcat(zeros(n_t÷4), sin.(range(0,π,length=n_t÷2)), zeros(n_t÷4))

## same signal for all subs
signal = repeat(signal,1,n_sub)


## add noise
data = randn(MersenneTwister(123),n_t,n_sub).+ snr .* signal

## by default assumes τ=2.3 (~alpha=0.05), and one-sample ttest
@time pvals = clusterdepth(data);

f = Figure()
ax = f[1,1] = Axis(f)


lines!(abs.(ClusterDepth.studentt(data)))
h1 = scatter(f[1,2],pvals;axis=(;yscale=log10),label="troendle")

pvals2 = clusterdepth(data;pval_type=:naive)
h2 = scatter!(1.2:40.2,pvals2,color="red",label="naive")
#hlines!(([0.05]))
axislegend()
f