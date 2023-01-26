using ClusterDepth
using Random
using CairoMakie


n_t =20 # timepoints
n_sub = 50
n_perm = 5000

snr = 0.5 # signal to nois

# add a signal to the middle
signal = vcat(zeros(n_t÷4), sin.(range(0,π,length=n_t÷2)), zeros(n_t÷4))

# same signal for all subs
signal = repeat(signal,1,n_sub)


# add noise
data = randn(n_t,n_sub).+ snr .* signal
pvals = clusterdepth(data)

f = Figure()
ax = f[1,1] = Axis(f)
# by default assumes τ=2.3 (~alpha=0.05), and one-sample ttest


lines!(ClusterDepth.studentt(data))
scatter(f[1,2],pvals;axis=(;yscale=log10))

pvals2 = clusterdepth(data;pval_type=:naive)
scatter!(pvals2,color="red")
hlines!([0.05,0.01])
f