I got asked by a collaborator on how to call `ClusterDepth` from R for a **between groups** design.

You should first install JuliaCall+ClusterDepth

```R
install.packages("JuliaCall")
library("JuliaCall")
install_julia() # installs Julia
```

Then make a reproducible environment

```R
library("JuliaCall")

julia <- julia_setup() # setup julia
path_to_env = '/tmp/my_julia_env' # could be any path to a Julia Project.toml
julia_eval(paste('import Pkg;','Pkg.activate("',path_to_env,'");Pkg.instantiate()'))

# if Unfold is not yet installed
julia_eval('Pkg.add("Unfold");')
```

Now load the packages we need in julia

```R
library("JuliaCall")
julia_install_package_if_needed("ClusterDepth")
julia_library("ClusterDepth")
julia_library("Random")
```

we create some random data and a group indicator
```R
R_data <- array(matrix(rnorm(90), nrow = 10, ncol = 9),dim=c(1,10,9)) # 1 channel 10 timepoints, 9 subjects
grp_indicator = c("A","A","A","A","B","B","B","B","B") == "B" # should be boolean
```

... push the data to julia ...
```R
julia_assign("jl_data",R_data) 
julia_assign("grp",grp_indicator)
```

(for now) we need to create our own permutation function - given we simulated 3D data, it needs to be 3D (I plan to move this function to ClusterDepth at a later point).
```R
julia_eval("permute_fun(rng,x) = @view x[:,:,randperm(rng,size(x,3))]")
```

Finally, there is no good interface yet to tell `ClusterDepth` who belongs to what group for the stats test, so we have to "bake" it into the function.

```R
julia_eval("mystat = (x) -> ClusterDepth.studentt_unpaired(x,grp)")
```

The last thing to do is to actually call the algorithm with our defined `statfun` and `permfun`. Note the fixed cluster-forming-threshold, this should be adapted to your study, e.g. popular is to use α = 0.05 in the two-sided t-distribution: `quantile(TDist(n_subjects - 1), 0.985)` which requires the `Distributions.jl` package (apparently `qt(0.975, df = n_subjects - 1)` should do the trick in R as well).
```R
julia_eval("pvals = clusterdepth(jl_data;statfun=mystat,permfun = permute_fun,nperm=200,τ=2.3)")
```
