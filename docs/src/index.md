```@meta
CurrentModule = ClusterDepth
```

# ClusterDepth

## Comparison to permuco R implementation
The implementation to Permuco is similar, but ClusterDepth.jl is more barebone - that is, we dont offer many permutation schemes, focus on the ClusterDepth Algorithm, and don't provide the nice wrappers like `clusterLM`.

Timing wise, a simple test on 50 subjects, 100 repetitions, 5000 permutations shows the following results:

|timepoints|ClusterDepth.jl|permuco|julia-speedup|
|---|---|---|---|
|40|0.03s|2.9s|~100x|
|400|0.14s|22s|~160x|
|4000|1.88s|240s|~120x|


```@index
```

```@autodocs
Modules = [ClusterDepth]
```
