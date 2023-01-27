```@meta
CurrentModule = ClusterDepth
```

# ClusterDepth

Documentation for [ClusterDepth](https://github.com/behinger/ClusterDepth.jl).

## Comparison to permuco R implementation
The implementation to Permuco is similar, but ClusterDepth.jl is more barebone - that is, we dont offer many permutation schemes, focus on the ClusterDepth Algorithm, and don't provide the nice wrappers like `clusterLM`.

Timing wise, a simple test on 50 subjects, 100 repetitions, 5000 permutations shows the following results:

|timepoints|ClusterDepth.jl|permuco|julia-speedup|
|---|---|---|---|
|40|0.1s|2.9s|29x|
|400|0.6s|22s|36x|
|4000|7s|240s|34x|


```@index
```

```@autodocs
Modules = [ClusterDepth]
```
