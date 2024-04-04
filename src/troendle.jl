"""

pvals_rankbased(perm::AbstractMatrix,stat::AbstractVector;kwargs...)
pvals_rankbased(cdm::ClusterDepthMatrix,stat::AbstractSparseVector;kwargs...)

Takes a matrix filled with permutation results and a observed matrix (size ntests) and calculates the p-values for each entry along 
the permutation dimension 

`perm`: Matrix of permutations with dimensions `(ntests x permutations)`
`stat`: Vector of observed statistics size `ntests`


For `cdm::ClusterDepthMatrix` we can easily trim the permutation matrix towards the end (as it is actually a ragged Matrix).
That is, the permutation matrix might look like:

perm = [
    x x x x x x x;
    x x . x . x x;
    x x . x . . x;
    x x . . . . x;
    . x . . . . .;
    . . . . . . .;
    . . . . . . .;
    . . . . . . .;
]

Then the last three rows can simply be removed. rowwise-Ranks/pvalues would be identical anyway

The same will be checked for the stat-vector, if the stat vector is only "3" depths long, but the permutation has been calculated for "10" depths, we do not need to check the last 7 depths
of the permutation matrix.

**Output** will always be Dense Matrix (length(stat),nperm+1) with the first column being the pvalues of the observations
"""

function pvals_rankbased(cdm::ClusterDepthMatrix, stat::AbstractSparseVector; kwargs...)
    perm = cdm.J
    if length(stat) > size(perm, 1) # larger cluster in stat than perms
        perm = sparse(findnz(perm)..., length(stat), size(perm, 2))
    elseif length(stat) < size(perm, 1)
        #stat = sparsevec(findnz(stat)...,size(perm,1))
        i, j, v = findnz(perm)
        ix = i .<= length(stat)
        perm = sparse(i[ix], j[ix], v[ix], length(stat), size(perm, 2))
    end
    pvals = pvals_rankbased(perm, stat; kwargs...)
    return pvals
end
function pvals_rankbased(perm::AbstractMatrix, stat::AbstractVector; type = :twosided)

    # add stat to perm
    d = hcat(stat, perm)

    # fix specific testing
    if type == :twosided
        d = .-abs.(d)
    elseif type == :greater
        d = .-d
    elseif type == :lesser
    else
        error("unknown type")
    end

    d = Matrix(d)

    # potential improvement, competerank 1224 -> but should be modified competerank 1334, then we could skip the expensive ceil below
    d = mapslices(tiedrank, d, dims = 2)
    # rank & calc p-val
    #@show(d[1:10,1:10])
    d .= ceil.(d) ./ (size(d, 2))
    return d

end


"""
in some sense: `argsort(argunique(x))`, returns the indices to get a sorted unique of x
"""
function ix_sortUnique(x)
    uniqueidx(v) = unique(i -> v[i], eachindex(v))
    un_ix = uniqueidx(x)

    sort_ix = sortperm(x[un_ix])
    sortUn_ix = un_ix[sort_ix]
    return sortUn_ix
end


"""
    calculates the minimum in `X` along `dims=2` in the columns specified by àrrayOfIndicearrays` which could be e.g. `[[1,2],[5,6],[3,4,7]]`
"""
function multicol_minimum(x::AbstractMatrix, arrayOfIndicearrays::AbstractVector)
    min = fill(NaN, size(x, 2), length(arrayOfIndicearrays))
    for to = 1:length(arrayOfIndicearrays)
        @views min[:, to] = minimum(x[arrayOfIndicearrays[to], :], dims = 1)
    end
    return min
end

"""
function troendle(perm::AbstractMatrix,stat::AbstractVector;type=:twosided)

Multiple Comparison Correction as in Troendle 1995

`perm` with size  ntests x nperms

`stat` with size ntests

`type` can be :twosided (default), :lesser, :greater

Heavily inspired by the R implementation in permuco from Jaromil Frossard

Note: While permuco is released under BSD, the author Jaromil Frossard gave us an MIT license for the troendle and the clusterdepth R-functions.


"""
function troendle(perm::AbstractMatrix, stat::AbstractVector; type = :twosided)
    pAll = pvals_rankbased(perm, stat; type = type)

    # get uncorrected pvalues of data
    #@show size(pAll)
    pD = pAll[:, 1] # get first observation
    # rank the pvalues
    pD_rank = tiedrank(pD)

    # test in ascending order, same p-vals will be combined later

    # the following two lines are 
    sortUn_ix = ix_sortUnique(pD_rank)

    # these two lines would be identical
    #testOrder = pD_rank[sortUn_ix]
    #testOrder = sort(unique(pD_rank))

    # as ranks can be tied, we have to take those columns, and run a "min" on them
    testOrder_all = [findall(x .== pD_rank) for x in pD_rank[sortUn_ix]]

    minPermPArray = multicol_minimum(pAll, testOrder_all)

    # the magic happens here, per permutation 
    resortPermPArray = similar(minPermPArray)

    # the following line can be made much faster by using views & reverse the arrays 
    #resortPermPArray = reverse(accumulate(min,reverse(minPermPArray),dims=2))
    @views accumulate!(
        min,
        resortPermPArray[:, end:-1:1],
        minPermPArray[:, end:-1:1],
        dims = 2,
    )

    pval_testOrder = pvals.(eachcol(resortPermPArray); type = :lesser)

    pval_testOrderMax = accumulate(max, pval_testOrder) # no idea why though

    uniqueToNonUnique = vcat([fill(x, length(v)) for (x, v) in enumerate(testOrder_all)]...)

    return pval_testOrderMax[uniqueToNonUnique][invperm(vcat(testOrder_all...))]

end
