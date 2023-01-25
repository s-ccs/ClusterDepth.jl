"""
Takes a matrix filled with permutation results (results x permutations) and a observed matrix (size results) and calculates the p-values for each entry along 
the permutation dimension 
"""
function pvals_rankbased(perm::AbstractSparseMatrixCSC,
                            stat::AbstractVector;type=:twosided) 
    if length(stat) > size(perm,1) # larger cluster in stat than perms
        perm = sparse(findnz(perm)...,length(stat),size(perm,2))
    elseif length(stat) < size(perm,1) 
        #stat = sparsevec(findnz(stat)...,size(perm,1))
        i,j,v = findnz(perm)
        ix = i .<= length(stat)
        perm = sparse(i[ix],j[ix],v[ix],length(stat),size(perm,2))
    end
    
    
    # add stat to perm
    d = hcat(stat,perm)
    
    # fix specific testing
    if type==:twosided
        d = .-abs.(d)
    elseif type==:greater
        d = .-d
    end

    d = Matrix(d)

    # potential improvement, competerank 1224 -> but should be modified competerank 1334, then we could skip the expensive ceil below
    d = mapslices(tiedrank,d,dims=2)
    # rank & calc p-val
    #@show(d[1:10,1:10])
    d .=  ceil.(d)./(size(d,2))
    return d
    
end

"""
calculates pvalues based on permutation results
"""
pvals(data;kwargs...) = pvals(data[2:end],data[1];kwargs...)
function pvals(data::AbstractVector,stat::Real;type=:twosided)

    data = vcat(stat,data)
    if type == :greater || type  == :twosided
        comp = >=
        if type == :twosided
            data = abs.(data)
        end
    elseif type == :lesser
        comp = <=
    else
        error("not implemented")
    end

    pvals = sum(comp(stat[1]),data)/length(data) 
    return pvals
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
function multicol_minimum(x::AbstractMatrix,arrayOfIndicearrays::AbstractVector)
    min = fill(NaN,size(x,2),length(arrayOfIndicearrays))
    for to = 1:length(arrayOfIndicearrays)	
        @views min[:,to] =minimum(x[arrayOfIndicearrays[to],:],dims=1)
    end
    return min
end
	
"""
Multiple Comparison Correction as in Troendle 1995

Heavily inspired by the R implementation in permuco from Jaromil Frossard

"""
function troendle(perm::AbstractMatrix,stat::AbstractVector)
	pAll = pvals_rankbased(perm,stat)

	# get uncorrected pvalues of data
	#@show size(pAll)
	pD = pAll[:,1] # get first observation
	# rank the pvalues
	pD_rank = tiedrank(pD)

    # test in ascending order, same p-vals will be combined later
	testOrder = sort(unique(pD_rank))

	sortUn_ix = ix_sortUnique(pD_rank)
	testOrder = pD_rank[sortUn_ix]
	testOrder_all = [findall(x.==pD_rank) for x in pD_rank[sortUn_ix]]

	minPermPArray = multicol_minimum(pAll,testOrder_all)
	
	# the magic happens here, per permutation 
	resortPermPArray = similar(minPermPArray)

    # the following line can be made much faster by using views & reverse the arrays 
    #resortPermPArray = reverse(accumulate(min,reverse(minPermPArray),dims=2))
	@views accumulate!(min,resortPermPArray[end:-1:1],minPermPArray[end:-1:1],dims=2)
	
	pval_testOrder = pvals.(eachcol(resortPermPArray);type=:lesser)
	
	pval_testOrderMax = accumulate(max,pval_testOrder) # no idea why though
	
	uniqueToNonUnique = vcat([fill(x,length(v)) for (x,v) in enumerate(testOrder_all)]...)
	
	return pval_testOrderMax[uniqueToNonUnique][invperm(vcat(testOrder_all...))]
	
end

end