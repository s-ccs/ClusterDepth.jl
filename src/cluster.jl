"""

clusterdepth(rng,data::AbstractArray;τ=2.3, statFun=x->abs.(studentt(x)),permFun=sign_permute!,nperm=5000,pval_type=:troendle)

calculate clusterdepth of given datamatrix. 


- `data`: `statFun` will be applied on second dimension of data (typically this will be subjects)

Optional
- `τ`: Cluster-forming threshold 
- `statFun`: default  the one-sample `studenttest`, can be any custom function on a Matrix returning a Vector
- `permFun`: default to sign-flip (for one-sample case)
- `nperm`: number of permutations, default 5000
- `pval_type`: how to calculate pvalues within each cluster, default `:troendle`, see `?pvals`

"""
clusterdepth(data::AbstractArray,args...;kwargs...) = clusterdepth(MersenneTwister(1),data,args...;kwargs...)
function clusterdepth(rng,data::AbstractArray;τ=2.3, statFun=x->abs.(studentt(x)),permFun=sign_permute!,nperm=5000,pval_type=:troendle)
	cdmTuple = perm_clusterdepths_both(rng,data,statFun,permFun,τ;nₚ=nperm)
	return pvals(statFun(data),cdmTuple,τ;type=pval_type)
end




function perm_clusterdepths_both(rng,data,statFun,permFun,τ;nₚ=1000)
	
	#Jₖ_head = ExtendableSparseMatrix(size(data,2),nₚ)
	#Jₖ_tail = ExtendableSparseMatrix(size(data,2),nₚ)
	data_perm = deepcopy(data)
	rows_h = Int[]
	cols_h = Int[]
	vals_h = Float64[]
	rows_t = Int[]
	cols_t = Int[]
	vals_t = Float64[]
	for i = 1:nₚ
		# permute	
		d0 = permFun(rng,data_perm,statFun)

		# get clusterdepth
		(fromTo,head,tail) = calc_clusterdepth(d0,τ)
		
		# save it
		if !isempty(head)
			append!(rows_h,fromTo)
			append!(cols_h,fill(i,length(fromTo)))
			append!(vals_h,head)

			#Jₖ_head[fromTo,i] .+=head
		end
		if !isempty(tail)
			append!(rows_t,fromTo)
			append!(cols_t,fill(i,length(fromTo)))
			append!(vals_t,tail)
			#Jₖ_tail[fromTo,i] .+=tail
		end		
	end
	
	Jₖ_head = sparse(rows_h,cols_h,vals_h)#SparseMatrixCSC(nₚ,maximum(rows_h), cols_h,rows_h,vals_h)
	Jₖ_tail =sparse(rows_t,cols_t,vals_t)#SparseMatrixCSC(nₚ,maximum(rows_t), cols_t,rows_t,vals_t)
	return ClusterDepthMatrix((Jₖ_head)), ClusterDepthMatrix((Jₖ_tail))
end


"""
calc_clusterdepth(data,τ)

returns tuple with three entries:
1:maxLength, maximal clustervalue per clusterdepth head, same for tail

We assume data and τ have already been transformed for one/two sided testing, so that we can do d0.>τ for finding clusters

"""
function calc_clusterdepth(d0::AbstractArray{<:Real,2},τ)
	nchan = size(d0,1)

	# save all the results from calling calc_clusterdepth on individual channels
	(allFromTo, allHead, allTail) = (Array{Vector{Integer}}(undef, nchan),Array{Vector{Float64}}(undef, nchan),Array{Vector{Float64}}(undef, nchan))
	fromTo = []
	for i = 1: nchan
		(a,b,c) = calc_clusterdepth(d0[i,:],τ);
		allFromTo[i] = a
		allHead[i] = b
		allTail[i] = c
		
		if(length(a)>length(fromTo)) # running check to find the length ('fromTo') of the largest cluster
			fromTo = a
		end
	end

	# for each clusterdepth value, select the largest cluster value found across all channels
	(head, tail) = (zeros(length(fromTo)),zeros(length(fromTo)));
	for i in 1:nchan
		for j in allFromTo[i]
			if allHead[i][j] > head[j]
				head[j] = allHead[i][j]
			end
			if allTail[i][j] > tail[j]
				tail[j] = allTail[i][j]
			end
		end
	end
	
	return fromTo, head, tail
end

function calc_clusterdepth(d0,τ)
	startIX,len = cluster(d0.>τ)
	if isempty(len) # if nothing above threshold, just go on
		return [],[],[]
	end
	
	maxL = 1 + maximum(len) # go up only to max-depth
	
	valCol_head = Vector{Float64}(undef,maxL)
	valCol_tail = Vector{Float64}(undef,maxL)
	
	fromTo = 1:maxL
	for j = fromTo
		# go over clusters implicitly
		# select clusters that are larger (at least one)
		selIX = len.>=(j-1)
		

		
		if !isempty(selIX)
			ix= startIX[selIX] .+ (j - 1)
			valCol_head[j]= maximum(d0[ix])

			# potential optimization is that for j = 0 and maxL = 0, tail and head are identical	
			ix= startIX[selIX] .+ (len[selIX]) .- (j-1)
			valCol_tail[j]= maximum(d0[ix])
		end
	end
	return fromTo,valCol_head,valCol_tail
end

"""
finds neighbouring clusters in the vector and returns start + length vectors

if the first and last cluster start on the first/last sample, we dont know their real depth

	Input is assumed to be a thresholded Array with only 0/1
"""
function cluster(data)
	label = label_components(data) 
	K = maximum(label)
	start = fill(0, K)
	stop = fill(0, K)
	for k= 1:K
		#length[k] = sum(label.==k)
		start[k]  = findfirst(label.==k)
		stop[k]  = findlast(label.==k)
		
		
	end
	len = stop .- start

	# if the first and last cluster start on the first/last sample, we dont know their real depth
	if length(start) > 0 && start[end]+len[end] == length(data)
		start = start[1:end-1]
		len = len[1:end-1]
	end
	if length(start) > 0 && start[1] == 1
		start = start[2:end]
		len = len[2:end]
	end
	return start,len
end
