"""
calculate clusterdepth of given datamatrix. 

clusterdepth([rng],data::AbstractMatrix,τ=2.3, statFun=twosided_studentt,nperm=5000;pval_type=:troendle)

	- `data`: `statFun` will be applied on second dimension of data (typically this will be subjects)

Optional
	- `τ`: Cluster-forming threshold 
	- `statFun`: default  `studenttest`, can be any custom function on a Matrix returning a Vector
	- `nperm`: number of permutations, default 5000
	- `pval_type`: how to calculate pvalues within each cluster, default `:troendle`, see `?pvals`

"""
clusterdepth(data::AbstractMatrix,args...;kwargs...) = clusterdepth(MersenneTwister(1),data,args...;kwargs...)
function clusterdepth(rng,data::AbstractMatrix;τ=2.3, statFun=x->abs.(studentt(x)),permFun=sign_permute,nperm=5000,pval_type=:troendle)
	cdmTuple = perm_clusterdepths_both(rng,data,statFun,permFun,τ;nₚ=nperm)
	return pvals(statFun(data),cdmTuple,τ;type=pval_type)
end




function perm_clusterdepths_both(rng,data,statFun,permFun,τ;nₚ=1000)

	Jₖ_head = spzeros(size(data,1),nₚ)
	Jₖ_tail = spzeros(size(data,1),nₚ)
	
	for i = 1:nₚ
		# permute	
		d0 = permFun(rng,data,statFun)

		# get clusterdepth
		(fromTo,head,tail) = calc_clusterdepth(d0,τ)
		
		# save it
		Jₖ_head[fromTo,i]=head
		Jₖ_tail[fromTo,i]=tail
		
	end
	# shrink J_k

		function shrink(x) 
			i,j,v = findnz(x)
			return sparse(i,j,v,maximum(i),nₚ)
		end
	#return all of it
	return ClusterDepthMatrix(shrink(Jₖ_head)), ClusterDepthMatrix(shrink(Jₖ_tail))
end



function calc_clusterdepth(d0,τ)
	startIX,len = cluster(d0,τ)
	if isempty(len) # if nothing above threshold, just go on
		return [],[],[]
	end
	
	maxL = maximum(len) # go up only to max-depth
	
	valCol_head = Vector{Float64}(undef,maxL)
	valCol_tail = Vector{Float64}(undef,maxL)
	
	fromTo = 1:maxL
	for j = fromTo
		# go over clusters implicitly
		# select clusters that are larger (at least one)
		selIX = len.>=j
		
		if !isempty(selIX)
			ix= startIX[selIX] .+ j
			valCol_head[j]= maximum(d0[ix])
			
			ix= len[selIX] .+ startIX[selIX] .- j
			valCol_tail[j]= maximum(d0[ix])
		end
	end
	return fromTo,valCol_head,valCol_tail
end

"""
finds neighbouring clusters in the vector and returns start + length vectors
"""
function cluster(data,τ)
	label = label_components(data.>τ) 
	K = maximum(label)
	start = fill(0, K)
	stop = fill(0, K)
	for k= 1:K
		#length[k] = sum(label.==k)
		start[k]  = findfirst(label.==k)
		stop[k]  = findlast(label.==k)
		
		
	end
	len = stop .- start
	
	return start,len
end
