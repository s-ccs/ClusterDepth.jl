"""
permute data and returns clusterdepth

returns for each cluster-depth, the maximum stat-value over all detected clusters, from both sides at the same time
    
Conversely, perm_clusterdepths_both, returns the maximum stat-value from head and tail separately (combined after calcuating pvalues)
"""
function perm_clusterdepths_combined(rng,data,statFun)

	Jₖ = spzeros(m-1,nₚ)
	
		
	for i = 1:nₚ
			
			d0 = perm(rng,data,statFun)
			startIX,len = get_clusterdepths(d0,τ)
			maxL = maximum(len)
			# go up to max-depth
			valCol = Vector{Float64}(undef,maxL)
			
			fromTo = 1:maxL
			for j = fromTo

				# go over clusters implicitly
				# select clusters that are larger (at least one)
				selIX = len.>=j
				
				if !isempty(selIX)
					ix_head= startIX[selIX] .+ j
					ix_tail= len[selIX] .+ startIX[selIX] .- j
					
					valCol[j]= maximum(d0[vcat(ix_head,ix_tail)])
					
				end
			end
			Jₖ[fromTo,i] = valCol
			
		end
		return Jₖ
end;