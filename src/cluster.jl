"""
Permutation via random sign-flip
Flips signs along the second dimension
"""
function sign_permute(rng,x::AbstractMatrix,fun) 
	binflip=repeat(((rand(rng,size(x,2)).>0.5).*2).-1, 1,size(x,1))'
	return fun(x.*binflip)
end


function perm_clusterdepths_both(rng,data,statFun,τ;nₚ=1000)

	Jₖ_head = spzeros(m-1,nₚ)
	Jₖ_tail = spzeros(m-1,nₚ)
	
	for i = 1:nₚ
		# permute	
		d0 = sign_permute(rng,data,statFun)

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
	return shrink(Jₖ_head), shrink(Jₖ_tail)
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
