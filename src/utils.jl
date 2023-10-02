
studentt(x::AbstractMatrix) = (mean(x,dims=2)[:,1])./(std(x,dims=2)[:,1]./sqrt(size(x,2)-1))

studentt(x::AbstractArray{<:Real,3}) = return dropdims(mapslices(studentt,x,dims=(2,3)),dims=3)

"""
Permutation via random sign-flip
Flips signs along the last dimension
"""
function sign_permute(rng,x::AbstractArray,fun) 
	n = ndims(x)
	@assert n > 1 "vectors cannot be permuted"
	
	fl = rand(rng,[-1,1],size(x,n))
	flipped = map((x,y)->x.*y,eachslice(x;dims=n),fl)

	y = similar(x)
	for (ix,sl) = enumerate(eachslice(y;dims=n))
		sl.=flipped[ix]
	end 

	return fun(y)
end
