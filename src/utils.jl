
studentt(x::AbstractMatrix) = (mean(x,dims=2)[:,1])./(std(x,dims=2)[:,1]./sqrt(size(x,2)-1))


"""
Permutation via random sign-flip
Flips signs along the second dimension
"""
function sign_permute(rng,x::AbstractMatrix,fun) 
	binflip=repeat(((rand(rng,size(x,2)).>0.5).*2).-1, 1,size(x,1))'
	return fun(x.*binflip)
end
