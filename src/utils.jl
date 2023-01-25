
twosided_studentt(x::AbstractMatrix) = abs.(mean(x,dims=2)[:,1])./(std(x,dims=2)[:,1]./sqrt(size(x,2)-1))