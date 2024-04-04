function studentt!(out::AbstractMatrix, x::AbstractArray{<:Real,3}; kwargs...)

    for (x_ch, o_ch) in zip(eachslice(x, dims = 1), eachslice(out, dims = 1))
        #@debug size(x_ch),size(o_ch)
        studentt!(o_ch, x_ch; kwargs...)
    end
    return out
end


"""
	studentt_test!(out,x;type=:abs)

strongly optimized one-sample t-test function.


Implements: t =  mean(x) / ( sqrt(var(x))) / sqrt(size(x,2)-1)

Accepts 2D or 3D matrices, always aggregates over the last dimension

"""
function studentt_test!(out, x::AbstractMatrix)
    mean!(out, x)
    df = 1 ./ sqrt(size(x, 2) - 1)
    #@debug size(out),size(x)
    tmp = [1.0]
    for k in eachindex(out)
        std!(tmp, @view(x[k, :]), out[k])
        @views out[k] /= (sqrt(tmp[1]) * df)
    end
    return out
end

function std!(tmp, x_slice, μ)

    @views x_slice .= (x_slice .- μ) .^ 2
    sum!(tmp, x_slice)
    tmp .= sqrt.(tmp ./ (length(x_slice) - 1))


end


function studentt!(out, x)
    #@debug size(out),size(x)
    mean!(out, x)
    out .= out ./ (std(x, mean = out, dims = 2)[:, 1] ./ sqrt(size(x, 2) - 1))
end
function studentt(x::AbstractMatrix)
    # more efficient than this one liner
    # studentt(x::AbstractMatrix) = (mean(x,dims=2)[:,1])./(std(x,dims=2)[:,1]./sqrt(size(x,2)-1))
    μ = mean(x, dims = 2)[:, 1]
    μ .= μ ./ (std(x, mean = μ, dims = 2)[:, 1] ./ sqrt(size(x, 2) - 1))
end

studentt(x::AbstractArray{<:Real,3}) =
    dropdims(mapslices(studentt, x, dims = (2, 3)), dims = 3)

"""
Permutation via random sign-flip
Flips signs along the last dimension
"""
function sign_permute!(rng, x::AbstractArray)
    n = ndims(x)
    @assert n > 1 "vectors cannot be permuted"

    fl = rand(rng, [-1, 1], size(x, n))

    for (flip, xslice) in zip(fl, eachslice(x; dims = n))
        xslice .= xslice .* flip
    end

    return x
end
