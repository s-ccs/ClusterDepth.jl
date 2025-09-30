function studentt!(out::AbstractMatrix, x::AbstractArray{<:Real,3}; kwargs...)

    for (x_ch, o_ch) in zip(eachslice(x, dims=1), eachslice(out, dims=1))
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
    out .= out ./ (std(x, mean=out, dims=2, corrected=true)[:, 1] ./ sqrt(size(x, 2)))
end
function studentt(x::AbstractMatrix)
    # more efficient than this one liner
    # studentt(x::AbstractMatrix) = (mean(x,dims=2)[:,1])./(std(x,dims=2)[:,1]./sqrt(size(x,2)-1))
    μ = mean(x, dims=2)[:, 1]
    μ .= μ ./ (std(x, mean=μ, dims=2, corrected=true)[:, 1] ./ sqrt(size(x, 2)))
end

studentt(x::AbstractArray{<:Real,3}) =
    dropdims(mapslices(studentt, x, dims=(2, 3)), dims=3)

"""
Permutation via random sign-flip
Flips signs along the last dimension
"""
function sign_permute!(rng, x::AbstractArray)
    n = ndims(x)
    @assert n > 1 "vectors cannot be permuted"

    fl = rand(rng, [-1, 1], size(x, n))

    for (flip, xslice) in zip(fl, eachslice(x; dims=n))
        xslice .= xslice .* flip
    end

    return x
end

"""
    studentt_unpaired(x::AbstractArray, group)

Implements a unpaired two groups t-test with unequal variances. 10x as fast as HypothesisTests because we don't allocate that much
Use like this:

```julia
studentt_unpaired([1,2,1,1,4,5,4],[false,false,false,false,true,true,true])
```

To use with ClusterDepth, you have to "bake-in" the group membership by defining your own method:
```julia
grp = ["bla","bla","bla","bla","blub","blub","blub"] .== "blub"
my_statfun = x->studentt_unpaired(x,grp)
```

"""
function studentt_unpaired(x, group)
    x_reshaped = reshape(x, :, size(x, ndims(x)))
    x₁ = x_reshaped[:, group]
    x₂ = x_reshaped[:, .!group]

    n₁, n₂ = size(x₁, 2), size(x₂, 2)
    μ₁, μ₂ = mean(x₁, dims=2), mean(x₂, dims=2)
    var₁, var₂ = var(x₁, dims=2, corrected=true), var(x₂, dims=2, corrected=true)

    se = sqrt.(var₁ ./ n₁ .+ var₂ ./ n₂)
    t_stat = (μ₁ .- μ₂) ./ se

    # Reshape t_stat to match input dimensions (excluding last dim)
    t_stat_reshaped = reshape(t_stat, size(x)[1:end-1])
    return t_stat_reshaped
end