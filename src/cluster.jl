"""

using Base: Stateful
clusterdepth(rng,data::AbstractArray;τ=2.3, statfun=x->abs.(studentt(x)),permfun=sign_permute!,nperm=5000,pval_type=:troendle)

calculate clusterdepth of given datamatrix. 


- `data`: `statfun` will be applied on last dimension of data (typically this will be subjects)

Optional
- `τ`: Cluster-forming threshold 
- `nperm`: number of permutations, default 5000
- `stat_type`: default  the one-sample `t-test`, custom function can be specified (see `statfun!` and `statfun`)
- `side_type`: default: `:abs` - what function should be applied after the `statfun`? could be `:abs`, `:square`, `:positive` to test positive clusters, `:negative` to test negative clusters. Custom function can be provided, see `sidefun``
- `perm_type`: default `:sign` for one-sample data (e.g. differences), performs sign flips. custom function can be provided, see  `permfun`
- `pval_type`: how to calculate pvalues within each cluster, default `:troendle`, see `?pvals`
- `statfun` / `statfun!` a function that either takes one or two arguments and aggregates over last dimension. in the two argument case we expect the first argument to be modified inplace and provide a suitable Vector/Matrix.
- `sidefun`: default `abs`. Provide a function to be applied on each element of the output of  `statfun`. 
- `permfun` function to permute the data, should accept an RNG-object and the data. can be inplace, the data is copied, but the same array is shared between permutations
- `show_warnings`: default `true` - function to suppress warnings, useful for simulations
"""
clusterdepth(data::AbstractArray, args...; kwargs...) =
    clusterdepth(MersenneTwister(1), data, args...; kwargs...)
function clusterdepth(
    rng,
    data::AbstractArray;
    τ=2.3,
    stat_type=:onesample_ttest,
    perm_type=:sign,
    side_type=:abs,
    nperm=5000,
    pval_type=:troendle,
    (statfun!)=nothing,
    statfun=nothing,
    permfun=nothing,
    show_warnings=true,
)
    if stat_type == :onesample_ttest && isnothing(statfun!) && isnothing(statfun)
        statfun! = studentt!
        statfun = studentt
    end
    if perm_type == :sign
        if isnothing(permfun)
            permfun = sign_permute!
        end
    end
    if side_type == :abs
        sidefun = abs
    elseif side_type == :square
        sidefun = x -> x^2
    elseif side_type == :negative
        sidefun = x -> -x
    elseif side_type == :positive
        sidefun = x -> x # the default :)
    else
        @assert isnothing(side_type) "unknown side_type ($side_type) specified. Check your spelling and ?clusterdepth"
    end
    data_obs = sidefun.(statfun(data))

    if show_warnings && (any(data_obs[:, 1] .> τ) || any(data_obs[:, end] .> τ))
        @warn "Your data shows a cluster that starts before the first sample, or ends after the last sample. There exists a fundamental limit in the ClusterDepth method, that the clusterdepth for such a cluster cannot be determined. Maybe you can extend the epoch to include more samples? This 'half'-cluster will be ignored and cannot become significant"
    end
    cdmTuple = perm_clusterdepths_both(
        rng,
        data,
        permfun,
        τ;
        nₚ=nperm,
        (statfun!)=statfun!,
        statfun=statfun,
        sidefun=sidefun,
    )

    return pvals(data_obs, cdmTuple, τ; type=pval_type)
end




function perm_clusterdepths_both(
    rng,
    data,
    permfun,
    τ;
    statfun=nothing,
    (statfun!)=nothing,
    nₚ=1000,
    sidefun=nothing,
)
    @assert !(isnothing(statfun) && isnothing(statfun!)) "either statfun or statfun! has to be defined"

    data_perm = deepcopy(data)
    rows_h = Int[]
    cols_h = Int[]
    vals_h = Float64[]
    rows_t = Int[]
    cols_t = Int[]
    vals_t = Float64[]

    if ndims(data_perm) == 2
        d0 = Array{Float64}(undef, size(data_perm, 1))
    else
        d0 = Array{Float64}(undef, size(data_perm)[[1, 2]])
    end
    #@debug size(d0)
    #@debug size(data_perm)
    #data_perm .-= mean(data_perm, dims=length(size(data_perm)))
    for i = 1:nₚ
        # permute	
        d_perm = permfun(rng, data_perm)

        if isnothing(statfun!)
            d0 = statfun(d_perm)
        else
            # inplace!
            statfun!(d0, d_perm)
        end

        d0 .= sidefun.(d0)

        # get clusterdepth
        (fromTo, head, tail) = calc_clusterdepth(d0, τ)

        # save it
        if !isempty(head)
            append!(rows_h, fromTo)
            append!(cols_h, fill(i, length(fromTo)))
            append!(vals_h, head)

            #Jₖ_head[fromTo,i] .+=head
        end
        if !isempty(tail)
            append!(rows_t, fromTo)
            append!(cols_t, fill(i, length(fromTo)))
            append!(vals_t, tail)
            #Jₖ_tail[fromTo,i] .+=tail
        end
    end

    Jₖ_head = sparse(rows_h, cols_h, vals_h, maximum(rows_h), nₚ)#SparseMatrixCSC(nₚ,maximum(rows_h), cols_h,rows_h,vals_h)
    Jₖ_tail = sparse(rows_t, cols_t, vals_t, maximum(rows_t), nₚ)#SparseMatrixCSC(nₚ,maximum(rows_t), cols_t,rows_t,vals_t)
    return ClusterDepthMatrix((Jₖ_head)), ClusterDepthMatrix((Jₖ_tail))
end


"""
calc_clusterdepth(data,τ)

returns tuple with three entries:
1:maxLength, maximal clustervalue per clusterdepth head, same for tail

We assume data and τ have already been transformed for one/two sided testing, so that we can do d0.>τ for finding clusters

"""
function calc_clusterdepth(d0::AbstractArray{<:Real,2}, τ)
    nchan = size(d0, 1)

    # save all the results from calling calc_clusterdepth on individual channels
    (allFromTo, allHead, allTail) = (
        Array{Vector{Integer}}(undef, nchan),
        Array{Vector{Float64}}(undef, nchan),
        Array{Vector{Float64}}(undef, nchan),
    )
    fromTo = []
    for i = 1:nchan
        (a, b, c) = calc_clusterdepth(d0[i, :], τ)
        allFromTo[i] = a
        allHead[i] = b
        allTail[i] = c

        if (length(a) > length(fromTo)) # running check to find the length ('fromTo') of the largest cluster
            fromTo = a
        end
    end

    # for each clusterdepth value, select the largest cluster value found across all channels
    (head, tail) = (zeros(length(fromTo)), zeros(length(fromTo)))
    for i = 1:nchan
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

function calc_clusterdepth(d0, τ)
    startIX, len = cluster(d0 .> τ)
    if isempty(len) # if nothing above threshold, just go on
        return [], [], []
    end

    maxL = 1 + maximum(len) # go up only to max-depth
    #    @debug maxL
    valCol_head = Vector{Float64}(undef, maxL)
    valCol_tail = Vector{Float64}(undef, maxL)

    fromTo = 1:maxL
    for j in fromTo
        # go over clusters implicitly
        # select clusters that are larger (at least one)
        selIX = len .>= (j - 1)



        if !isempty(selIX)
            ix = startIX[selIX] .+ (j - 1)
            valCol_head[j] = maximum(d0[ix])

            # potential optimization is that for j = 0 and maxL = 0, tail and head are identical	
            ix = startIX[selIX] .+ (len[selIX]) .- (j - 1)
            valCol_tail[j] = maximum(d0[ix])
        end
    end
    return fromTo, valCol_head, valCol_tail
end

"""
finds neighbouring clusters in the vector and returns start + length vectors

if the first and last cluster start on the first/last sample, we dont know their real depth

	Input is assumed to be a thresholded Array with only 0/1
"""
function cluster(data::BitVector)


    start = []
    stop = []

    state = false
    for i in eachindex(data)
        if data[i] == 1
            if !state
                append!(start, i)
                state = true
            end
        else
            if state
                append!(stop, i - 1)
                state = false
            end
        end

    end





    # if the first and last cluster start on the first/last sample, we dont know their real depth
    if length(start) == length(stop) + 1
        start = start[1:end-1]
    end
    len = stop .- start

    if length(start) > 0 && start[1] == 1
        start = start[2:end]
        len = len[2:end]
    end
    return start, len
end
