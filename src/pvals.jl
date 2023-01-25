"""
Calculate pvals from cluster-depth permutation matrices
"""

function pvals(stat,Jₖ;type=:troendle)

    start,len = get_clusterdepths(stat,τ) # get observed clusters
    p = fill(1.,size(stat,1))
    if type == :troendle
        if isa(Jₖ,Tuple)
            twoside = true
            Jₖ_head = Jₖ[1]
            Jₖ_tail = Jₖ[2]
        else 
            twoside = false
            Jₖ_head = Jₖ
        end
        
        for k = 1:length(start) # go over clusters
            s = start[k]
            l = len[k]
            forwardIX = s:(s+l)
                @views t_head = troendle(Jₖ_head,sparsevec(1:(l+1), stat[forwardIX],l+1))
            @views p[forwardIX] = t_head[1:(l+1)]
            if twoside
                backwardIX = (s+l):-1:s
                @views t_tail = troendle(Jₖ_tail,sparsevec(1:(l+1),stat[backwardIX],l+1))
                @views p[backwardIX] = max.(p[backwardIX],t_tail[1:(l+1)])
            end
        end
    elseif type == :naive
        for k = 1:length(start) # go over clusters
            for ix = start[k]: (start[k]+len[k]) 
                if len[k]>= size(Jₖ,1)
                    valsNull = 0
                else
                    valsNull = @view Jₖ[len[k]+1,:]
                end
                p[ix] = mean(stat[ix] .<= valsNull)
            end
        end
    else 
        error("unknown type")
    end
    p = p .+ 1/nₚ
    return p
    
end