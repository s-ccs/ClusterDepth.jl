
"""
pvals(data;kwargs...) = pvals(data[2:end],data[1];kwargs...)

pvals(data::AbstractVector,stat::Real;type=:twosided)

calculates pvalues based on permutation results

if called with `stat`, first entry is assumed to be the observation 
"""
pvals(data;kwargs...) = pvals(data[2:end],data[1];kwargs...)
function pvals(data::AbstractVector,stat::Real;type=:twosided)

    data = vcat(stat,data)
    if type == :greater || type  == :twosided
        comp = >=
        if type == :twosided
            data = abs.(data)
            stat = abs.(stat)
        end
    elseif type == :lesser
        comp = <=
    else
        error("not implemented")
    end

    pvals = (sum(comp(stat[1]),data))/(length(data))
    return pvals
end



"""
Calculate pvals from cluster-depth permutation matrices
"""
pvals(stat,J::ClusterDepthMatrix,args...;kwargs...) = pvals(stat,(J,),args...;kwargs...)
function pvals(stat,Jₖ::NTuple{T,ClusterDepthMatrix},τ;type=:troendle) where{T}

    start,len = cluster(stat,τ) # get observed clusters
    p = fill(1.,size(stat,1))

    if type == :troendle
               
        for k = 1:length(start) # go over clusters
            s = start[k]
            l = len[k]
            forwardIX = s:(s+l)
            

            @views t_head = troendle(Jₖ[1],sparsevec(1:(l+1), stat[forwardIX],l+1))
            @views p[forwardIX] = t_head[1:(l+1)]
            if length(Jₖ) == 2
                backwardIX = (s+l):-1:s
                @views t_tail = troendle(Jₖ[2],sparsevec(1:(l+1),stat[backwardIX],l+1))
                @views p[backwardIX] = max.(p[backwardIX],t_tail[1:(l+1)])
            end
        end
    elseif type == :naive
        function getJVal(Jₖ,l)
            if l>= size(Jₖ,1)
                valsNull = 0
            else
                valsNull = @view Jₖ[l+1,:]
            end 
            return valsNull
        end
        
        for k = 1:length(start) # go over clusters
            for ix = start[k]: (start[k]+len[k]) 
             
                    p[ix] = (1+sum(stat[ix] .<= getJVal(Jₖ[1].J,len[k])))/(size(Jₖ[1].J,2)+1)
                if length(Jₖ) == 2
                    tail_p = (1+sum(p[ix],mean(stat[ix] .<= getJVal(Jₖ[2].J,len[k])))) ./(size(Jₖ[2].J,2)+1)
                    p[ix] = max(p[ix],tail_p)
                end
            end
        end
        #p = p .+ 1/size(Jₖ[1].J,2)
    else 
        error("unknown type")
    end

    # add conservative permutation fix
    
    return p
    
end