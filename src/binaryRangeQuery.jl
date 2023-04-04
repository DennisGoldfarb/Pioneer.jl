function binaryGetNearest(arr::Vector{Union{Missing, T}}, query::T, low_tol::T, high_tol::T) where T<:Real

    #Check special cases (is the answer on the boundary or is the array empty?)
    n = length(arr)
    if n == 0 return 0 end
    if query < arr[1] - low_tol return Missing end
    if query > arr[n] + high_tol return Missing end

    function getNearest(arr::Vector{Union{Missing, T}}, lo::Int, hi::Int, query::T) where T<:Real
        if hi - lo>1
            smallest_distance = abs(query - arr[lo])
            best_idx = 1
            for (i, mass) in enumerate(@view(arr[lo:hi]))
                if abs(query - mass) < smallest_distance
                    smallest_distance = abs(query - mass)
                    best_idx = i
                end
            end
            return best_idx
        else
            return 1
        end
    end
    lo, hi = 1, n
    while lo <= hi
        mid = (lo + hi) ÷ 2
        if arr[mid] < (query - low_tol)
             lo = mid + 1
        elseif arr[mid] > (query + high_tol)
            hi = mid - 1
        else
            return getNearest(arr,lo, hi, query) + lo - 1
        end
    end

    return Missing

end

function binaryRangeQuery()
function getPrecursors(window_center::Float32, precursorList::Vector{Precursor}, params)
    l_bnd, u_bnd = window_center - params[:lower_tol], window_center + params[:upper_tol]
    start, stop = searchsortedfirst(precursorList, l_bnd,lt=(t,x)->getMZ(t)<x), searchsortedlast(precursorList, u_bnd,lt=(x,t)->getMZ(t)>x)
    return @view(precursorList[start:stop])
end


#searchsortedlast(bestPSMs[!,:retentionTime],  x, lt=(t, x)->t<x)
searchsortedfirst(bestPSMs[!,:retentionTime], 50 ,lt=(t,x)->t<x)
searchsortedlast(bestPSMs[!,:retentionTime], 50 ,lt=(x, t)->t>x)
searchsortedfirst(bestPSMs[!,:retentionTime],  x, ltx->50>x)

bestPSMs[170:180,:]