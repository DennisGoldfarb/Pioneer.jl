abstract type SpectralScores{T<:AbstractFloat} end

struct SpectralScoresComplex{T<:AbstractFloat} <: SpectralScores{T}
    spectral_contrast::T
    fitted_spectral_contrast::T
    gof::T
    max_matched_residual::T
    max_unmatched_residual::T
    fitted_manhattan_distance::T
    matched_ratio::T
    #entropy_score::T
end

struct SpectralScoresSimple{T<:AbstractFloat} <: SpectralScores{T} 
    scribe::T
    city_block::T
    spectral_contrast::T
    matched_ratio::T
    entropy_score::T
    num_scored_peaks::T
end

function getDistanceMetrics(H::SparseArray{Ti,T}, spectral_scores::Vector{SpectralScoresSimple{U}}) where {Ti<:Integer,T,U<:AbstractFloat}

    for col in range(1, H.n)
        h_sqrt_sum = zero(T)
        x_sqrt_sum = zero(T)

        h2_norm = zero(T)
        x2_norm = zero(T)

        dot_product = zero(T)
        scribe_score = zero(T)
        city_block_dist = zero(T)

        matched_sum = zero(T)
        unmatched_sum = zero(T)

        N = 0
        #@turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        @inbounds @fastmath for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if iszero(H.isotope[i])==false
                continue
            end

            h_sqrt_sum += sqrt(H.nzval[i])
            x_sqrt_sum += sqrt(H.x[i])#/Xsum

            h2_norm += H.nzval[i]^2
            x2_norm += H.x[i]^2

            dot_product += H.nzval[i]*H.x[i]

            if H.matched[i]
                matched_sum += H.nzval[i]
            else
                unmatched_sum += H.nzval[i]
            end

            N += 1
        end
          
        #Sqrt of sum of squares
        h2_norm = sqrt(h2_norm)
        x2_norm = sqrt(x2_norm)

        @inbounds @fastmath for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if iszero(H.isotope[i])==false
                continue
            end
            scribe_score += ((sqrt(H.nzval[i])/h_sqrt_sum) - (sqrt(H.x[i])/x_sqrt_sum))^2

            city_block_dist += abs(
                (H.nzval[i]/h2_norm) - (H.x[i]/h2_norm)
            )   
        end
        spectral_scores[col] = SpectralScoresSimple(
            Float16(-log((scribe_score)/N)), #scribe_score
            Float16(-log((city_block_dist/N))),
            #Float16(dotp_sqrt/(h_norm_m0*x_norm_m0)),#Float16(-log((city_block_dist)/N_M0)), #city_block
            Float16(dot_product/(h2_norm*x2_norm)), #spectral angle
            Float16(log2(matched_sum/unmatched_sum)), #matched_ratio
            Float16(Float16(-1.0)*getEntropy(H, col)), #entropy
            Float16(0.0f0)
        )
        #println(spectral_scores[col])
    end
end

function getDistanceMetrics(w::Vector{T}, r::Vector{T}, H::SparseArray{Ti,T}, spectral_scores::Vector{SpectralScoresComplex{U}}) where {Ti<:Integer,T,U<:AbstractFloat}


    ########
    #Get residuals
    #########
    @turbo for i in range(1, H.m)
        r[i] = zero(T)
    end

    for n in range(1, H.n_vals)
        if iszero(r[H.rowval[n]])
            r[H.rowval[n]] = -H.x[n]
        end
    end

    for col in range(1, H.n)
        start = H.colptr[col]
        stop = H.colptr[col+1] - 1
        for n in start:stop
            r[H.rowval[n]] += w[col]*H.nzval[n]
        end
    end

    for col in range(1, H.n)
        
        h2_sum = zero(T)
        x2_sum = zero(T)
        x_sum = zero(T)
        dot_product = zero(T)
        fitted_dotp = zero(T)
        matched_sum = zero(T)
        unmatched_sum = zero(T)    

        manhattan_distance = zero(T)
        max_matched_residual = zero(T)
        max_unmatched_residual = zero(T)
        sum_of_residuals = zero(T)
        fitted_dotp = zero(T)        
        fitted_dotp_norm1 = zero(T)
        sum_of_fitted_peaks_matched = zero(T)
        sum_of_fitted_peaks_unmatched = zero(T)

        @inbounds @fastmath for i in range(H.colptr[col], H.colptr[col + 1]-1)

            #Fitted Manhattan Distance
            x_sum += H.x[i]
            manhattan_distance += abs(w[col]*H.nzval[i] - H.x[i])

            #Normalized Dot Product 
            dot_product += H.nzval[i]*H.x[i]
            x2_sum += (H.x[i])^2
            h2_sum += (H.nzval[i])^2 
    

            sum_of_residuals += abs(r[H.rowval[i]])

            fitted_peak = w[col]*H.nzval[i]
            r_abs = abs(r[H.rowval[i]])
            sum_of_residuals += r_abs

            if H.matched[i]
                matched_sum += H.nzval[i]
                fitted_dotp += sqrt(fitted_peak + r_abs)*sqrt(fitted_peak)
                fitted_dotp_norm1 += fitted_peak + r_abs
                sum_of_fitted_peaks_matched += fitted_peak
                if r_abs > max_matched_residual
                    max_matched_residual = r_abs
                end
            else
                unmatched_sum += H.nzval[i]
                sum_of_fitted_peaks_unmatched += fitted_peak
                if r_abs > max_unmatched_residual
                    max_unmatched_residual = r_abs
                end
            end
        end
        sum_of_fitted_peaks =  sum_of_fitted_peaks_matched +  sum_of_fitted_peaks_unmatched

        fitted_dotp_norm = fitted_dotp/(sqrt(fitted_dotp_norm1)*sqrt(sum_of_fitted_peaks))
        fitted_spectral_contrast = fitted_dotp_norm#1 - 2*acos(fitted_dotp_norm)/π
        dot_product_norm = dot_product/(sqrt(h2_sum)*sqrt(x2_sum))
        spectral_contrast = dot_product_norm#1 - 2*acos(dot_product_norm)/π

        spectral_scores[col] = SpectralScoresComplex(
            Float16(spectral_contrast), #spectral_contrast
            Float16(fitted_spectral_contrast), #fitted_spectral_contrast
            Float16(-log2(sum_of_residuals/sum_of_fitted_peaks)), #gof
            Float16(-log2(max_matched_residual/sum_of_fitted_peaks_matched)),#max_matched_residual
            Float16(-log2(max_unmatched_residual/sum_of_fitted_peaks + 1e-10)), #max_unmatched_residual
            Float16(-log2(manhattan_distance/x_sum)), #fitted_manhattan_distance
            Float16(log2(matched_sum/unmatched_sum)), #matched_ratio
            #Float16(-1.0*getEntropy(H, r, col)) #entropy
        )
    end
end
    
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        #iso = (one(UInt8)-H.isotope[i])
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += xp*log(xp + Float32(1e-10))
        Hentropy += hp*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #iso = (one(UInt8)-H.isotope[i])
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += xp*log(xp + Float32(1e-10))
            Xsum += xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            #iso = (one(UInt8)-H.isotope[i])
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, r::Vector{T}, w::Vector{T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @inbounds @fastmath for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        if H.matched[i]
            hp = H.nzval[i]*w[col]
            xp = H.nzval[i]*w[col] + abs(r[H.rowval[i]])
        else
            hp = zero(T)#H.nzval[i]*w[col]
            xp =  abs(r[H.rowval[i]])
        end
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += (xp)*log(xp + Float32(1e-10))
        Hentropy += (hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    #=
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @inbounds @fastmath  for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if H.matched[i]
                xp = ((H.nzval[i]*w[col] + abs(r[H.rowval[i]]))^Xw)#/Xsum
                Xentropy += xp*log(xp + Float32(1e-10))
                Xsum += xp
            end
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @inbounds @fastmath  for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise    
            if H.matched[i]
                hp = H.nzval[i]*w[col]
                xp = H.nzval[i]*w[col] + abs(r[H.rowval[i]])
            else
                hp = H.nzval[i]*w[col]
                xp = zero(T)
            end
            hxp = (hp + xp)^HXw
            hp = (hp^Hw) #/Hsum
            xp = (xp^Xw)#/Xsum
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    =#
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, r::Vector{T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += (xp)*log(xp + Float32(1e-10))
        Hentropy += (hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += xp*log(xp + Float32(1e-10))
            Xsum += xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise    
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
#=

function getEntropy(H::SparseArray{Ti, T}, r::Vector{T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        iso = one(UInt8)-H.isotope[i]
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += iso*xp
        Hsum += iso*hp
        HXentropy += iso*(hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += iso*(xp)*log(xp + Float32(1e-10))
        Hentropy += iso*(hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            iso = (one(UInt8)-H.isotope[i])
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            Xsum += iso*xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            iso = (one(UInt8)-H.isotope[i])
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += iso*hp*log(hp + Float32(1e-10))
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            HXentropy += iso*hxp*log(hxp + Float32(1e-10))
            Xsum += iso*xp
            Hsum += iso*hp
            HXsum += iso*hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end

function getEntropy2(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        iso = one(UInt8)#H.isotope[i]
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += iso*xp
        Hsum += iso*hp
        HXentropy += iso*(hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += iso*(xp)*log(xp + Float32(1e-10))
        Hentropy += iso*(hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            iso = one(UInt8)#H.isotope[i]
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            Xsum += iso*xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            iso = one(UInt8)#H.isotope[i]
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += iso*hp*log(hp + Float32(1e-10))
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            HXentropy += iso*hxp*log(hxp + Float32(1e-10))
            Xsum += iso*xp
            Hsum += iso*hp
            HXsum += iso*hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
=#
#=
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end

function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
=#
function getEntropy(H::SparseArray{Ti, T}) where {Ti<:Integer,T<:AbstractFloat}
    entropy_sim = zeros(T, H.n)
    for col in range(1, H.n)

        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]#/Hsum
            xp = H.x[i]#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            #HXentropy += (hp + xp)*log(hp + xp)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
            #HXsum += xp + hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]/Hsum
            xp = H.x[i]/Xsum
            HXentropy += (hp + xp)*log(hp + xp + 1e-10)
            HXsum += xp + hp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        
        if (Xentropy < 3) | (Hentropy < 3)

            Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
            Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
            HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

            Hentropy = zero(T)
            Xentropy = zero(T)
            HXentropy = zero(T)
            Hsum = zero(T)
            Xsum = zero(T)
            HXsum = zero(T)
            @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
                #MASK is true for selected ions and false otherwise
                hp = H.nzval[i]^Hw#/Hsum
                xp = H.x[i]^Xw#/Xsum
                Hentropy += hp*log(hp + 1e-10)
                Xentropy += xp*log(xp + 1e-10)
                Xsum += xp
                Hsum += hp
            end

            @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
                #MASK is true for selected ions and false otherwise
                hp = (H.nzval[i]^Hw)/Hsum
                xp = (H.x[i]^Xw)/Xsum
                hxp = (hp + xp)^HXw
                HXentropy += (hxp)*log(hxp + 1e-10)
                HXsum += hxp
            end

            Xentropy = log(Xsum) - Xentropy/Xsum
            HXentropy = log(HXsum) - HXentropy/HXsum 
            Hentropy = log(Hsum) - Hentropy/Hsum
            HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
            HXentropy = zero(Float32)
            HXsum = zero(Float32)
            @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
                #MASK is true for selected ions and false otherwise
                hp = (H.nzval[i]^Hw)/Hsum
                xp = (H.x[i]^Xw)/Xsum
                hxp = (hp + xp)^HXw
                HXentropy += (hxp)*log(hxp + 1e-10)
                HXsum += hxp
            end

            HXentropy = log(HXsum) - HXentropy/HXsum
        end

        entropy_sim[col] = (1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
    end
    return entropy_sim
end


#= function getDistanceMetrics(
    H::SparseArray{Ti,T},
    spectral_scores::Vector{SpectralScoresSimple{U}},
    spectral_angle_threshold::Float32,
    min_peaks::Int64
) where {Ti<:Integer, T, U<:AbstractFloat}

    # We'll iterate over each column (mass spectrum)
    for col in 1:H.n
        # Gather the indices of relevant peaks in this spectrum
        included = Int[]  # We'll store the indices from H.colptr[col] : H.colptr[col+1]-1
        for i in H.colptr[col]:(H.colptr[col+1]-1)
            # Skip if it's an isotope or something we don't consider
            if !iszero(H.isotope[i])
                continue
            end
            push!(included, i)
        end

        # Repeatedly compute metrics, check spectral angle, possibly remove worst peak
        done = false
        while !done
            # 1) Compute metrics on the current set of included peaks
            scribe_val, city_val, spec_angle, matched_ratio, ent, num_scored_peaks = computeMetricsFor(H, included)

            # 2) Check the spectral angle threshold
            if spec_angle < spectral_angle_threshold && length(included) > min_peaks
                # Remove the "worst matching" peak
                worst_idx = findWorstMatchingIndex(H, included)
                deleteat!(included, worst_idx)
            else
                # Either the angle is good enough or we can't remove more
                # => store final metrics and exit loop
                #println(scribe_val, " ", city_val, " ", spec_angle, " ", matched_ratio, " ", ent)
                spectral_scores[col] = SpectralScoresSimple(
                    Float16(scribe_val),
                    Float16(city_val),
                    Float16(spec_angle),
                    Float16(matched_ratio),
                    Float16(ent),
                    UInt8(num_scored_peaks)
                )
                done = true
            end
        end
    end
end =#


function getDistanceMetrics(
    H::SparseArray{Ti,T},
    spectral_scores::Vector{SpectralScoresSimple{U}},
    cosine_threshold::Float32,
    max_drop::Int64
) where {Ti<:Integer, T, U<:AbstractFloat}

    # We'll iterate over each column (mass spectrum)
    for col in 1:H.n
        # Gather the indices of relevant peaks in this spectrum
        included = Int[]  # We'll store the indices from H.colptr[col] : H.colptr[col+1]-1
        tot_pred_signal = 0.0
        for i in H.colptr[col]:(H.colptr[col+1]-1)
            # Skip if it's an isotope or something we don't consider
            if !iszero(H.isotope[i])
                continue
            end
            push!(included, i)
            tot_pred_signal += H.nzval[i]
        end


        done = false
        dropped = 0

        scribe_val_current, city_val_current, cosine_current, matched_ratio_current, ent_current, num_scored_peaks_current = computeMetricsFor(H, col, included)
        #println("dropped: ", dropped, " scribe: ", scribe_val_current, " city: ", city_val_current, " cosine: ", spec_angle_current, " matched_ratio: ", matched_ratio_current, " ent: ", ent_current, " num_scored_peaks: ", num_scored_peaks_current)
        num_scored_peaks_current = 0
        while !done
            worst_idx, pred_signal_lost = findWorstMatchingIndex(H, included)
            deleteat!(included, worst_idx)
            scribe_val, city_val, cosine, matched_ratio, ent, num_scored_peaks = computeMetricsFor(H, col, included)
            dropped +=1

            #println("dropped: ", dropped, " scribe: ", scribe_val, " city: ", city_val, " cosine: ", spec_angle, " matched_ratio: ", matched_ratio, " ent: ", ent, " num_scored_peaks: ", num_scored_peaks)

            if (scribe_val_current > (scribe_val - 1)) || isinf(scribe_val) || isnan(scribe_val) || (length(included) < 3) || (dropped > max_drop)#|| (dropped > 3) # || (length(included) < 20)
                break
            end

            scribe_val_current = scribe_val
            city_val_current = city_val#
            cosine_current = cosine#
            matched_ratio_current = matched_ratio#
            ent_current = ent#
            num_scored_peaks_current += pred_signal_lost / tot_pred_signal
        end

        spectral_scores[col] = SpectralScoresSimple(
                    Float16(scribe_val_current),
                    Float16(city_val_current),
                    Float16(cosine_current),
                    Float16(matched_ratio_current),
                    Float16(ent_current),
                    Float16(num_scored_peaks_current)
                )

        if dropped > 4
            println("scribe: ", scribe_val_current, " city: ", city_val_current, " cosine: ", cosine_current, " matched_ratio: ", matched_ratio_current, " ent: ", ent_current, " num_scored_peaks: ", num_scored_peaks_current, " dropped: ", dropped-1)
        end
        #println(num_scored_peaks_current)

        #exit()


        # Repeatedly compute metrics, check spectral angle, possibly remove worst peak
        ##done = false
        ##dropped = 0
        ##while !done
            # 1) Compute metrics on the current set of included peaks
        ##   scribe_val, city_val, spec_angle, matched_ratio, ent, num_scored_peaks = computeMetricsFor(H, col, included)

        ##   worst_idx = findWorstMatchingIndex(H, included)
        ##    deleteat!(included, worst_idx)
            
        ##    scribe_val2, city_val2, spec_angle2, matched_ratio2, ent2, num_scored_peaks2 = computeMetricsFor(H, col, included)

        ##    spectral_scores[col] = SpectralScoresSimple(
        ##            Float16(max(scribe_val, scribe_val2)),
        ##            Float16(min(city_val, city_val2)),
        ##            Float16(max(spec_angle, spec_angle2)),
        ##            Float16(max(matched_ratio, matched_ratio2)),
        ##            Float16(min(ent, ent2)),
        ##            UInt8(num_scored_peaks)
        ##        )
        ##    done = true
            
            
           #=  # 2) Check the spectral angle threshold
            if spec_angle < spectral_angle_threshold && length(included) > 3 && dropped < max_drop
                # Remove the "worst matching" peak
                worst_idx = findWorstMatchingIndex(H, included)
                deleteat!(included, worst_idx)
                dropped += 1
            else
                # Either the angle is good enough or we can't remove more
                # => store final metrics and exit loop
                spectral_scores[col] = SpectralScoresSimple(
                    Float16(scribe_val),
                    Float16(city_val),
                    Float16(spec_angle),
                    Float16(matched_ratio),
                    Float16(ent),
                    UInt8(num_scored_peaks)
                )
                done = true
            end =#
        #end
    end
end

function computeMetricsFor(H, col, included_indices)
    # We'll accumulate partial sums and compute the same metrics as your snippet.
    # (For clarity, we skip inlining optimizations like @inbounds, @fastmath here.)

    N = length(included_indices)
    if N == 0
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0)
    end

    # Sums for numerator/denominator
    # We also want the "normalized" predicted and observed, so we sum intensities
    total_h = 0.0
    total_x = 0.0

    for i in included_indices
        total_h += H.nzval[i]
        total_x += H.x[i]
    end

    # Protect against zero total intensities
    total_h = max(total_h, eps()) 
    total_x = max(total_x, eps())

    # Sums for scribe, city-block, dot-product, norms, etc.
    h_sqrt_sum = 0.0
    x_sqrt_sum = 0.0

    h2_norm = 0.0
    x2_norm = 0.0

    h2_norm_v2 = 0.0
    x2_norm_v2 = 0.0

    dot_product = 0.0
    matched_sum = 0.0
    unmatched_sum = 0.0

    for i in included_indices

        # "Normalized" predicted, observed
        h_val_v2 = H.nzval[i] / total_h
        x_val_v2 = H.x[i] / total_x
        dot_product += h_val_v2 * x_val_v2
        
        h2_norm_v2 += h_val_v2^2
        x2_norm_v2 += x_val_v2^2

        # Regular predicted, observed
        h_val = H.nzval[i]
        x_val = H.x[i]

        h_sqrt_sum += sqrt(h_val)
        x_sqrt_sum += sqrt(x_val)

        h2_norm += h_val^2
        x2_norm += x_val^2

        if H.matched[i]
            matched_sum += h_val
        else
            unmatched_sum += h_val
        end
    end

    # sqrt of sums
    h2_norm = sqrt(h2_norm)
    x2_norm = sqrt(x2_norm)

    h2_norm_v2 = sqrt(h2_norm_v2)
    x2_norm_v2 = sqrt(x2_norm_v2)

    # Scribe and city-block 
    # (using the approach from your snippet, but now we've pre-normalized intensities)
    scribe_val = 0.0
    city_val   = 0.0
    for i in included_indices
        h_val = H.nzval[i]
        x_val = H.x[i]

        scribe_val += ( (sqrt(h_val)/h_sqrt_sum) - (sqrt(x_val)/x_sqrt_sum) )^2
        city_val   += abs( (h_val/h2_norm) - (x_val/h2_norm) )
    end
    scribe_val = -log2(scribe_val / N)
    city_val   = -log2(city_val / N)

    spectral_angle = dot_product / (h2_norm_v2 * x2_norm_v2)

    # matched_ratio = log2(matched_sum/unmatched_sum)
    # protect from zero unmatched
    unmatched_sum = max(unmatched_sum, eps())
    matched_ratio = log2(matched_sum / unmatched_sum)

    # Entropy – if your original code uses `getEntropy(H, col)` you might
    # want a version that only considers `included_indices`. We'll do a simple
    # placeholder for now:
    #ent_val = -1.0 * getEntropy(H, col, included_indices)  # you’d define getEntropy to handle included
    ent_val = Float16(-1.0) * getEntropy(H, col) 

    return (scribe_val, city_val, spectral_angle, matched_ratio, ent_val, length(included_indices))
end


function findWorstMatchingIndex(H, included_indices)
    # We'll pick the i that yields the max difference between (nzval[i], x[i])
    # Return the *position* within included_indices (the list index to delete).
    # That way we can do `deleteat!(included, worst_idx)`.
    worst_val = -Inf
    worst_pos = 0
    worst_idx = 0

    # First, we could compute total predicted, total observed to normalize
    total_h = 0.0
    total_x = 0.0
    for idx in included_indices
        total_h += H.nzval[idx]
        total_x += H.x[idx]
    end
    total_h = max(total_h, eps())
    total_x = max(total_x, eps())

    for (local_pos, idx) in enumerate(included_indices)
        h_val = H.nzval[idx] / total_h
        x_val = H.x[idx] / total_x
        diff  = abs(h_val - x_val)
        if diff > worst_val
            worst_val = diff
            worst_pos = local_pos
            worst_idx = idx
        end
    end

    return (worst_pos, H.nzval[worst_idx])
end

function getEntropy(H, included_indices)
    # e.g. Shannon entropy over normalized intensities in included
    total = 0.0
    for idx in included_indices
        total += H.nzval[idx]
    end
    total = max(total, eps())
    entropy = 0.0
    for idx in included_indices
        p = H.nzval[idx] / total
        entropy -= p*log2(p + eps())
    end
    return entropy
end

function getEntropy(H::SparseArray{Ti, T}, col::Int64, included_indices) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    for i in range(H.colptr[col], H.colptr[col + 1]-1)
        if i ∉ included_indices
            continue
        end
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        #iso = (one(UInt8)-H.isotope[i])
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += xp*log(xp + Float32(1e-10))
        Hentropy += hp*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if i ∉ included_indices
                continue
            end
            #iso = (one(UInt8)-H.isotope[i])
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += xp*log(xp + Float32(1e-10))
            Xsum += xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if i ∉ included_indices
                continue
            end
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            #iso = (one(UInt8)-H.isotope[i])
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end