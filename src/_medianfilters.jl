using Statistics

export thresholdmedian, thresholdmedian_threads, medianfilter

function thresholdmedian(img::Matrix{T}, meanlevel::Number, ksize::Int = 3) where T<:Number
    h, w = size(img)
    nn = (ksize-1)>>1
    result = copy(img)
    ml = convert(T, meanlevel)
    @inbounds for j in 1+nn:1:w-nn
        @inbounds for i in 1+nn:1:h-nn
            neighbors = img[i-nn:i+nn, j-nn:j+nn]
            nd = (sum(neighbors)-img[i, j])/(length(neighbors)-1)
            if ml > zero(T)
                if (abs(img[i, j] - nd) > ml*sqrt(nd) )
                    result[i, j] = median(neighbors)
                end
            end
        end
    end
    return result
end

function thresholdmedian_threads(img::Matrix{T}, meanlevel::Number, ksize::Int=3) where T<:Number
    h, w = size(img)
    nn = (ksize-1)>>1
    result = copy(img)
    ml = convert(T, meanlevel)
    @Threads.threads for j in 1+nn:1:w-nn
        for i in 1+nn:1:h-nn
            # @views 를 앞에 쓰면 memory allocation 은 작아지지만 속도가 느려진다.
            neighbors = img[i-nn:i+nn, j-nn:j+nn]
            nd = (sum(neighbors)-img[i, j])/(length(neighbors)-1)
            if meanlevel > zero(T)
                if (abs(img[i, j] - nd) > ml*sqrt(nd) )
                    result[i, j] = median(neighbors)
                end
            end
        end
    end
    return result
end



function medianfilter(img::Matrix{T}, ksize = 3) where T<:Number
    h, w = size(img)
    nn = (ksize-1)>>1
    
    result = zero(img)
    for j = nn+1:w-nn, i = nn+1:h-nn
        @views result[i, j] = median(img[i-nn:i+nn, j-nn:1:j+nn])
    end
       
    return result
end
