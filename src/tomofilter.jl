using Statistics, Images
 
abstract type AbstractTomoFilter end

struct IdentityFilter <: AbstractTomoFilter
    ksize::Integer

    function IdentityFilter(ksize::Integer)
        @assert isodd(ksize) && ksize > 1
        return new(ksize)
    end
end

struct MedianFilter <: AbstractTomoFilter
    ksize::Integer
    use_threads::Bool 

    function MedianFilter(ksize::Integer=3, use_threads::Bool = false)
        @assert isodd(ksize) && ksize > 1
        return new(ksize, use_threads)
    end
end


struct ThreasholdMedianFilter <: AbstractTomoFilter
    ksize::Integer
    meanlevel::Real
    use_threads::Bool

    function ThreasholdMedianFilter(ksize::Integer, meanlevel::Real, use_threads::Bool=false)
        @assert isodd(ksize) && ksize > 1
        @assert meanlevel > 1
        return new(ksize, meanlevel, use_threads)
    end
        
end

function (p::IdentityFilter)(img::Matrix{<:Real})
    return img[:, :]    
end

function (p::MedianFilter)(img::Matrix{<:Real})
    result = mapwindow(median, img, (p.ksize, p.ksize))
    # h, w = size(img)
    # nn = (p.ksize-1)>>1
    
    # result = copy(img)
    # if p.use_threads 
    #     Threads.@threads for j = 1:w
    #         for i = 1:h
    #             @inbounds result[i, j] = median(img[max(1,i-nn):min(h, i+nn), max(1, j-nn):min(w, j+nn)])
    #         end
    #     end
    # else
    #     for j = 1:w, i = 1:h
    #         @inbounds result[i, j] = median(img[max(1,i-nn):min(h, i+nn), max(1, j-nn):min(w, j+nn)])
    #     end
    # end
    return result
end

function (p::ThreasholdMedianFilter)(img::Matrix{<:Real})
    h, w = size(img)
    nn = (p.ksize-1)>>1
    result = copy(img)
    Elt = eltype(img)
    ml = convert(Elt, p.meanlevel)
    if p.use_threads
        @Threads.threads for j in 1:w
            for i in 1:h
                # @views 를 앞에 쓰면 memory allocation 은 작아지지만 속도가 느려진다.
                neighbors = img[max(1, i-nn):min(h, i+nn), max(1, j-nn):min(w, j+nn)]
                nd = (sum(neighbors)-img[i, j])/(length(neighbors)-1)
                if (abs(img[i, j] - nd) > ml*sqrt(nd) )
                    result[i, j] = median(neighbors)
                end
            end
        end
    else 
        @inbounds for j in 1:w, i in 1:h
            neighbors = img[max(1, i-nn):min(h, i+nn), max(1, j-nn):min(w, j+nn)]
            nd = (sum(neighbors)-img[i, j])/(length(neighbors)-1)
                if (abs(img[i, j] - nd) > ml*sqrt(nd) )
                result[i, j] = median(neighbors)
            end
        end

    end
    return result
end
