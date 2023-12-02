using Images, ColorSchemes, CUDA

"""
convert matrix to gray image of Images.jl
"""
function mat2gray(mat)
    mv, Mv = extrema(mat)
    if typeof(mat) <: CuMatrix
        return Gray.((Array(mat).-mv)./(Mv-mv))
    else
        return Gray.((mat.-mv)./(Mv-mv))
    end
end

function colorize(img::Union{Matrix{T}, CuMatrix{T}}, color::Symbol, scale=:linear, minval::Union{Nothing, Real}=nothing, maxval::Union{Nothing, Real}=nothing) where T<:Real
    m, n = size(img)
    mv, Mv = extrema(img)

    if typeof(img) <: CuMatrix
        img = Array(img)
    end

    colorscheme = colorschemes[color]

    L = length(colorscheme)

    if minval !== nothing && maxval !== nothing
        mv, Mv = minval, maxval
    end
    if minval !== nothing
        mv = minval
    end
    if maxval !== nothing
        Mv = maxval
    end    

    @assert mv < Mv

    dd = (Mv-mv)/L

    result = zeros(RGB, m, n)

    @Threads.threads for I in eachindex(img)
        ll = round(Int64, (img[I]-mv)/dd)

        if ll < 1
            @inbounds result[I] = colorscheme[ll+1]
        elseif ll > (L-1)
            @inbounds result[I] = colorscheme[end]
        else 
            @inbounds result[I] = colorscheme[ll]
        end
    end
    return result
end

function colorize(img::Matrix{Gray}, colorscheme::ColorScheme, scale=:linear, minval::Union{Nothing, Real}=nothing, maxval::Union{Nothing, Real}=nothing)
    return colorize(Float64.(img), colorscheme, scale, minval, maxval)
end


function _thist(img::Matrix{T}, dedge::S) where {T<:Real, S<:Real}
    @assert dedge > zero(dedge)
    
    M, N = size(img)
    etp = eltype(img)
    
    mm, MM = extrema(img)
    edges = mm:dedge:(MM+dedge)
    counts = zeros(Int, length(edges)-1)
    for j in 1:N, i in 1:M
        counts[floor(Int64, (img[i, j]-mm)/dedge)+1]+=1
    end
    return edges, counts
end

