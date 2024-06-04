using Images, ColorSchemes, CairoMakie

"""
    mat2gray(mat::Matrix{<:Real})

Convert matrix to gray image of Images.jl
"""
function mat2gray(mat::Matrix{<:Real})
    mv, Mv = extrema(mat)
    
    return Gray.((mat.-mv)./(Mv-mv))
    
end

"""
    colorize(img, color, scale, minval, maxval)

Add color to 2D image. Color schemes uses ColorSchemes.jl[1].  

Reference
---------
[1] https://juliagraphics.github.io/ColorSchemes.jl/stable/
"""
function colorize(
    img::Matrix{<:Real}, 
    color::Symbol, 
    scale=:linear, 
    minval::Union{Nothing, Real}=nothing, 
    maxval::Union{Nothing, Real}=nothing)

    m, n = size(img)
    mv, Mv = extrema(img)

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

function colorize(img::Matrix{Gray}, 
    colorscheme::ColorScheme, 
    scale=:linear, minval::Union{Nothing, Real}=nothing, 
    maxval::Union{Nothing, Real}=nothing)
    return colorize(Float64.(img), colorscheme, scale, minval, maxval)
end

"""
    _thist(img::Matrix{T}, dedge::S) where {T<:Real, S<:Real}

Get histogram of matrix. dedge is the inteval of bins
"""
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


"""
    image_profile(img, figsize, binsize)


Generate makie figure with two plots. The left one is image and the right one is histogram of it
"""
function image_profile(img::Matrix; figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing, binsize=2)
    
    if isa(figsize, Nothing) == false
        @assert (400 ≤ figsize[1] ≤ 10000) && (400 ≤  figsize[2] ≤ 6000)
    else         
        figsize = (800, 400)
    end

    f = Figure(size = figsize)

    Mv = maximum(img)
    v1 = image(f[1, 1:2], img', axis=(aspect = DataAspect(), yreversed=true, xticklabelsvisible=false, yticklabelsvisible=false))
    v2 = hist(f[1, 3], vec(img), bins = 0:binsize:Mv)

    return f
end

"""
vscode jupyter notebook 에서 각 함수, 객체 등에 관한 docstring 을 markdown 으로 출력한다.

## Example

@h sin
"""
macro h(x)
    quote
        display("text/markdown", @doc $x)
    end    
end
