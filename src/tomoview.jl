include("hanaro.jl")


function image_profile(tomo::TomoReader, kind = :data, index::Int64 = 1; figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing, binsize=2, xlog::Bool=false, ylog::Bool=false)

    img = get_data_by_index(tomo, kind, index)
    
    return image_profile(img, figsize=figsize, binsize=blisize, xlog=xlog, ylog=ylog)
end


function get_data_by_index(tomo::TomoReader, kind = :data, index::Integer=0)
    @assert kind ∈ (:data, :white, :dark)    
    if kind == :data
        @assert 1 ≤ index ≤ length(tomo.data_files)
        img = read_nrimage(joinpath(tomo.data_dir, (tomo.data_files[index])[3]), tomo.scale_down)
    elseif  kind == :dark
        @assert 1 ≤ index ≤ length(tomo.dark_files)
        img = read_nrimage(joinpath(tomo.dark_dir, (tomo.dark_files[index])), tomo.scale_down)
    elseif kind == :white 
        @assert 1 ≤ index ≤ length(tomo.white_files)
        img = read_nrimage(joinpath(tomo.white_dir, (tomo.white_files[index])), tomo.scale_down)
    end

    return img
end

function get_image_by_index(tomo::TomoReader, kind = :data, index::Integer=0)
    img = get_data_by_index(tomo, kind, index)
    return mat2gray(Float64.(img))
end

    

function get_image(tomo::TomoReader, angle::Real=0.0, contrast::Real = 1.0; 
    threashold::Union{Real, Nothing}=nothing)
    
    @assert 0 < contrast <= 1.0
    @assert 0.1 ≤ resize_scale ≤ 10.0

    contrast = Float32(contrast)
    
    if threashold ≠ nothing
        threashold = Float32(threashold)
    end
    
    ths = [th for (obj, th, fn) ∈ tomo.data_files if obj == 1]
    if length(ths) == 0 
        @error "No images for object=$object, angle=$angle"
        return 
    end
    
    m = argmin(abs.(angle .- ths))
    img = Float32.(read_nrimage(joinpath(tomo.data_dir, (tomo.data_files[m])[3])), tomo.scale_down)
    thv = 1.0f0
    if threashold === nothing
        img /= maximum(img)/contrast
    else 
        img /= threashold/contrast
        img[img.>1.0f0] .= 1.0f0
    end
    

    if tomo.to_be_transposed
        img = Gray.(img')
    else 
        img = Gray.(img)
    end

    if abs(resize_scale-1.0)>0.1 
        resized = round.(Int64, [r*resize_scale for r in size(img)])
        img = imresize(img, (resized[1], resized[2]))
    end
    return img
end 

function get_image(tomo::TomoReader, angles::AbstractVector, contrast::Real = 1.0) 
    imgs = [get_image(tomo, th, contrast) for th ∈ angles]
    return (imgs |> sum)
end


function get_image(tomo::TomoData, i::Integer , contrast=0.99)
    @assert 1 ≤ i ≤ size(tomo.data)[1]
    return mat2gray(tomo.data[i, :, :], contrast)
end

function get_histogram(tomo::TomoData, kind = :data, index = nothing, dbin=1.0)
    @assert kind ∈ (:data, :white, :dark)
    @assert index === nothing || typeof(index) <:Integer 
    if index === nothing 
        @assert kind ∈ (:white, :dark)
    end

    if kind == :data 
        p = tomo.data[index, :, :]
    elseif kind == :white && index ≠ nothing
        p = tomo.white[index, :, :]
    elseif kind == :white 
        p = sum(tomo.white, dims=1)[1, :, :]
    elseif kind == :dark && index ≠ nothing
        p = tomo.dark[index, :, :]
    else 
        p = sum(tomo.dark, dims=1)[1, :, :]
    end

    return _thist(p, dbin)
end


function get_sinogram(tomo::TomoData, i=1, color=:gray)
    @assert :sinogram ∈ tomo.process
    @assert color ∈ (:gray, :rgb)

    p = tomo.data[:, i, :]
    mp, Mp = extrema(p)
    r = (p .- mp)./(Mp - mp)
    if color == :gray
        return Gray.(r)
    elseif  color == :rgb
        return RGB.(r, r, r)
    end
end



function add_rect(img::Matrix{Gray{T}}, p1, p2, intensity = 1.0, lw=3) where T<:Real
    @assert 0.0 ≤ intensity ≤ 1.0
    lc = Gray(intensity)
    M, N = size(img)
    lw = Int(lw) >>1
    @assert 0<lw<10
    x1, x2 = minmax(Int(p1[1]), Int(p2[1]))
    y1, y2 = minmax(Int(p1[2]), Int(p2[2]))
    
    @assert 1 ≤ x1 ≤ N
    @assert 1 ≤ x2 ≤ N
    @assert 1 ≤ y1 ≤ M
    @assert 1 ≤ y2 ≤ M

    result = copy(img)
    
    ly1 = max(1, y1-lw)
    ly2 = min(M, y2+lw)
    lx1 = max(1, x1-lw)
    lx2 = min(N, x2+lw)


    result[ly1:y1+lw, x1:x2] .= lc
    result[y1:y2, x2-lw:lx2] .= lc
    result[y2-lw:ly2, x1:x2] .= lc
    result[y1:y2, lx1:x1+lw] .= lc
    return result
end

function cor_image(tomocor::TomoCor, width::Int64 = 5)
    @assert width > 3
    kw = width>>1
    m, M = extrema(tomocor.img)
    p0 = (tomocor.img .- m) ./ (M -m)
    p = RGB.(p0, p0, p0)

    for ys in 1:(size(p)[1])
        xs = round.(Int64, tomocor(ys))
        xi, xf = max(1, xs-kw), min(xs+kw, size(p)[2])
        p[ys, xi:xf] .= RGB(1.0, 0.0, 0)
    end

    for i in eachindex(tomocor.x)
         p[tomocor.y[i], tomocor.x[i]-kw:tomocor.x[i]+kw] .= RGB(0, 1, 0)
    end
    return p
end

