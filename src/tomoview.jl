include("hanaro.jl")

function img2gray(img::Matrix, contrast::Real = 1.0, resize_scale::Real = 1.0)
    @assert 0.0 ≤ contrast ≤ 1.0
    @assert 0.1 ≤ resize_scale ≤ 10.0
    mv, Mv = extrema(img)
    contrast = convert(eltype(img), contrast)
    img =  Gray.((img .- mv)./(Mv-mv) .* contrast)
    if abs(resize_scale-1.0)>0.1 
        resized = [r*resize_scale for r in size(img)]
        img = imresize(img, resized)
    end
end

function get_image(tomo::TomoReader, angle::Real=0.0, contrast::Real = 1.0; 
    threashold::Union{Real, Nothing}=nothing, resize_scale::Real = 1.0)
    
    @assert 0 < contrast <= 1.0
    @assert 0.1 ≤ resize_scale ≤ 10.0

    contrast = Float32(contrast)
    threashold = Float32(threashold)
    ``
    ths = [th for (obj, th, fn) ∈ tomo.data_files if obj == 1]
    if length(ths) == 0 
        @error "No images for object=$object, angle=$angle"
        return 
    end
    
    m = argmin(abs.(angle .- ths))
    img = Float32.(read_nrimage(joinpath(tomo.data_dir, (tomo.data_files[m])[3])))
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

function get_image(tomo::TomoReader, angles::AbstractVector, contrast::Real = 1.0, resize_scale::Real = 1.0) 
    imgs = [get_image(tomo, th, contrast, resize_scale) for th ∈ angles]
    return (imgs |> sum)

end


function get_image(tomo::TomoData, i::Integer , contrast=0.99)
    @assert 1 ≤ i ≤ size(tomo.data)[1]
    return img2gray(tomo.data[i, :, :], contrast)
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

function filtered_back_projection(
    tomo::TomoData, 
    cor::TomoCor,
    index = 1,
    filtername = "hann", 
    interpolation = "nearest", 
    with_threads = true)
    
    nAngles, Nheight, Ndet = size(tomo.data)

    @assert 1 ≤ index ≤ Nheight

    return filtered_back_projection(tomo.data[:, index, :], cor(index), tomo.ths, filtername, interpolation, with_threads)
end


function recon_fbp(
    tomo::TomoData,
    filtername = "hann",
    interpolation = "nearest",
    with_threads = true
    )

    nAngles, Nheight, Ndet = size(tomo.data)

    img_reconed = Array{Float32, 3}(undef, (Ndet, Ndet, Nheight))

    Threads.@threads for i in 1:Nheight
        @inbounds img_reconed[:, :, i] = filtered_back_projection(tomo.data[:, i, :], tomo.ths, tomo_cor(tomo, i), filtername) 
    end
    
    return img_reconed
end

function recon_fbp!(
    tomo::TomoData,
    filtername = "hann",
    interpolation = "nearest",
    with_threads = true
    )

    nAngles, Nheight, Ndet = size(tomo.data)

    img_reconed = Array{Float32, 3}(undef, (Ndet, Ndet, Nheight))

    Threads.@threads for i in 1:Nheight
        @inbounds img_reconed[:, :, i] = filtered_back_projection(tomo.data[:, i, :], tomo.ths, tomo_cor(tomo, i), filtername) 
    end
    
    push!(tomo.process, :recon_fbp)

    tomo.data = img_reconed
end

# function recon_fbp_parts(
#     tomo::TomoData,
#     cor::TomoCor,
#     filtername = "hann",
#     interpolation = "nearest",
#     with_threads = true,
#     nparts::Integer = 2,
#     workingdir = ""
#     )

#     @assert nparts ≥ 2 
#     nAngles, Nheight, Ndet = size(tomo.data)

#     starts_at = [(i*(Nheight ÷ nparts) +1) for i in 0:(nparts-1)]
#     partsinfo = []

#     for j in 1:nparts
#         fn = "recon" * "_" * lpad(j, 4, "0") * ".bin"
#         fpath = joinpath(workingdir, fn)
#         if j ≠ nparts
#             pheight = starts_at[j+1]-starts_at[h]
#             img_reconed = Array{Float32, 3}(undef, (Ndet, Ndet, pheight))
#         else 
#             pheight = Nheights - starts_at[nparts] +1
#         end
        
#         Threads.@threads for i in 1:pheight
#             img_reconed[:,:, i] = filtered_back_projection(tomo.data[:, i, :], cor(i), tomo.ths, filtername, interpolation, with_threads) 
#         end
        
#         write(fpath, fbpimg)
#         push!(partsinfo, (Ndet, Ndet, pheight, starts_at[j], fpath))
#     end
    
#     return partsinfo
#     #push!(tomo.process, :recon_fbp)

#     #tomo.data = img_reconed
# end
