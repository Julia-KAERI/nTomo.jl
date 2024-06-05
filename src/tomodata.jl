# include("medianfilters.jl")
include("tomofilter.jl")
tomodata_process = (:dataread, :medianfiltering, :normalization, :sinogram, :pixel_normalization, :cor, :recon_fbp)
reconstruction_method = (:recon_fbp, )
# image_load=(:auto, :memory, :hdf5)


_Version = (0, 1)

abstract type AbstractTomoData end

mutable struct TomoData <: AbstractTomoData
    version
    ths::Vector{Float32}
    data::Array{Float32, 3}
    white::Array{Float32, 3}
    dark::Array{Float32, 3}
    norm_data
    norm_white
    norm_dark
    process
    cor



    function TomoData(
        tomo::TomoReader,
        θ₁::Real = 0.0, 
        θ₂::Real = 180.0,
        imgfilters::Union{AbstractVector{<:AbstractTomoFilter}, Nothing} = nothing;
        )

        isfiltered = false

        if imgfilters ≠ nothing 
            kk = [f.ksize for f in imgfilters]
            if all(x->x==kk[1], kk) === false 
                error("Each of filters have the same size")
            else 
                noise_filter_size = kk[1]
                isfiltered = true
            end
        else 
            imgfilters = [IdentityFilter(3), ]
        end

        @assert isa(tomo.crop_region, Vector)
        @assert isa(tomo.norm_region, Vector)
        
        cx1, cx2 = tomo.crop_region[1], tomo.crop_region[3]
        cy1, cy2 = tomo.crop_region[2], tomo.crop_region[4]

        nx1, nx2 = tomo.norm_region[1], tomo.norm_region[3]
        ny1, ny2 = tomo.norm_region[2], tomo.norm_region[4]

        @assert (cx2-cx1) > noise_filter_size &&  (cy2-cy1) > noise_filter_size
        @assert (nx2-nx1) > noise_filter_size &&  (ny2-ny1) > noise_filter_size
        
        # ksize ≠ filters size. 
        ksize = (noise_filter_size >> 1)

        ths = []
        fns = []
        
        for (obj, θ, fn) in tomo.data_files
            if θ₁ <= θ < θ₂
                push!(ths, θ)
                push!(fns, fn)
            end
        end

        ths = Float32.(ths)

        norm_data = zero(ths)
        norm_white = 0.0f0
        norm_dark = 0.0f0

        # filter size 를 고려하여 잘라내는 이미지 사이즈
        Nx, Ny = cx2-cx1+1 + 2* ksize, cy2-cy1+1 + 2*ksize

        data = Array{Float32}(undef, (length(fns), Ny, Nx))
        white = Array{Float32}(undef, (length(tomo.white_files), Ny, Nx))
        dark = Array{Float32}(undef, (length(tomo.dark_files), Ny, Nx))

        Threads.@threads for i ∈ eachindex(fns)
            fn = fns[i]
            if tomo.to_be_transposed
                pdata = (Float32.(read_nrimage(joinpath(tomo.data_dir, fn))))'
            else 
                pdata = Float32.(read_nrimage(joinpath(tomo.data_dir, fn)))
            end 
            

            ndata = pdata[ny1-ksize:ny2+ksize, nx1-ksize:nx2+ksize]
            pdata = pdata[cy1-ksize:cy2+ksize, cx1-ksize:cx2+ksize]
            
            for i in eachindex(imgfilters)
                ndata = ndata |> imgfilters[i]
                pdata = pdata |> imgfilters[i]
            end

            norm_data[i] = sum(ndata[ksize+1:end-ksize])
            
            data[i, :, :] = pdata[:,:]
        
        end

        for (i, fn) in enumerate(tomo.white_files)
            if tomo.to_be_transposed
                pdata = (Float32.(read_nrimage(joinpath(tomo.white_dir, fn))))'
            else 
                pdata = Float32.(read_nrimage(joinpath(tomo.white_dir, fn)))
            end 
            
            ndata = pdata[ny1-ksize:ny2+ksize, nx1-ksize:nx2+ksize]
            pdata = pdata[cy1-ksize:cy2+ksize, cx1-ksize:cx2+ksize]
            
            for i in eachindex(imgfilters)
                ndata = ndata |> imgfilters[i]
                pdata = pdata |> imgfilters[i]
            end

            norm_white += sum(ndata[ksize+1:end-ksize])
            
            white[i, :, :] = pdata[:,:]
            
        end
            
        norm_white /= length(tomo.white_files)

        for (i, fn) in enumerate(tomo.dark_files)
            if tomo.to_be_transposed
                pdata = (Float32.(read_nrimage(joinpath(tomo.dark_dir, fn))))'
            else 
                pdata = Float32.(read_nrimage(joinpath(tomo.dark_dir, fn)))
            end 
 
            ndata = pdata[ny1-ksize:ny2+ksize, nx1-ksize:nx2+ksize]
            pdata = pdata[cy1-ksize:cy2+ksize, cx1-ksize:cx2+ksize]
            
            for i in eachindex(imgfilters)
                ndata = ndata |> imgfilters[i]
                pdata = pdata |> imgfilters[i]
            end
            
            norm_dark += sum(ndata[ksize+1:end-ksize])
       
            dark[i, :, :] = pdata[:,:]
        end
            
        norm_dark /= length(tomo.white_files)

        
        if isfiltered
            prc = [:dataread, :medianfiltering, ]
        else
            prc = [:dataread]
        end


        return new(_Version, ths, data, white, dark, norm_data, norm_white, norm_dark, prc, [0.0, 0.0])
    end


    function TomoData(_Version, ths, data, white, dark, norm_data, norm_white, norm_dark, prc, cor=[0.0, 0.0])
        return new(_Version, ths, data, white, dark, norm_data, norm_white, norm_dark, prc, cor)
    end
end

"""
    tomosave(filepath, tomo::TomoData)

TomoData can be saved as two format : hdf5(.h5 extension) or NRRD (.nrrd extionsion). Although hdf5 format file
can be restored to original TomoData object, NRRD format file can't be restored. Informations generated in data
processing is stored only in hdf5 format. 
"""
function tomosave(filepath, tomo::TomoData)
    @assert !ispath(filepath)
    _ext = splitext(filepath)[2]

    @assert _ext ∈ (".h5", ".nrrd")

    if _ext == ".h5"
        fid = h5open(filepath, "w")
        create_group(fid, "TomoData")
        g = fid["TomoData"]
        g["version"] = [t for t in tomo.version]
        g["ths"] = tomo.ths
        g["data"] =  tomo.data
        g["white"] =  tomo.white
        g["dark"] = tomo.dark
        g["norm_data"] = tomo.norm_data 
        g["norm_white"] =  tomo.norm_white 
        g["norm_dark"] =  tomo.norm_dark

        g["process"] = [String(t) for t in tomo.process]
        g["cor"] = tomo.cor

        close(fid)
    else 
        # save as nrrd format
        save(filepath, tomo.data)
    end
end

"""
    tomoload(filepath)

load TomoData from hdf5 format. The filepath must be ended with .h5
"""
function tomoload(filepath)
    _ext = splitext(filepath)[2]
    
    @assert _ext ∈ [".h5", ]

    fid = h5open(filepath)
    _version = Tuple(read(fid, "TomoData/version"))
    _ths = read(fid, "TomoData/ths")
    _data = read(fid, "TomoData/data")
    _white = read(fid, "TomoData/white")
    _dark = read(fid, "TomoData/dark")
    _norm_data = read(fid, "TomoData/norm_data")
    _norm_white = read(fid, "TomoData/norm_white")
    _norm_dark = read(fid, "TomoData/norm_dark")
    _proc = [Symbol(t) for t in read(fid, "TomoData/process")]
    _cor = read(fid, "TomoData/cor")
    return TomoData(_version, _ths, _data, _white, _dark, _norm_data, _norm_white, _norm_dark, _proc, _cor)
end



function resize!(tomo::TomoData, factor::Integer = 1)
    @assert :normalization ∈ tomo.process
    @assert :sinogram ∉ tomo.process
    @assert factor ∈ 1:5
    L, M, N = size(tomo.data)
    
    if factor == 1
        return
    end


    _data = Array{eltype(tomo.data)}(undef, (L, M÷3, N÷3))

    @Threads.threads for l ∈ 1:L
        @inbounds _data[l, :, :] = rescale(tomo.data[l, :, :], factor)
    end

    tomo.data = _data

    L, M, N = size(tomo.white)
    _white = Array{eltype(tomo.white)}(undef, (L, M÷3, N÷3))
    @Threads.threads for l ∈ 1:L
        @inbounds _white[l, :, :] = rescale(tomo.white[l, :, :], factor)
    end

    tomo.white = _white

    L, M, N = size(tomo.dark)
    _dark = Array{eltype(tomo.dark)}(undef, (L, M÷3, N÷3))
    @Threads.threads for l ∈ 1:L
        @inbounds _dark[l, :, :] = rescale(tomo.dark[l, :, :], factor)
    end

    tomo.dark = _dark
end

function area_normalize!(tomo::TomoData)
    @assert :medianfiltering ∈ tomo.process
    @assert :normalization ∉ tomo.process

    @Threads.threads for i in 1:size(tomo.data)[1]
        tomo.data[i, :, :] = tomo.data[i, :, :]./(tomo.norm_data[i]/tomo.norm_white)
    end
    push!(tomo.process, :normalization)
end

function sinogram!(tomo::TomoData, remove_zero_fluct = true, zero_fluct_factor = 0.01)
    @assert :normalization ∈ tomo.process
    @assert :sinogram ∉ tomo.process

    # mhiwte, mdark : mean of white and dark datas
    mwhite = mean(tomo.white, dims=1)[1, :, :]
    mdark = mean(tomo.dark, dims=1)[1, :, :]

    @Threads.threads for i in 1:(size(tomo.data)[1])
        tomo.data[i,:,:] = -log.((tomo.data[i,:,:] .- mdark)./(mwhite .- mdark))
    end

    if remove_zero_fluct
        zero_fluct_factor = Float32(zero_fluct_factor)
        @assert  0.0f0< zero_fluct_factor < 1.0f0
        factor = zero_fluct_factor * maximum(tomo.data)
        
        tomo.data[tomo.data .< factor] .= 0.0
    end

    push!(tomo.process, :sinogram)

end

function pixel_normalize(sino::Matrix{Float32}, Ns=2, lv=1.2)
    psum = sum(sino, dims = 1)
    sino_corrected = zero(sino)
    for i in (1+Ns):(size(sino)[2]-Ns)
        nmean = mean(psum[i-Ns:i+Ns])
        if abs(psum[i] * nmean ) < 1.0e-1

        elseif abs(psum[i] - nmean) > lv 
            sino_corrected[:, i] = sino[:, i] .* (nmean/psum[i])
        else 
            sino_corrected[:, i] = sino[:, i]
        end
    end
    return sino_corrected
end

function pixel_normalize!(tomo::TomoData, Ns=2, lv=1.2)
    @assert :pixel_normalize ∉ tomo.process
    @assert :sinogram ∈ tomo.process
    for index in 1:size(tomo.data)[2]
        sino = tomo.data[:, index, :]
        psum = sum(sino, dims = 1)
        sino_corrected = zero(sino)


        for i in (1+Ns):(size(sino)[2]-Ns)
            nmean = mean(psum[i-Ns:i+Ns])
            if abs(psum[i] * nmean ) < 1.0e-1

            elseif abs(psum[i] - nmean) > lv 
                sino_corrected[:, i] = sino[:, i] .* (nmean/psum[i])
            else 
                sino_corrected[:, i] = sino[:, i]
            end
        end
        tomo.data[:, index, :] = sino_corrected
    end
    push!(tomo.process, :pixel_normalization)
end


function calc_cor!(tomo::TomoData, y1=nothing, y2=nothing, return_img = false)
    @assert :sinogram ∈ tomo.process
    
    M, N = size(tomo.data)[2:3]

    c1, c2 = ceil(Int64, N/3), floor(Int64, N*2/3)
    
    if y1 === nothing
        y1 = 1
    end

    if y2 === nothing
        y2 = M
    elseif y2 < 0
        y2 = M-y2
    end

    y1, y2 = minmax(y1, y2)
    yi, yf = ceil(Int64, y1), floor(Int64, y2)
    @assert 1 ≤ yi < yf ≤ M

    p = tomo.data[1, yi:yf, :] + tomo.data[end, yi:yf, :]
    cor = []
    
    for y in yi:1:yf
        tt = []      
        for c in c1:1:c2
            ts = 0
            if (c <= (c1+c2)*0.5)
                i0, i1 = 1, 2*c-1
            else 
                i0, i1 = 2*c-N, N  
            end

            for x in (i0+1):1:(i1-1)
                ts += (p[y-yi+1, x]-p[y-yi+1, i1-x+i0])^2
            end
            push!(tt, ts) 
        end
        push!(cor, [y, findmin(tt)[2]+c1-1])  
    end
    
    xyarr = reduce(hcat, cor)'
    y, x = xyarr[:,1], xyarr[:,2]
    
    A = [sum(y.^2) sum(y) ; sum(y) length(y)]
    B = [sum(y.*x) ; sum(x)]
    C = inv(A)*B
    tomo.cor = [C[1], C[2]]

    if :cor ∉ tomo.process
        push!(tomo.process, :cor)
    end

    if return_img
        
        kw = 2 # image 에서의 중앙을 나타내는 선의 굵기

        img = tomo.data[1, :, :] + tomo.data[end, :, :]
        m, M = extrema(img)
        p0 = (img .- m) ./ (M -m)
        p = RGB.(p0, p0, p0)

        for ys in 1:(size(p)[1])
            xs = round.(Int64, C[1]*ys+C[2])
            xi, xf = max(1, xs-kw), min(xs+kw, size(p)[2])
            p[ys, xi:xf] .= RGB(1.0, 0.0, 0) # Green color
        end

        for i in eachindex(x)
            p[y[i], x[i]-kw:x[i]+kw] .= RGB(0, 1, 0) #Red color
        end
        return p
    end
end

"""
    cor(tomo::TomoData, y::Real)

return center of rotation at y (vertical axis)
"""
function tomo_cor(tomo::TomoData, y::Real)
    if sum(abs.(tomo.cor)) < 1.0 
        return nothing

    else 
        return tomo.cor[1]*y+tomo.cor[2]
    end
end


function tomo_cor(cor::Vector, y::Real)
    @assert length(cor)==2
    return cor[1]*y + cor[2]
end


function export_tif(tomo::TomoData, pixeltype, targetdir::String="", filenamebase="recon")
    @assert :recon_fbp ∈ tomo.process
    @assert pixeltype ∈ (UInt8, UInt16) 
    if length(targetdir) > 0
        @assert isdir(targetdir)
    end

    mv, Mv = extrema(tomo.data)
    for i in 1:size(tomo.data)[3]
        p = round.(pixeltype, (tomo.data[:,:, i] .- mv) .* (typemax(pixeltype)/(Mv - mv)))
        fn = filenamebase * "_" * lpad(i, 4, "0") * ".tif"
        save(joinpath(targetdir, fn), p)
    end
end