using Images, HDF5

tomo_types = (:parallel, :fan, :cone)
tomo_image_types = (:undef, :tif, :tiff, :fits)

instruments = (:undef, :hanaro_nr, :hanaro_enf, :hanaro_enf2)

image_types = Dict(
    :undef => nothing,
    :hanaro_nr => :tif,
    :hanaro_enf => :fits,
    :hanaro_enf2 => :tif,
    )


"""
    get_file_list_hanaro(instrument, image_type, dir_data, dir_white, dir_dark, angle_step)

각각의 디렉토리로부터 데이터, whitebeam, dark 영상 파일의 목록을 작성한다.

하나로 NR 의 경우 파일이름에 회전각이 포함되지 않기 때문에 angle_step 을 이용하여 0.0 도부터 회전각을
부여한다. 하나로 enf 의 경우 파일이름에 회전각 정보가 포함되기 때문에 angle_step 값은 무시한다.
"""
function get_file_list_hanaro(instrument, image_type, dir_data, dir_white, dir_dark, angle_step)    
    
    datafiles = []
    darkfiles= []
    whitefiles = []
    extstr = String(image_type)
    len_ext = length(extstr)
    for fn in readdir(dir_data)
        if (last(fn, len_ext) == extstr) && (~startswith(fn, "._"))
            if instrument == :hanaro_enf
                angle = _rotation_angle_from_filename_enf(fn)
                push!(datafiles, (1, angle, fn))
            elseif instrument == :hanaro_enf2
                angle = _rotation_angle_from_filename_enf2(fn)
                push!(datafiles, (1, angle, fn))
            
            elseif instrument == :hanaro_nr 
                angle = angle_step * parse(Int64, split(fn, "_")[3])
                objstr = split(fn, "_")[1][end]
                object = parse(Int64, objstr)
                #push!(objset, object)
                push!(datafiles, (angle, fn))
            end
        
        end
    end


    for fn in readdir(dir_dark)        
        if (last(fn, len_ext) == extstr) && (~startswith(fn, "._"))
            push!(darkfiles, fn)
        end
    end
    
    for fn in readdir(dir_white)
        
        if (last(fn, len_ext) == extstr) && (~startswith(fn, "._"))
            push!(whitefiles, fn)
        end
    end
    sort!(datafiles, by=(x->x[2]))
    return datafiles, whitefiles, darkfiles
end

"""
    _rotation_angle_from_filename_enf2(fn)

HANARO ENF for high resolution mode uses `tif` format and rotation angle information is written in filename. For example, 
"CICC_20_Camera(0p00)_Rotation(27p0)_1.fits". 27p0 means the rotation angle is 2.70 deg (not 27.0 deg).

It returns the rotation angle in degree.
"""
function _rotation_angle_from_filename_enf2(fn)
    if last(fn, 3) == "tif" 
        pp = split(fn, "_")[end-1][10:end-1]
        angle = 0.1*parse(Int64, split(pp, "p")[1])
        return angle
    else 
        return nothing
    end
end

"""
    _rotation_angle_from_filename_enf(fn)

HANARO ENF's image file use `fits` format and rotation angle information is written in filename. For example, 
"CICC_20_Camera(0p00)_Rotation(27p0)_1.fits". 27p0 means the rotation angle is 2.70 deg (not 27.0 deg).

It returns the rotation angle in degree.
"""
function _rotation_angle_from_filename_enf(fn)
    # println(last(fn, 4))
    if last(fn, 4) == "fits" 
        pp = split(fn, "_")[4][10:end-1]
        angle = 0.1*parse(Int64, split(pp, "p")[1])
        return angle
    else 
        return nothing
    end
end

mutable struct TomoReader
    tomo_type::Symbol
    instrument::Symbol
    data_dir::String
    data_files::AbstractVector
    white_dir::String
    white_files::AbstractVector    
    dark_dir::String
    dark_files::AbstractVector
    working_dir::Union{String, Nothing}
    image_type::Union{Symbol, Nothing}
    to_be_transposed::Bool
    crop_region::Union{Vector{Integer}, Nothing}
    norm_region::Union{Vector{Integer}, Nothing}
    

    function TomoReader(;
        data_dir::String,
        white_dir::String,
        dark_dir::String,
        working_dir::Union{String, Nothing} = nothing,
        tomo_type::Symbol = :parallel,
        instrument::Symbol = :hanaro_nr,
        image_type::Union{Symbol, Nothing} = nothing,
        to_be_transposed::Bool = false,
        angle_step::Real = 0.3,
        crop_region = nothing,
        norm_region = nothing
        )

        @assert tomo_type ∈ tomo_types
        @assert instrument ∈ instruments
        if working_dir === nothing
            @assert prod(isdir.([data_dir, white_dir, dark_dir]))
        else
            @assert prod(isdir.([data_dir, white_dir, dark_dir, working_dir]))
        end
    
        if image_type === nothing
            image_type = image_types[instrument]
        end
    
        datafiles, whitefiles, darkfiles = get_file_list_hanaro(instrument, image_type, data_dir, white_dir, dark_dir, angle_step)
    
        if length(datafiles) == 0
            @error "no data file"
        end
        if length(whitefiles) == 0
            @error "no white file"
        end
        if length(darkfiles) == 0
            @error "no dark file"
        end
    
        return new(tomo_type, instrument, data_dir, datafiles, white_dir, whitefiles, dark_dir, darkfiles, working_dir, image_type, to_be_transposed)
    end
end

function Base.show(io::IO, ::MIME"text/plain", tomo::TomoReader)
    r = "Tomograophy Reader\n"
    r *= "  - Instrumnent : $(String(tomo.instrument))\n"
    r *= "  - Tomography type : $(String(tomo.tomo_type))\n"
    r *= "  - To be transposed : $(tomo.to_be_transposed)\n"
    r *= "  - Number of images : $(length(tomo.data_files))"
    println(io, r)
end

function tomoreaderview(
    tomo::TomoReader, 
    angle_step::Real = 10.0, 
    contrast::Real = 1.0, 
    lw::Integer=5, 
    resize_scale=1.0; 
    threashold::Union{Real, Nothing}=nothing)
    
    @assert 0.0 < contrast <= 1.0

    contrast = Float32(contrast)
    θ, _fn = tomo.data_files[1][2:end]
    img = Float32.(read_nrimage(joinpath(tomo.data_dir, _fn)))
    Nimg = 1
    θ += angle_step
    for obj in tomo.data_files[2:end]
        th, _fn = obj[2], obj[3]
        if th > 180.0
            break
        end
        if th >=θ 
            img += Float32.(read_nrimage(joinpath(tomo.data_dir, _fn)))
            Nimg += 1
            θ += angle_step
            println(_fn, size(img))
        end
    end

    if threashold === nothing 
        img /= (maximum(img)/contrast)
    else 
        img /= (Float32(threashold)/contrast)
        img[img .> 1.0f0] .= 1.0f0
    end
        
    if tomo.to_be_transposed
        img = RGB.(img', img', img')
    else 
        img = RGB.(img, img, img)
    end

    M, N = size(img)
    if tomo.crop_region ≠ nothing
        x1, y1, x2, y2 = tomo.crop_region[1:4]
        ly1 = max(1, y1-lw)
        ly2 = min(M, y2+lw)
        lx1 = max(1, x1-lw)
        lx2 = min(N, x2+lw)

        lc = RGB(1.0, 0.0, 0.0)
        img[ly1:y1+lw, x1:x2] .= lc
        img[y1:y2, x2-lw:lx2] .= lc
        img[y2-lw:ly2, x1:x2] .= lc
        img[y1:y2, lx1:x1+lw] .= lc
    end
    
    if tomo.norm_region ≠ nothing
        x1, y1, x2, y2 = tomo.norm_region[1:4]
        ly1 = max(1, y1-lw)
        ly2 = min(M, y2+lw)
        lx1 = max(1, x1-lw)
        lx2 = min(N, x2+lw)

        lc = RGB(1.0, 1.0, 0.0)
        img[ly1:y1+lw, x1:x2] .= lc
        img[y1:y2, x2-lw:lx2] .= lc
        img[y2-lw:ly2, x1:x2] .= lc
        img[y1:y2, lx1:x1+lw] .= lc
    end

    if abs(resize_scale-1.0)>0.1 
        resized = round.(Int64, [r*resize_scale for r in size(img)])
        img = imresize(img, (resized[1], resized[2]))
    end

    return img

end

"""
        set_crop_region(tomo::TomoReader, x1::Real, y1::Real, x2::Real, y2::Real)

Select the region to be reconstructed. The region is a rectangular defined by two points (x1, y1)
and (x2, y2).
"""
function set_crop_region(tomo::TomoReader, x1::Real, y1::Real, x2::Real, y2::Real)

    @assert 1 ≤ x1 < x2
    @assert 1 ≤ y1 < y2
    tomo.crop_region = [x1, y1, x2, y2]
end


"""
        set_norm_region(tomo::TomoReader, x1::Real, y1::Real, x2::Real, y2::Real)

Select the region for normalization. The region is a rectangular defined by two points (x1, y1)
and (x2, y2).
"""
function set_norm_region(tomo::TomoReader, x1::Real, y1::Real, x2::Real, y2::Real)
    @assert 1 ≤ x1 < x2
    @assert 1 ≤ y1 < y2
    tomo.norm_region = [x1, y1, x2, y2]
end
