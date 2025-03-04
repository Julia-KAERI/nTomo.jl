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
                angle = angle_step * (parse(Int64, split(fn, "_")[3]) -1)
                objstr = split(fn, "_")[1][end]
                object = parse(Int64, objstr)
                #push!(objset, object)
                push!(datafiles, (1, angle, fn))
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

"""
    TomoReader

Tomography 데이터를 읽고 처리하기 위한 자료구조. 데이터, 다크, 화이트 이미지는 각각 다른 디렉토리에 분리되어 있어야 한다. 
본격적으로 데이터를 처리하기 전에 아래와 같은 일을 수행 할 수 있다.
    
- `set_crop_region` 함수를 통해 전체 이미지에서 reconstruction 에 사용될 이미지 영역을 정한다.
- `set_norm_region` 함수를 통해 normalization 할 영역을 정한다.
- 이미지 일부를 이용하여 median filter 를 어떤 크기로 얼마나 사용할지를 정한다.

Fields
======

    tomo_type           : (:parallel, :pan, :cone) 중의 하나.
    instrument          : 장치 이름
    data_dir            : 데이터 디렉토리 경로
    data_files          : 각 이미지 파일이 (object, 각도(degree), 파일이름) 들의 벡터로 저장된다.
    white_dir           : 화이트빔 데이터 디렉토리 경로
    white_files         : 화이트빔 이미지 파일 이름들의 벡터
    dark_dir            : 다크 데이터 디렉토리 경로
    dark_files          : 다크 이미지 파일 이름들의 벡터
    working_dir         : 작업 디렉토리, 현재로서는 별 기능을 하지 않는다.
    image_type          : 이미지 형식. 현재로서는 기능을 하지 않는다. 하나로 장치의 겨우 instrument 에 따라 
                            정해진다.
    to_be_transposed    : 계산을 위해 회전축은 이미지의 세로방향이어야 한다. 가로방향일경우 `true`.
    crop_region         : crop 영역
    norm_region         : normalizaion 영역
    scale_down          : 1, 2, 3, 4, 5 중의 하나. 1/scale_down 으로 가로/세로 크기가 작아진다. 
"""
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
    scale_down::Int64
    

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
        norm_region = nothing,
        scale_down::Int64 = 1
        )

        @assert tomo_type ∈ tomo_types
        @assert instrument ∈ instruments
        @assert scale_down ∈ 1:5
        
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
        
        return new(tomo_type, instrument, data_dir, datafiles, white_dir, whitefiles, dark_dir, darkfiles, working_dir, image_type, 
                    to_be_transposed, nothing, nothing, scale_down)
    end
end

function Base.show(io::IO, ::MIME"text/plain", tomo::TomoReader)
    r = "Tomograophy Reader\n"
    r *= "  - Instrumnent : $(String(tomo.instrument))\n"
    r *= "  - Tomography type : $(String(tomo.tomo_type))\n"
    r *= "  - To be transposed : $(tomo.to_be_transposed)\n"
    r *= "  - Number of images : $(length(tomo.data_files))\n"
    r *= "  - Scale down factor : $(tomo.scale_down)"
    println(io, r)
end

function tomo_range_view(
    tomo::TomoReader, 
    angle_step::Real = 10.0;
    contrast::Real = 1.0,
    lw::Integer=2,  
    threashold::Union{Real, Nothing}=nothing)
    
    @assert 0.0 < contrast <= 1.0

    contrast = Float32(contrast)
    θ, _fn = tomo.data_files[1][2:end]
    img = Float32.(read_nrimage(joinpath(tomo.data_dir, _fn), tomo.scale_down ))
    Nimg = 1
    θ += angle_step
    for obj in tomo.data_files[2:end]
        th, _fn = obj[2], obj[3]
        if th > 180.0
            break
        end
        if th >=θ 
            img += Float32.(read_nrimage(joinpath(tomo.data_dir, _fn), tomo.scale_down))
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

"""
    set_scale_down(tomo::TomoReader, factor::Integer=1)

Change ths scale_down factor
"""
function set_scale_down(tomo::TomoReader, factor::Integer=1)
    @assert factor ∈ 1:5
    tomo.scale_down = factor
end

