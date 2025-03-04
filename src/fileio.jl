"""
`sinosave(filepath, sinogram, ths, center)`

Save sinogram and theta angles to hdf5 format. The first axis of sinogram is 
the rotation angle and the second axis is the detector position. The ths can
be Vector or StepRangeLen.

  - `filepath` is the path to save the sinogram.
  - `sinogram` is the sinogram data.
  - `ths` is the theta angles.
  - `center` is the center of rotation. If the center is not given, it 
    is set to -1.0f0 and ignored when reading the file if the center is 
    negative.
""" 

function sinosave(path, sinogram::Matrix, ths::AbstractVector, center::Union{Number, Nothing} = nothing)
    @assert !ispath(path) "File already exists"
    _ext = splitext(path)[2]

    @assert _ext ∈ (".h5", ".hdf5")


    @assert size(sinogram)[1] == length(ths)

    h5open(path, "w") do fid
        create_group(fid, "Sinogram")
        g = fid["Sinogram"]
        g["sinogram"] = sinogram
        g["version"] = [0, 1, 0]
        g["center"] = (center == nothing) ? -1.0f0 : Float32(center)
        if typeof(ths) <:Vector
            g["theta"] = ths
        elseif typeof(ths) <: AbstractRange
            g["theta"] = collect(ths)
        else
            error("Invalid type of ths")
        end
            

        close(fid)
    end
end

"""
`sinoload(filepath)`

Read sinogram file of hdf5 format and return sinogram and ther angles.
"""
function sinoload(path)
    @assert ispath(path)
    _ext = splitext(path)[2]
    @assert _ext ∈ (".h5", ".hdf5")
    
    fid = h5open(path)
    @assert keys(read(fid)) == Set(["Sinogram",])

    _version = read(fid, "Sinogram/version")
    _data = read(fid, "Sinogram/sinogram")
    _theta = read(fid, "Sinogram/theta")
    _center = read(fid, "Sinogram/center")

    if _version[1] == 0
        if _center < 0
            return _data, _theta, nothing
        else 
            return _data, _theta, _center
        end
    else 
        error("Invalid version")
    end
end

"""
`tomosave(filepath, tomo::TomoData)`

Save the TomoData. TomoData can be saved as two format : hdf5(.h5 extension) or 
NRRD (.nrrd extionsion). Although hdf5 format file can be restored to original 
TomoData object, NRRD format file can't be restored. Informations generated in 
data processing is stored only in hdf5 format. 
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
`tomoload(filepath)`

Read TomoData from hdf5 format. The filepath must be ended with .h5
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



"""
`export_tif(tomo, pixeltype, targetdir, filenamebase)`

Export reconstructed tomography data to tif images. `pixeltype` can be `UInt8`
or `UInt16`. `targetdir` is the directory to save the tif files. `filenamebase` 
is the base name of the tif files.

  - `tomo` must be a TomoData object.
  - `pixeltype` must be `UInt8` or `UInt16`.
  - `targetdir` must be a valid directory path of existing directory.
  - `filenamebase` must be a string.
"""
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