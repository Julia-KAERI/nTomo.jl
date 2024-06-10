using FileIO, CFITSIO, Images

function bin2str(arr, endian)
    return String(Char.(arr))
end

function read_tif(filepath)
    #print("read_tif")
    p1 = load(filepath)
    return broadcast(q->UInt16(q.val.i),p1)
end

"""
    read_nrimage(filepath, scald_down)

Rean Hanaro NR or ENF file to return Matrix{UInt16}. The image is resized to the factor of 1/scale_down
"""
function read_nrimage(filepath, scale_down::Integer = 1)
    @assert scale_down ∈ 1:5
    extension = lowercase(last(split(filepath, ".")))
    if extension ∈ ("tif", "tiff")
        # return readtif(filepath)[1]["data"]
        return rescale(read_tif(filepath), scale_down)

    elseif extension ∈ ("fits",)
        f1 = fits_open_file(filepath)
        fsize = fits_get_img_size(f1)
        dd = zeros(UInt16, (fsize[1], fsize[2], fsize[3]))
        cc= fits_read_pix(f1, dd)
        return rescale(dd[:,:,1]', scale_down)
    else
        error("Input should be tif or fits file")
    end
end

