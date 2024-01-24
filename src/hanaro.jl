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
    read_nrimage(filepath)

rean Hanaro NR or ENF file to return Matrix{UInt16}
"""
function read_nrimage(filepath)
    extension = lowercase(last(split(filepath, ".")))
    if extension ∈ ("tif", "tiff")
        # return readtif(filepath)[1]["data"]
        return read_tif(filepath)

    elseif extension ∈ ("fits",)
        f1 = fits_open_file(filepath)
        fsize = fits_get_img_size(f1)
        dd = zeros(UInt16, (fsize[1], fsize[2], fsize[3]))
        cc= fits_read_pix(f1, dd)
        return dd[:,:,1]'
    else
        error("Input should be tif or fits file")
    end
end

