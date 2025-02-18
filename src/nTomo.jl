# using Images

module nTomo

using Images

include("fourier_filter.jl")
include("radon.jl")
include("hanaro.jl")
include("utils.jl")
include("phantoms.jl")
include("filtered_back_projection.jl")
include("tomofilter.jl")
include("tomoreader.jl")
include("tomodata.jl")
include("tomocor.jl")
include("tomoview.jl")
include("art.jl")
include("recon.jl")



export read_nrimage,
    fourier_filter2,
    phantom_shepp_logan, 
    _rotimg,
    radon, 

    mat2gray, 
    colorize, 
    image_profile,
    rescale, 
    @h,
    radon, 
    filtered_back_projection,
    filtered_back_projection2,
    get_ordered_indice_by_golden_ratio,
    iradon_sart,
    iradon_fbp,
    
    get_file_list_hanaro,
    
    IdentityFilter,
    MedianFilter,
    CVMedianFilter,
    ThresholdMedianFilter,

    TomoReader,
    set_crop_region,
    set_norm_region,
    tomo_range_view,
    TomoData,
    TomoDataRaw,
    TomoCor,
    get_image,
    get_data_by_index,
    get_image_by_index,
    get_overlapped_image, 
    add_rect,
    area_normalize!,
    detector_normalize!,
    pixel_normalize!,
    sinogram!,
    tomo_cor,
    calc_cor!,
    calc_cor2, # to be deleted
    cor_image,
    rotation_shift!,
    recon_fbp,
    recon_fbp!,
    recon_fbp_parts,
    tomosave,
    tomoload, 
    export_tif 

end # module nTomo
