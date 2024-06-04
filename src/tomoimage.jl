abstract type AbstractTomoImage end

struct TomoSlice <:AbstractTomoImage
    theta::Real
    data::Matrix{Real}
end

struct TomoSinogram <:AbstractTomoImage
    data::Matrix{Real}
end

struct TomoVolumeSlice <:AbstractTomoImage
    data::Matrix{Real}
end

function tomoview_makie(p::AbstractTomoImage; size::Tuple{Int64, Int64}=(800, 600))
    f=Figure()
    image(f[1, 1], p.data', axis=(aspect = DataAspect(), yreversed=true, xticklabelsvisible=false, yticklabelsvisible=false))
    return f
end


function makieview(img::Matrix{Real}; size::Tuple{Int64, Int64}=(800, 600))
    f=Figure()
    image(f[1, 1], img', axis=(aspect = DataAspect(), yreversed=true, xticklabelsvisible=false, yticklabelsvisible=false))
    return f
end