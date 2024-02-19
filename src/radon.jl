using Images, Rotations, ImageTransformations, CoordinateTransformations, Interpolations

"""
    _rotimg(img, θ, c)

Rotate the image around c by θ. θ is in degrees.
"""
function _rotimg(img::Matrix{T}, θ::AbstractFloat , c=Union{Real, Nothing}=nothing) where T<:AbstractFloat
    m, n = size(img)
    if c === nothing
        img_center = (m>>1, n>>1)
    else 
        img_center = (c[1], c[2])
    end
    θ = θ/180.0*π
    mv, Mv = extrema(img)
    timg = Gray.((img.-mv)/(Mv-mv))
    trfm = recenter(RotMatrix(-θ), img_center)
    img1 = warp(timg, inv(trfm), method=BSpline(Linear()), fillvalue = Flat(), axes(timg))
    result = T.((img1[1:m, 1:n]).*(Mv-mv) .+ mv)
    return result
end

"""
    radon(img, angles, center)

Output sinogram by radon transform. The axes of sinograms are (angles, pixels).
"""
function radon(img::Matrix{T}, angles, center::Union{Real, Tuple, Nothing}=nothing) where T<:AbstractFloat
    result = zeros(T, (length(angles), size(img)[1]))
    Threads.@threads for i in eachindex(angles)
        rr = _rotimg(img, angles[i], center)
        @inbounds result[i, :] = sum(_rotimg(img, angles[i], center), dims=1)[1,:]
    end
    return result
end