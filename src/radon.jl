using Images, Rotations, ImageTransformations, CoordinateTransformations, Interpolations

function _rotimg(img, θ, c=nothing)
    m, n = size(img)
    if c === nothing
        img_center = (m>>1, n>>1)
    else 
        img_center = (c[1], c[2])
    end
    θ = θ/180.0*π
    mv, Mv = extrema(img)
    timg = Gray.((img.-mv)/(Mv-mv))
    trfm = recenter(RotMatrix(θ), img_center)
    img1 = warp(timg, inv(trfm), method=BSpline(Linear()), fillvalue = Flat(), axes(timg))
    result = Float32.(img1[1:m, 1:n]).*(Mv-mv) .+ mv
    return result
end

function radon(img, angles, c=nothing)
    result = zeros(Float32, (length(angles), size(img)[1]))
    Threads.@threads for i in eachindex(angles)
        rr = _rotimg(img, angles[i], c)
        @inbounds result[i, :] = sum(_rotimg(img, angles[i], c), dims=1)[1,:]
    end
    return result
end