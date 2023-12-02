using Images

"""
    phantom_shepp_logan(n, padding=0)

Output modified shepp logan phantom with zero padding.  Shepp logan phantom is taken from the shepp_logan() function 
in `Images.jl`. The total image size is (n+2*padding, n+2*padding)
"""
function phantom_shepp_logan(n, padding=0)
    img0 = Float32.(shepp_logan(n))
    img = zeros(Float32, (n+2*padding, n+2*padding))
    img[padding+1:n+padding, padding+1:n+padding] = img0[:, :]
    return img
end

