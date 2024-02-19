"""
Compute the projection of an image along a ray.

Parameters
----------
image : Matrix, shape (M, N)
    Image to project.
θ : Radian angle of the projection. 
ray_position : float
    Position of the ray within the projection.

Returns
-------
projected_value : float
    Ray sum along the projection.
norm_of_weights :
    A measure of how long the ray's path through the reconstruction circle was.
"""
function bilinear_ray_sum(
    image::Matrix{<:Real}, 
    θ, 
    x0::Integer, 
    center)
    
    cx, cy = center[1:2]
    M, N = size(image)
    radius = (M >>1) -1
    
    # (s, t) is the (x, y) system rotated by θ
    t = x0 - cx
    # s0 is the half-length of the ray's path in the reconstruction circle
    s0 = sqrt(max(radius^2 - t^2, 0.0f0))
    Ns = 2 * ceil(Int64, 2 * s0)  # number of steps along the ray
    ray_sum = 0.0f0
    weight_norm = 0.0f0

    cosθ, sinθ = cos(θ), sin(θ)
    if Ns > 0
        # step length between samples
        ds = 2 * s0 / Ns
        dy = -ds * cosθ
        dx = -ds * sinθ
        # point of entry of the ray into the reconstruction circle
        y0 = s0 * cosθ - t * sinθ
        x0 = s0 * sinθ + t * cosθ
        for k in 0:Ns
            x = x0 + k * dx
            y = y0 + k * dy
            index_i = y + cy
            index_j = x + cx
            i = floor(Int64, index_i)
            j = floor(Int64, index_j)
            di = index_i - i
            dj = index_j - j
            # Use linear interpolation between values
            # Where values fall outside the array, assume zero
            
            if 1 < i ≤ M  && 1 < j  ≤ N
                weight = (1. - di) * (1. - dj) * ds
                ray_sum += weight * image[i, j]
                weight_norm += weight * weight
            end
            if 1 < i ≤ M && 1 ≤ j < N
                weight = (1. - di) * dj * ds
                ray_sum += weight * image[i, j+1]
                weight_norm += weight * weight
            end
            if 1≤ i < M && 1 < j ≤ N
                weight = di * (1 - dj) * ds
                ray_sum += weight * image[i+1, j]
                weight_norm += weight * weight
            end
            if 1 ≤ i < M && 1 ≤ j < N
                weight = di * dj * ds
                ray_sum += weight * image[i+1, j+1]
                weight_norm += weight * weight
            end
        end
    end

    return ray_sum, weight_norm
end

"""
    bilinear_ray_update!(image, image_update, θ, x0, center, projected_value)

    Compute the update along a ray using bilinear interpolation.

Parameters
----------
image   : ndarray of float, shape (M, N)
    Current reconstruction estimate.
image_update : ndarray of float, shape (M, N)
    Array of same shape as ``image``. Updates will be added to this array.
θ   : Radian angle of the projection.
x0  : Position of the ray within the projection.
projected_value : float
    Projected value (from the sinogram).

Returns
-------
deviation :
    Deviation before updating the image.
"""

function bilinear_ray_update!(
    image::Matrix{<:Real},
    image_update::Matrix{<:Real},
    θ::Real, 
    x0,
    center, 
    projected_value)

    cx, cy = center[1:2]
    M, N = size(image)

    ray_sum, weight_norm = bilinear_ray_sum(image, θ, x0, center)
    
    if weight_norm > 0
        deviation = -(ray_sum - projected_value) / weight_norm
    else
        deviation = 0.0f0
    end
    radius = M>>1
    # (s, t) is the (x, y) system rotated by theta
    t = x0 - cx
    # s0 is the half-length of the ray's path in the reconstruction circle
    
    s0 = (radius*radius ≥ t*t) ? sqrt(radius*radius - t*t) : 0.0f0
    
    Ns = 2 * ceil(Int64, 2 * s0)
    
    # beta for equiripple Hamming window
    hamming_beta = 0.46164f0
    
    if Ns > 0
        # Step length between samples
        ds = 2 * s0 / Ns
        dy = -ds * cos(θ)
        dx = -ds * sin(θ)
        # Point of entry of the ray into the reconstruction circle
        y0 = s0 * cos(θ) - t * sin(θ)
        x0 = s0 * sin(θ) + t * cos(θ)
        for k ∈ 0:Ns
            x = x0 + k * dx
            y = y0 + k * dy
            index_i = y + cy
            index_j = x + cx
            i = floor(Int64, index_i)
            j = floor(Int64, index_j)
            di = index_i - i
            dj = index_j - j
            hamming_window = ((1 - hamming_beta)
                                - hamming_beta * cos(2 * π * k / (Ns - 1)))
            if 1 < i ≤ M && 1 < j ≤ N
                image_update[i, j] += (deviation * (1. - di) * (1. - dj) * ds * hamming_window)
            end
            if 1 < i ≤ M && 1 < j < N
                image_update[i, j+1] += (deviation * (1. - di) * dj * ds * hamming_window)
            end
            if 1 ≤ i < M && 1 < j ≤ N
                image_update[i+1, j] += (deviation * di * (1 - dj) * ds * hamming_window)
            end
            if 1 ≤ i < M &&  1 ≤ j < N
                image_update[i+1, j+1] += (deviation * di * dj * ds * hamming_window)
            end
        end
    end

    return deviation
end 

"""
    sart_projection_update(image, θ, projection, center)

Compute update to a reconstruction estimate from a single projection
using bilinear interpolation.

Parameters
----------
image : Matrix for current reconstruction estimate
θ : Radian angle of the projection
projection : ndarray of float, shape (P,)
    Projected values, taken from the sinogram

Returns
-------
image_update : ndarray of float, shape (M, N)
    Array of same shape as ``image`` containing updates that should be
    added to ``image`` to improve the reconstruction estimate
"""
function sart_projection_update(
    image,
    θ,
    projection::Array{Float32}, 
    center)

    M = size(projection)[1]
    image_update = zero(image)
    for i in 1:(size(projection)[1])
        bilinear_ray_update!(image, image_update, θ, i, center, projection[i])
    end
    return image_update
end

"""
    get_ordered_indice_by_golden_ratio(θs, interval = 180.0)

Order angles to reduce the amount of correlated information in
subsequent projections.

Parameters
----------
ths : Projection angles in degrees. Duplicate angles are not allowed.
interval : 

Returns
-------
indices_generator : generator yielding unsigned integers
    The returned generator yields indices into ``theta`` such that
    ``theta[indices]`` gives the approximate golden ratio ordering
    of the projections. In total, ``len(theta)`` indices are yielded.
    All non-negative integers < ``len(theta)`` are yielded exactly once.

Notes
-----
The method used here is that of the golden ratio introduced by T. Kohler[1].

References
----------
[1] Kohler, T. "A projection access scheme for iterative
    reconstruction based on the golden section." Nuclear Science
    Symposium Conference Record, 2004 IEEE. Vol. 6. IEEE, 2004.

[2] Winkelmann, Stefanie, et al. "An optimal radial profile order
    based on the Golden Ratio for time-resolved MRI."
    Medical Imaging, IEEE Transactions on 26.1 (2007): 68-76.

"""
function get_ordered_indice_by_golden_ratio(θs, interval = 180.0)

    γ = (2.0/(√5-1))^2
    # indices into theta yield an arbitrary angle to start things off
    ths = Array(θs)
    indices_remained = sortperm(ths)
    angle = ths[1]
    ordered_indices = Array{Int64}([])
    first_idx = 1
    append!(ordered_indices, indices_remained[first_idx])
    popat!(indices_remained, first_idx)
    # determine subsequent angles using the golden ratio method
    angle_increment = interval / γ
    
    while length(indices_remained)>0
        remaining_angles = ths[indices_remained]
        tangle = (angle + angle_increment) % interval
        
        v, idx = findmin(abs.(ths[indices_remained] .- tangle))
        angle = ths[indices_remained[idx]]
        append!(ordered_indices, indices_remained[idx])
        popat!(indices_remained, idx)
        
    end
    return ordered_indices
end


"""
    iradon_sart(sinogram, ths; center=nothing, relaxation=0.15)

Reconstruct an image from the radon transform, using a single iteration of
the Simultaneous Algebraic Reconstruction Technique (SART) algorithm.

Parameters
----------
sinogram : Matrix with (angles, pixels) size.
ths : angles of siniogram.
center(optional) : center of rotation. If not given, center of pixels is automatically 
    asigned.
relaxation(optional) : Relaxation parameter for the update step. A higher value can improve 
    the convergence rate, but one runs the risk of instabilities. Values close to or 
    higher than 1 are not recommended.


Notes
-----
Algebraic Reconstruction Techniques are based on formulating the tomography
reconstruction problem as a set of linear equations. Along each ray,
the projected value is the sum of all the values of the cross section along
the ray. A typical feature of SART (and a few other variants of algebraic
techniques) is that it samples the cross section at equidistant points
along the ray, using linear interpolation between the pixel values of the
cross section. The resulting set of linear equations are then solved using
a slightly modified Kaczmarz method.

When using SART, a single iteration is usually sufficient to obtain a good
reconstruction. Further iterations will tend to enhance high-frequency
information, but will also often increase the noise.

References
----------
[1] AC Kak, M Slaney, "Principles of Computerized Tomographic
    Imaging", IEEE Press 1988.
[2] AH Andersen, AC Kak, "Simultaneous algebraic reconstruction
    technique (SART): a superior implementation of the ART algorithm",
    Ultrasonic Imaging 6 pp 81--94 (1984)
[3] S Kaczmarz, "Angenäherte auflösung von systemen linearer
    gleichungen", Bulletin International de l’Academie Polonaise des
    Sciences et des Lettres 35 pp 355--357 (1937)
[4] Kohler, T. "A projection access scheme for iterative
    reconstruction based on the golden section." Nuclear Science
    Symposium Conference Record, 2004 IEEE. Vol. 6. IEEE, 2004.
[5] Kaczmarz' method, Wikipedia,
    https://en.wikipedia.org/wiki/Kaczmarz_method
"""
function iradon_sart(
    sinogram::Matrix{<:Real},
    ths::AbstractVector{<:Real};
    image::Union{Matrix{<:Real}, Nothing} = nothing,
    center::Union{Real, Nothing} = nothing,
    relaxation = 0.15)

    Mtype = eltype(sinogram)

    M, N = size(sinogram)
    if center === nothing 
        center = (N>>1, N>>1)
    else 
        center = (center, N>>1)
    end
    @assert M == length(ths)

    if image === nothing
        image = zeros(Mtype, (N, N))
    else 
        @assert size(image) == (N, N)
    end

    ordered_indices = get_ordered_indice_by_golden_ratio(ths)

    for angle_index ∈ ordered_indices
        # println("angle = $(ths[angle_index]), index = $angle_index")
        image_update = sart_projection_update(
            image, 
            ths[angle_index] * π / 180.0,
            sinogram[angle_index, :], 
            center)
        image .= image .+ (relaxation .* image_update)
    end
    
    return image
end