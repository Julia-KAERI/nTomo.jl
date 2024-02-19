using FFTW

include("fourier_filter.jl")


function fbp_preproc(sinogram, center, Ndet)
    nAngles, Ndet = size(sinogram)

    N = max(64, 2^ceil(Int64, log2(2 * Ndet)))

    # angles = Float32.(-  .* (pi/180.0) .+ pi/2.0)
    S = zeros(Float32, (N, nAngles))
    xshift = round(Int32, center-Ndet/2)
    #xshift = 0.0f8
    if xshift == 0
        S[1:Ndet, :] = sinogram'
    elseif xshift>0
        S[1:Ndet-xshift, :] = sinogram[:, xshift+1:Ndet]'
    else 
        S[-xshift+1:Ndet, :] = sinogram[:, 1:Ndet+xshift]'
    end

    return S
end

function maincal_cpu(S, angles, ffilter, Ndet, center)
    I = zeros(Float32,(Ndet, Ndet))
    x = Float32.(range(-0.5, stop=0.5, length = Ndet) .+ (center/Ndet - 0.5))
    dx = (center - Ndet/2)
    Threads.@threads for t in eachindex(angles)
        A = real(ifft(fft(S[:, t]).*ffilter))[1:Ndet]
        cosθ, sinθ = cos(angles[t]), sin(angles[t])
        @inbounds for n in 1:Ndet
            yp = x[n]
            @inbounds for m in 1:Ndet
                xp = x[m] 
                rn = (xp*cosθ+yp*sinθ + 0.5)*Ndet 

                ind = round(Int64, rn)
                if 1 ≤ ind ≤ Ndet
                    m2 = round(Int64, m+dx)
                    #m2=m
                    n2 = round(Int64, n+2*dx)
                    if (1 ≤ m2 ≤ Ndet) && (1≤ n2 ≤ Ndet)
                        @inbounds I[m2, n2] += A[ind]
                    end
                    
                end
            end
        end
    end
    return I/pi/ (2 * length(angles))
end


# function maincal_cuda!(Sa, angles, x, Ndet, Im)

#     t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     i = (blockIdx().y - 1) * blockDim().y + threadIdx().y
#     j = (blockIdx().z - 1) * blockDim().z + threadIdx().z

#     nAngles = length(angles)
#     cosθ, sinθ = cos(angles[t]), sin(angles[t])

#     rn = (x[i]*cosθ + x[j]*sinθ + 0.5)*(Ndet)

#     ind = round(Int64, rn)
#     if (ind > 0) && (ind<=Ndet) && (i≤ Ndet) && (j ≤ Ndet) && (t ≤ nAngles) 
#         CUDA.@atomic Im[i, j] += Sa[ind, t]
#     end
#     return nothing
# end

# function filtered_back_projection(sinogram, angles, center, filtername, mode=:cpu)
#     @assert mode ∈ (:cpu, :cuda)
#     if filtername ∉ filter_name
#         error("filtername shouled be one of $(filter_name)")
#     end
#     nAngles, Ndet = size(sinogram)
    
#     S = fbp_preproc(sinogram, center, Ndet)
#     N = size(S)[1]
#     ffilter = fourier_filter2(N, filtername)
#     angles = Float32.(-angles .* (pi/180.0) .+ pi/2.0)

#     if mode == :cpu
#         img = maincal_cpu(S, angles, ffilter, Ndet, center)
#     else 
#         error("no cuda yet")

#         # Sa = real(ifft(fft(CuArray(S), 1) .* repeat(CuArray(ffilter), 1, length(angles)), 1))[1:Ndet, :]


#         # Im = CUDA.zeros(Float32,(Ndet, Ndet))

#         # xx = CuArray(Float32.(range(-0.5, stop=0.5, length = Ndet) .+ (center/Ndet - 0.5)))

#         # nthreads = (8, 8, 8)
#         # nsize = (length(angles), Ndet, Ndet)
#         # nblocks = ceil.(Int, nsize .÷ nthreads)     
#         # @cuda threads=nthreads blocks=nblocks maincal_cuda!(Sa, CuArray(angles), xx, Ndet, Im)
#         # img = Im./(π*2*length(angles))
#     end
# end

function iradon_fbp(
    sinogram::Matrix{<:Real},
    angles::AbstractVector{<:Real},
    center::Union{Real, Nothing} = nothing,
    filtername::String = "hann")

    if filtername ∉ filter_name
        error("filtername shouled be one of $(filter_name)")
    end
    nAngles, Ndet = size(sinogram)
    
    S = fbp_preproc(sinogram, center, Ndet)
    N = size(S)[1]
    ffilter = fourier_filter2(N, filtername)
    angles = Float32.(angles .* (pi/180.0) .+ pi/2.0)

    
    return maincal_cpu(S, angles, ffilter, Ndet, center)
    
end


function filtered_back_projection2(
    sinogram::Matrix{Float32},
    ths,
    center,
    filtername = "hann",
    interpolation = "nearest",
    with_threads = true)
    
    nAngles, Ndet = size(sinogram)

    N = max(64, 2^ceil(Int64, log2(2 * Ndet)))

    angles = Float32.(- ths .* (pi/180.0) .+ pi/2.0)
    S = zeros(Float32, (N, nAngles))

    N, nAngles = size(S)

    radius = Float32(max(center, N-center))

    xshift = round(Int32, center-Ndet/2)
    # shift sinogram along spatial axis to coincide the center of rotation and 
    # the center of sinogram.
    if xshift == 0
        S[1:Ndet, :] = sinogram'
    elseif xshift>0
        S[1:Ndet-xshift, :] = sinogram[:, xshift+1:Ndet]'
    else 
        S[-xshift+1:Ndet, :] = sinogram[:, 1:Ndet+xshift]'
    end

    if interpolation ∉ ("nearest", "linear")
        error("$(interpolation) is not supperted interpolation method")
    end

    I = zeros(Float32,(Ndet, Ndet))

    filtername = lowercase(filtername)

    filter = fourier_filter2(N, filtername)
    if filter === nothing
        error("$filtername is not supported")
    end

    x = Float32.(range(-0.5, stop=0.5, length = Ndet) .+ (center/Ndet - 0.5))

    if with_threads
        Threads.@threads for t in eachindex(angles)
            A = real(ifft(fft(S[:, t]).*filter))[1:Ndet]
            cosθ, sinθ = cos(angles[t]), sin(angles[t])
            @inbounds for n in 1:Ndet
                yp = x[n]
                @inbounds for m in 1:Ndet
                    xp = x[m]
                    rn = (xp*cosθ+yp*sinθ + 0.5)*Ndet

                    if interpolation == "nearest"
                        ind = round(Int64, rn)
                        if (ind > 0) && (ind<=Ndet)
                            @inbounds I[m, n] += A[ind]
                        end
                    elseif interpolation == "linear"
                        ind = floor(Int64, rn)
                        if (ind>0) && (ind<Ndet)
                            @inbounds I[m, n] += (rn - ind)*(A[ind+1]-A[ind])  + A[ind]
                        end
                    end
                end
            end
        end

    else 
        for t in eachindex(angles)
            A = real(ifft(fft(S[:, t]).*filter))[1:Ndet]
            @inbounds for t in 1:nAngles
                cosθ, sinθ = cos(angles[t]), sin(angles[t])
                @inbounds for n in 1:Ndet
                    yp = x[n]
                    @inbounds for m in 1:Ndet
                        xp = x[m]
                        rn = (xp*cosθ - yp*sinθ + 0.5)*Ndet              

                        if interpolation == "nearest"
                            ind = round(Int64, rn)
                            if (ind > 0) && (ind<=Ndet)
                                I[m, n] += A[ind]
                            end
                        elseif interpolation == "linear"
                            ind = floor(Int64, rn)
                            if (ind>0) && (ind<Ndet)
                                I[m, n] += (rn - ind)*(A[ind+1]-A[ind, t])  + A[ind, t]
                            end
                        end
                    end
                end
            end
        end
    end

    return I/pi/ (2 * nAngles)
end

function filtered_back_projection_cuda()
    # to be done

end