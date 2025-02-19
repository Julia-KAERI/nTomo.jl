function filtered_back_projection(
    tomo::TomoData, 
    cor::TomoCor,
    index = 1,
    filtername = "hann", 
    interpolation = "nearest", 
    with_threads = true)
    
    nAngles, Nheight, Ndet = size(tomo.data)

    @assert 1 ≤ index ≤ Nheight

    return filtered_back_projection(tomo.data[:, index, :], cor(index), tomo.ths, filtername, interpolation, with_threads)
end


"""
`recon_fbp(tomo::TomoData, filtername = "hann", interpolation = "nearest", with_threads = true)`

Reconstruct image by filtered back projection method.

parameters
==========
tomo : tomodata
filtername : one of ("ramp", "shepp-logan", "cosine", "hamming", "hann")
interpolation : one of ("nearst", "linear")
"""
function recon_fbp(
    tomo::TomoData,
    filtername = "hann",
    interpolation = "nearest",
    with_threads = true
    )

    nAngles, Nheight, Ndet = size(tomo.data)

    img_reconed = Array{Float32, 3}(undef, (Ndet, Ndet, Nheight))

    Threads.@threads for i in 1:Nheight
        # @inbounds img_reconed[:, :, i] = filtered_back_projection(tomo.data[:, i, :], tomo.ths, tomo_cor(tomo, i), filtername) 
        @inbounds img_reconed[:, :, i] = iradon_fbp(tomo.data[:, i, :], tomo.ths, tomo_cor(tomo, i), filtername) 
    end
    
    return img_reconed
end

function recon_fbp!(
    tomo::TomoData,
    filtername = "hann",
    interpolation = "nearest",
    with_threads = true
    )

    nAngles, Nheight, Ndet = size(tomo.data)

    img_reconed = Array{Float32, 3}(undef, (Ndet, Ndet, Nheight))

    Threads.@threads for i in 1:Nheight
        @inbounds img_reconed[:, :, i] = iradon_fbp(tomo.data[:, i, :], tomo.ths, tomo_cor(tomo, i), filtername) 
    end
    
    push!(tomo.process, :recon_fbp)

    tomo.data = img_reconed
end

# function recon_fbp_parts(
#     tomo::TomoData,
#     cor::TomoCor,
#     filtername = "hann",
#     interpolation = "nearest",
#     with_threads = true,
#     nparts::Integer = 2,
#     workingdir = ""
#     )

#     @assert nparts ≥ 2 
#     nAngles, Nheight, Ndet = size(tomo.data)

#     starts_at = [(i*(Nheight ÷ nparts) +1) for i in 0:(nparts-1)]
#     partsinfo = []

#     for j in 1:nparts
#         fn = "recon" * "_" * lpad(j, 4, "0") * ".bin"
#         fpath = joinpath(workingdir, fn)
#         if j ≠ nparts
#             pheight = starts_at[j+1]-starts_at[h]
#             img_reconed = Array{Float32, 3}(undef, (Ndet, Ndet, pheight))
#         else 
#             pheight = Nheights - starts_at[nparts] +1
#         end
        
#         Threads.@threads for i in 1:pheight
#             img_reconed[:,:, i] = filtered_back_projection(tomo.data[:, i, :], cor(i), tomo.ths, filtername, interpolation, with_threads) 
#         end
        
#         write(fpath, fbpimg)
#         push!(partsinfo, (Ndet, Ndet, pheight, starts_at[j], fpath))
#     end
    
#     return partsinfo
#     #push!(tomo.process, :recon_fbp)

#     #tomo.data = img_reconed
# end
