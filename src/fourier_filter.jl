filter_name = ("ramp", "shepp-logan", "cosine", "hamming", "hann")

function ramp(N::Int64)
    result = zeros(N)
    if iseven(N)
        center = round(Int64, N/2)
    else 
        center = round(Int64, (N+1)/2)
    end
    result[center]=0.25
    result[center+1:2:N] = (-1.0/pi/pi) ./( center .- collect(center+1:2:N)) .^2
    result[center-1:-2:1] = (-1.0/pi/pi) ./( center .- collect(center-1:-2:1)) .^2
    
    return real(fftshift(fft(result)))
end

function hamming(N::Int64)
    a = 0.54
    return a .- (1.0-a) .* cos.((2*pi/(N-1)) .* (Array(0:1:N-1)))
end

function hanning(N::Int64)
    a = 0.5
    return a .- (1.0-a) .* cos.((2*pi/(N-1)) .* (Array(0:1:N-1)))
end

function fourier_filter(N::Int, filtername = "ramp")
    filtername = lowercase(filtername)
    
    if filtername ∉ filter_name 
        return nothing
    end
    
    filter::Array{Float32, 1} = 1.0 .- abs.(range(-1, stop=1, length=N))
    
    if filtername == "ramp"
    elseif filtername == "shepp-logan"
        omega::Array{Float32, 1} = Float32.(pi .* fftfreq(N)[2:N])
        filter[2:N] .*= sin.(omega) ./ omega
        
    elseif filtername == "cosine"
        freq::Array{Float32, 1} = range(0, stop = np.pi, length = N+1)[1:N]
        cosine_filter::Array{Float32, 1} = fftshift(sin.(freq))
        filter .*= cosine_filter
    elseif filtername == "hamming"
        filter .*= fftshift(hamming(N))
    elseif filtername == "hann"
        filter .*= fftshift(hanning(N))
    end
    
    return filter
end

function fourier_filter2(N::Int, filtername = "ramp")
    filtername = lowercase(filtername)
    
    if filtername ∉ filter_name 
        return nothing
    end
    
    n = [collect(1:2:(N/2)); collect((N/2-1):-2:1)]
    f = zeros(Float32, N)
    f[1]=0.25
    @inbounds for j in 2:2:N
        f[j] = -1/(π*n[j÷2])^2
    end

    filter = 2*real(fft(f))

    if filtername == "ramp"
    elseif filtername == "shepp-logan"
        omega::Array{Float32, 1} = Float32.(pi .* fftfreq(N)[2:N])
        filter[2:N] .*= sin.(omega) ./ omega
        
    elseif filtername == "cosine"
        freq::Array{Float32, 1} = range(0, stop = np.pi, length = N+1)[1:N]
        cosine_filter::Array{Float32, 1} = fftshift(sin.(freq))
        filter .*= cosine_filter
    elseif filtername == "hamming"
        filter .*= fftshift(hamming(N))
    elseif filtername == "hann"
        filter .*= fftshift(hanning(N))
    end
    
    return filter
end