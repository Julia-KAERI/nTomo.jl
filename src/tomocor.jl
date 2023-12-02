struct TomoCor
    img
    c1
    c2
    x
    y
end


function calc_cor2(tomo::TomoData, y1=nothing, y2=nothing)

    M, N = size(tomo.data)[2:3]

    c1, c2 = ceil(Int64, N/3), floor(Int64, N*2/3)
    
    if y1 === nothing
        y1 = 1
    end

    if y2 === nothing
        y2 = M
    elseif y2 < 0
        y2 = M-y2
    end

    y1, y2 = minmax(y1, y2)
    yi, yf = ceil(Int64, y1), floor(Int64, y2)
    @assert 1 ≤ yi < yf ≤ M

    p = tomo.data[1, yi:yf, :] + tomo.data[end, yi:yf, :]
    cor = []
    
    for y in yi:1:yf
        tt = []      
        for c in c1:1:c2
            ts = 0
            if (c <= (c1+c2)*0.5)
                i0, i1 = 1, 2*c-1
            else 
                i0, i1 = 2*c-N, N  
            end

            for x in (i0+1):1:(i1-1)
                ts += (p[y-yi+1, x]-p[y-yi+1, i1-x+i0])^2
            end
            push!(tt, ts) 
        end
        push!(cor, [y, findmin(tt)[2]+c1-1])  
    end
    
    xyarr = reduce(hcat, cor)'
    y, x = xyarr[:,1], xyarr[:,2]
    
    A = [sum(y.^2) sum(y) ; sum(y) length(y)]
    B = [sum(y.*x) ; sum(x)]
    C = inv(A)*B
    retimg = tomo.data[1, :, :] + tomo.data[end, :, :]
    return TomoCor(retimg, C[1], C[2], x, y)
end

function (f::TomoCor)(y::Integer)
    return f.c1*y+f.c2 
end
