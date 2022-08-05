module ExperimentsStandardMap

export F, orbit

F(u::Vector; L = 1.0) = [mod(u[1]+(L/2π)*sin(2π*u[2]), 1); mod(u[2]+u[1]+(L/2π)*sin(2π*u[2]), 1)] 

F(u::Matrix; L = 1.0) = vcat(mod.(u[1,:]+(L/2π)*sin.(2π*u[2,:]), 1)', mod.(u[2, :]+u[1, :]+(L/2π)*sin.(2π*u[2, :]), 1)')

function orbit(v::Matrix{T}, n, L) where {T}
    i, j = size(v)
    out_orbit = zeros(T, (2*n, j))
    out_orbit[1:2, :] = v
    for i in 2:n
        v = F(v; L = L)
        out_orbit[2*i-1:2*i, :] = v
    end
    return out_orbit
end

end