module ExperimentsStandardMap

include("FromDS.jl")

export F, orbit

F(u::Vector; L = 2π) = [mod(u[1]+(L/2π)*sin(2π*u[2]), 1); mod(u[2]+u[1]+(L/2π)*sin(2π*u[2]), 1)] 

F(u::Matrix; L = 2π) = vcat(mod.(u[1,:]+(L/2π)*sin.(2π*u[2,:]), 1)', mod.(u[2, :]+u[1, :]+(L/2π)*sin.(2π*u[2, :]), 1)')

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

"""
    initial_cond_prop_v(k, v)

This function extracts `k` random numbers in [0,1] and returns
a matrix whose columns are given by the product of these random numbers with v
"""
function initial_cond_prop_v(k, v, ϵ = 0.1)
    x0 = zeros(2, k)
    for (i, λ) in enumerate(range(0, ϵ, length = k))
        x0[:, i] = λ*v
    end
    return x0 
end

using Interpolations

function compute_unstable(L, n, k, v, ϵ)
    x = initial_cond_prop_v(k, v, ϵ)
    @info length(x)
    for _ in 1:n
        x = F(x, L = L)
    end

    t = range(0, ϵ, length = k)
    itp = Interpolations.scale(interpolate(x', (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)
    return itp
end

DF(u::Vector; L = 2π) = [1.0 2π*(L/2π)*cos(2π*u[2]); 1.0 1.0+2π*(L/2π)*cos(2π*u[2])]

# ugly, hacky stuff!!!
function invF(v; L = 2π)
    x = zeros(2)
    for _ in 1:10
        #@info "x" x
        #@info DF(x, L = L)
        #@info F(x, L= L)-v
        x = x - DF(x, L = L)\(F(x, L = L)-v) 
        x = mod1.(x, 1)
    end
    return x
end

function compute_stable(L, n, k, v, ϵ)
    x = initial_cond_prop_v(k, v, ϵ)
    @info length(x)
    for _ in 1:n
        x = F(x, L = L)
    end

    t = range(0, ϵ, length = k)
    itp = Interpolations.scale(interpolate(x', (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)
    return itp
end






end