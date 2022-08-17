using DynamicalSystems

ds = Systems.standardmap(; k = 1.2)
xs = range(0, stop = 2π, length = 11); ys = copy(xs)
ics = [SVector{2}(x,y) for x in xs for y in ys]

# All permutations of [±1, ±1]:
singss = lambdaperms(2)[2] # second entry are the signs

# I know from personal research I only need this `inds`:
indss = [[1,2]] # <- must be container of vectors!

λs = 0.005 # <- only this allowed to not be vector (could also be vector)

orders = [3, 4]
ALLFP = Dataset{2, Float64}[];

for o in orders
    FP = periodicorbits(ds, o, ics, λs, indss, singss)
    push!(ALLFP, FP)
end

using StaticArrays

function FSVector( (θ, p); k = 1.0 ) 
    p1 = p + k*sin(θ)
    θ1 = θ + p1

    while p1 >= 2π 
        p1-=2π
    end
    while p1 < 0 
        p1+=2π
    end

    while θ1 >= 2π 
        θ1 -=2π
    end
    while θ1 < 0 
        θ1 +=2π
    end
    return SVector(θ1, p1)
end

DFSVector( (θ, p); k = 1.0 ) = SMatrix{2,2}([1.0+k*cos(θ) 1.0; k*cos(θ) 1.0])

function Period3Root((u1, u2, u3, u4, u5, u6); k = 1.0)
    v1, v2 = FSVector((u1, u2); k = k)-SVector(u3, u4)
    v3, v4 = FSVector((u3, u4); k = k)-SVector(u5, u6)
    v5, v6 = FSVector((u5, u6); k = k)-SVector(u1, u2)
    return SVector(v1, v2, v3, v4, v5, v6)
end

function DPeriod3Root((u1, u2, u3, u4, u5, u6); k = 1.0)
    T = typeof(u1)
    M = zeros(T, (6, 6))
    M[1:2, 1:2] = DFSVector((u1, u2); k = k)
    M[1:2, 3:4] = [-1 0; 0 -1]
    M[3:4, 3:4] = DFSVector((u3, u4); k = k)
    M[3:4, 5:6] = [-1 0; 0 -1]
    M[5:6, 5:6] = DFSVector((u5, u6); k = k)
    M[5:6, 1:2] = [-1 0; 0 -1]
    return SMatrix{6,6}(M)
end

function Period4Root((u1, u2, u3, u4, u5, u6, u7, u8); k = 1.0)
    v1, v2 = FSVector((u1, u2); k = k)-SVector(u3, u4)
    v3, v4 = FSVector((u3, u4); k = k)-SVector(u5, u6)
    v5, v6 = FSVector((u5, u6); k = k)-SVector(u7, u8)
    v7, v8 = FSVector((u7, u8); k = k)-SVector(u1, u2)
    return SVector(v1, v2, v3, v4, v5, v6, v7, v8)
end

function DPeriod4Root((u1, u2, u3, u4, u5, u6, u7, u8); k = 1.0)
    T = typeof(u1)
    M = zeros(T, (8, 8))
    M[1:2, 1:2] = DFSVector((u1, u2); k = k)
    M[1:2, 3:4] = [-1 0; 0 -1]
    M[3:4, 3:4] = DFSVector((u3, u4); k = k)
    M[3:4, 5:6] = [-1 0; 0 -1]
    M[5:6, 5:6] = DFSVector((u5, u6); k = k)
    M[5:6, 7:8] = [-1 0; 0 -1]
    M[7:8, 7:8] = DFSVector((u7, u8); k = k)
    M[7:8, 1:2] = [-1 0; 0 -1]
    return SMatrix{8,8}(M)
end

using IntervalArithmetic

function buildbox3(z1, z2, z3; ϵ = 0.01)
    return Interval(z1[1]-ϵ, z1[1]+ϵ) × Interval(z1[2]-ϵ, z1[2]+ϵ) × 
            Interval(z2[1]-ϵ, z2[1]+ϵ) × Interval(z2[2]-ϵ, z2[2]+ϵ) ×
            Interval(z3[1]-ϵ, z3[1]+ϵ) × Interval(z3[2]-ϵ, z3[2]+ϵ)
end

function buildbox4(z1, z2, z3, z4; ϵ = 0.01)
    return Interval(z1[1]-ϵ, z1[1]+ϵ) × Interval(z1[2]-ϵ, z1[2]+ϵ) × 
            Interval(z2[1]-ϵ, z2[1]+ϵ) × Interval(z2[2]-ϵ, z2[2]+ϵ) ×
            Interval(z3[1]-ϵ, z3[1]+ϵ) × Interval(z3[2]-ϵ, z3[2]+ϵ) ×
            Interval(z4[1]-ϵ, z4[1]+ϵ) × Interval(z4[2]-ϵ, z4[2]+ϵ)
end

function NewtonStep(f, df, x::IntervalBox{N, T}) where {N, T}
    m = IntervalBox{N, T}(mid.(x))
    return m-IntervalBox{N, T}(df(x)\f(m))
end

function NewtonStep(f, df, x::Array{Interval{T}, 1}) where {N, T}
    m = Interval{T}.(mid.(x))
    DF = Matrix(df(x))
    return m-Interval{T}.(DF\f(m))
end

