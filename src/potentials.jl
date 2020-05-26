abstract type AbstractPotential end

struct LennardJones <: AbstractPotential
    C6::Float64
    C12::Float64
end

(p::LennardJones)(r) = p.C12*r^-12 - p.C6*r^-6

(p::LennardJones)(r::AbstractArray) = sum(p.C12.*r.^-12 .- p.C6.*r.^-6)
