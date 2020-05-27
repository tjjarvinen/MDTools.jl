abstract type AbstractPotential end

struct LennardJones <: AbstractPotential
    C12::Float64
    C6::Float64
end

function LennardJones(ϵ::Unitful.Energy, σ::Unitful.Length)
    s = auconvert(σ)
    e = auconvert(ϵ)
    A = 4*e*s^12
    B = 4*e*s^6
    return LennardJones(austrip(A),austrip(B))
end

function LennardJones(ϵ::Unitful.Temperature, σ)
    return LennardJones(ϵ*u"k", σ)
end


#(p::LennardJones)(r) = p.C12*r^-12 - p.C6*r^-6

#(p::LennardJones)(r::AbstractArray) = sum(p.C12.*r.^-12 .- p.C6.*r.^-6)

function (p::LennardJones)(r)
    @avx r6 = r.^-6
    return sum(muladd.(p.C12, r6, -p.C6) .* r6)
end
