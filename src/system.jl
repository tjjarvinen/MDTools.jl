
abstract type AbstractBox end

mutable struct CubicBox <: AbstractBox
    a::Float64
end

struct NonPeriodixBox <: AbstractBox
    nothing
end


CubicBox(a::Unitful.Length) = CubicBox(austrip(a))






mutable struct MoleculeSystem
    atoms::Vector{String}
    mass::Vector{Float64}
    function MoleculeSystem()
        new(Vector{String}(), Vector{Float64}())
    end
    function MoleculeSystem(atoms::AbstractVector{String})
        new(atoms, [MDTools.atom_masses[name] for name in atoms])
    end
end



mutable struct SimulationSystem
    mol::MoleculeSystem
    simbox::CubicBox
    integrator::LeapFrog
    coordinates::Matrix{Float64}
    velocities::Matrix{Float64}
    potential
    function SimulationSystem()
        new(MoleculeSystem(), CubicBox(Inf), LeapFrog(0.5u"fs"),
            Matrix{Float64}(undef,0,0), Matrix{Float64}(undef,0,0), Nothing
        )
    end
end

Base.length(m::MoleculeSystem) = length(m.atoms)

function add!(m::MoleculeSystem, name::AbstractString)
    push!(m.atoms, name)
    push!(m.mass, MDTools.atom_masses[name])
    return m
end


function (*)(n::Integer, m::MoleculeSystem)
    return MoleculeSystem(repeat(m.atoms, n))
end

(*)(m::MoleculeSystem, n::Integer) = n*m

function (+)(m1::MoleculeSystem, m2::MoleculeSystem)
    return MoleculeSystem(vcat(m1.atoms, m2.atoms))
end

function add!(ss::SimulationSystem, m::MoleculeSystem; mindis=3.0u"Å", maxtries=100)
    L = PeriodicEuclidean([ss.simbox.a, ss.simbox.a, ss.simbox.a])
    for i in 1:length(m)
        if size(ss.coordinates,1) == 3
            for j in 1:maxtries
                r = ss.simbox.a .* rand(3,1)
                d=pairwise(L, ss.coordinates, r; dims=2)
                if ! any(d .< austrip(mindis))
                    ss.coordinates = hcat(ss.coordinates, r)
                    ss.mol += MoleculeSystem([m.atoms[i]])
                    break
                end
                j == maxtries && error("Out of tries")
            end

        else
            ss.coordinates = ss.simbox.a .* rand(3,1)
            ss.mol += MoleculeSystem([m.atoms[i]])
        end
    end
    return ss
end
