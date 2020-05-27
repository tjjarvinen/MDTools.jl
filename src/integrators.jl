
abstract type AbstractMDIntegrator end


struct LeapFrog <: AbstractMDIntegrator
    Δt::Float64
end

mutable struct VelocityVerlet <: AbstractMDIntegrator
    Δt::Float64
end


LeapFrog(Δt::Unitful.Time) = LeapFrog(austrip(Δt))


function propagate!(coordinates, velocities, potential, mass, integrator::LeapFrog)
    velocities .+= force(potential, coordinates).*integrator.Δt ./mass
    coordinates .+= velocities.*integrator.Δt
    return  coordinates, velocities
end

function propagate!(r1, v1, f1, r0, v0, f0, potential, mass, integrator::VelocityVerlet)
    r1 .= r0 .+ v0.*integrator.Δt .+ integrator.Δt^2*f1./mass
    v1 .= v0 .+ integrator.Δt.*(f0 .+ force!(f1,r1))./(2*mass)
    return r1, v1
end
