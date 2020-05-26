
abstract type AbstractMDIntegrator end


struct LeapFrog <: AbstractMDIntegrator
    Δt::Float64
end

struct VelocityVerlet <: AbstractMDIntegrator
    Δt::Float64
end


function propagate!(coordinates, velocities, potential, mass, integrator::LeapFrog)
    velocities .+= force(potential, coordinates).*integrator.Δt
    coordinates .+= velocities.*integrator.Δt
    return  coordinates, velocities
end

function propagate!(coordinates, velocites, potential, integrator::VelocityVerlet)
    body
end
