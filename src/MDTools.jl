module MDTools

using Zygote
using Unitful
using UnitfulAtomic
using Distances
using LoopVectorization

import Base.*
import Base.+

include("integrators.jl")
include("atoms.jl")
include("potentials.jl")
include("system.jl")



export AbstractMDIntegrator,
       add!,
       CubicBox,
       energy,
       force,
       LeapFrog,
       propagate!


function force!(f, potential, coordinates)
    f .= - gradient(potential,coordinates)[1]
end

force(potential, coordinates) = - gradient(potential, coordinates)[1]


function energy(coordinates, velocity, potential, mass)
    Ek = 0.5*sum( mass.*velocity.^2 )
    Ep = potential(coordinates)
    return Ek, Ep, Ek+Ep
end

greet() = print("Hello World!")

end # module
