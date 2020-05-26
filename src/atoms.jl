
const atom_names = Dict(
    "H"=>1,
    "He"=>2,
    "Li"=>3,
    "Be"=>4,
    "B"=>5,
    "C"=>6,
    "N"=>7,
    "O"=>8,
    "F"=>9,
    "Ar"=>18
)

const m_electron = 9.10938291E-31
const m_u = 1.660539040E-27
const m_au = m_u/m_electron
const proton_mass = 1836.0

const atom_masses = Dict(
    "H"=>1.008*m_au,
    "He"=>4.0026*m_au,
    "Li"=>6.94*m_au,
    "Be"=>9.0122*m_au,
    "B"=>10.81*m_au,
    "C"=>12.011*m_au,
    "N"=>14.007*m_au,
    "O"=>15.999*m_au,
    "F"=>18.998*m_au,
    "Ne"=>20.180*m_au,
    "Ar"=>39.948*m_au
)
