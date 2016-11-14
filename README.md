Utility for porcupine plots
---------

NOTE: MIToS no está actualizado p/ Julia v0.5, así q hay q usar el v0.47.
---
Apenas se actualice, aviso.
---
Hay q tener Julia 0.5 instalada. Estos repos suelen tener Julia actualizado.
Seguir las instrucciones de la página p/ poder agregarlos e instalar Julia:
`https://launchpad.net/~staticfloat/+archive/ubuntu/juliareleases`

Luego, hay q instalar los siguientes paqutes en Julia: DataFrames, MIToS.PDB y Distributions
P/ hacer eso:

julia> Pkg.add("DataFrames")
julia> Pkg.add("MIToS")
julia> Pkg.add("Distributions")


Usage:
---
`./porky.jl <input PDB> <input vector> <multiplier> <output PDB> "AMBER mode number"`

porky.jl lee modos de Calpha o modos all atom y desplaza la estructura original
a lo largo del modo. El 5to argumento es opcional y solo debe ser ingresado
cuando `<input vector>` sea un archivo de modos de PCA salido del `cpptraj` de
AMBER.
