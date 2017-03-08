Utility for porcupine plots
---------

NOTE: MIToS no está actualizado p/ Julia v0.5, así q hay q usar el v0.47.
---
Apenas se actualice, aviso.
---
Hay q tener Julia 0.5 instalada. Estos repos suelen tener Julia actualizado.
Seguir las instrucciones de la página p/ poder agregarlos e instalar Julia:
`https://launchpad.net/~staticfloat/+archive/ubuntu/juliareleases`

Luego, hay q instalar los siguientes paqutes en Julia: DataFrames, MIToS.PDB, Distributions y ArgParse.

P/ hacer eso:

```
julia> Pkg.add("DataFrames")
julia> Pkg.add("MIToS")
julia> Pkg.add("Distributions")
julia> Pkg.add("ArgParse")
```

Usage
---
`./porky.jl -p INPDB -v VECTOR -m MULTIPLIER -o OUTPDB [-i INDEX] [--script]`

`porky.jl` lee modos de Calpha o modos all atom y desplaza la estructura original
a lo largo del modo. El archivo `VECTOR` puede ser:

- un archivo de modos de PCA salido del `cpptraj` de AMBER. En tal caso se debe incluir el argumento `INDEX` p/ especificar el índice del modo de interés
- un archivo de texto con el modo en 1 sola columna.


Recomiendo usar el flag `--script` p/ q `porky.jl` escriba un script p/ Pymol llamado `script_porky.py`.
De modo tal q luego de correr `porky.jl`, uno haga:
```
pymol script_porky.py
```
y obtenga el porcupine plot automáticamente.

---

modevectors.py by Sean Law & Srinivasa
