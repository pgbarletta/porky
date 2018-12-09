#!/usr/bin/env julia
###############################################################################
###############   utility to displace a PDB along a vector   ##################
#   code by pgbarletta
#########################################################
using Chemfiles, JUMD
using StaticArrays
using ArgParse, LinearAlgebra, DelimitedFiles

function displaceAA(in_frm, aa, aa_3, in_vec)
    # Preparo variables
    in_top = Topology(in_frm)
    natoms = convert(Int64, size(in_top))
    in_xyz = positions(in_frm)

    # Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(undef, aa)
    ids = Array{Int64}(undef, aa)
    [ ids[i + 1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa - 1 ]
    idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido
    [ tmp[i + 1] = size(Residue(in_top, i)) for i = 0:aa - 1 ]
    natom_aa = tmp[idx]

    # Paso el vector columna de tamaño 1xaa_3 a 3xaa
    vector = reshape(in_vec, 3, aa)
    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    sum_mat = Array{Float64}(undef, 3, natoms)
    cursor = 0
    for i = 1:aa
        if i == 1
            sum_mat[:, 1:natom_aa[i]] = repeat(vector[:, 1], 1, natom_aa[i])
            cursor = natom_aa[i]
            continue
        end
        rango = collect(cursor + 1:cursor + natom_aa[i])
        sum_mat[:, rango] = repeat(vector[:, i], 1, natom_aa[i])
        cursor += natom_aa[i]
    end

    # Listo, ahora puedo mover el pdb
    out_frm = deepcopy(in_frm)
    out_xyz = positions(out_frm)

    # Tengo q hacer esto por ahora. Hasta q arreglemos Chemfiles.
    for i = 1:size(in_xyz)[1]
        for j = 1:size(in_xyz)[2]
            out_xyz[i, j]  = round(in_xyz[i, j] + sum_mat[i, j], digits = 3)
        end
    end
    return out_frm
end
# Arg Parse settings
s = ArgParseSettings()
@add_arg_table s begin
    "--inpdb", "-p"
    help = "Input PDB"
    arg_type = String
    required = true
    "--vector", "-v"
    help = "Input vector"
    arg_type = String
    required = true
    "--outpdb", "-o"
    help = "Output PDB"
    arg_type = String
    required = true
    "--multiplier", "-m"
    help = "Constant to scale the normalized vector. Defaults to 1"
    arg_type = Float64
    default = 1.
    "--index", "-i"
    help = "Mode number when reading cpptraj PCA modes"
    arg_type = Int
    default = 0
    "--script", "-s"
    help = "Only write output PDB"
    action = :store_true
end

# Read arguments from console
parsed_args = parse_args(ARGS, s)
args = Array{Any,1}(undef, 0)
for (arg, val) in parsed_args
    arg = Symbol(arg)
    @eval (($arg) = ($val))
end
# Append ".pdb" to output pdb
outpdb = outpdb * ".pdb"

println("Input parameters:")
println("INPDB          ", inpdb)
println("VECTOR         ", vector)
println("OUTPDB         ", outpdb)
println("MULTIPLIER     ", multiplier)
println("INDEX          ", index)
println("SCRIPT         ", script)

# Get ready
in_vec = Array{Float64,1}
# Read PDB
const in_trj = Trajectory(inpdb)
const in_frm = read(in_trj)
const in_top = Topology(in_frm)
const aa = convert(Int64, count_residues(in_top))
const aa3 = aa * 3

if (index == 0)
    try
        global in_vec = convert(Array{Float64,1}, readdlm(vector)[:, 1])
    catch e
        println(vector, " error:")
        error(e)
    end
else
    try
        in_modes, in_evals = JUMD.readPtrajModes(vector, index)
        global in_vec = convert(Array{Float64,1}, in_modes[:, index])
    catch e
        println(vector, " error:")
        error(e)
    end
end

in_vec = in_vec ./ norm(in_vec) .* multiplier
const out_frm = displaceAA(in_frm, aa, aa3, in_vec);
# Y guardo
out_trj = Trajectory(outpdb, 'w')
write(out_trj, out_frm)

const script_filename = string("script_porky_", splitext(outpdb)[1], ".py")
# Finalmente, hago el script
if script == true
    load = "cmd.load(\""
    f = open(script_filename, "w")
    write(f, "from pymol.cgo import *\n")
    write(f, "from pymol import cmd\n\n")
    write(f, load, inpdb, "\")\n")
    write(f, load, outpdb, "\")\n")
    write(f, load, "modevectors.py\")\n")
    write(f, "rgb=\"1.0, .5, .5   \"\n")
    write(f, "modevectors(\"", inpdb[1:end - 4], "\", \"", outpdb[1:end - 4], "\", ")
    write(f, "outname=\"", string(splitext(outpdb)[1], "_porky"), "\", head=0.5, tail=0.3, headrgb = rgb, tailrgb = rgb, cutoff=3.0)\n")
    close(f)
end
