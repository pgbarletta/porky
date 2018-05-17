#!/usr/bin/env julia
###############################################################################
###############   utility to displace a PDB along a vector   ##################
#   code by pgbarletta
#########################################################
# Safety pig included:
#
#    _._ _..._ .-',     _.._(`))
#   '-. `     '  /-._.-'    ',/
#      )         \            '.
#     / _    _    |             \
#    |  a    a    /              |
#    \   .-.                     ;
#     '-('' ).-'       ,'       ;
#        '-;           |      .'
#           \           \    /
#           | 7  .__  _.-\   \
#           | |  |  ``/  /`  /
#          /,_|  |   /,_/   /
#             /,_/      '`-'
###############################################################################
using Chemfiles
using StaticArrays
using ArgParse
##########
# functions
##########
function read_ptraj_modes(filename, nmodes::Int64=0, norma::Bool=true)
    modes_text = readdlm(filename, skipstart=0, skipblanks=true, comments=true,
        comment_char='\*')

    if nmodes == 0
        nmodes = modes_text[1, 5]
    end
    modes_elements = modes_text[2, 1]

    ncoords = convert(Int64, modes_elements)
    lines = ceil(Int64, ncoords/7)
    rest = convert(Int64, ncoords % 7)

    eval = Array{Float64}(nmodes);
    mode = Array{Float64}(ncoords, nmodes);
    temp1 = Array{Float64}(ncoords, 1);
    temp2 = Array{Float64}(ncoords+(7-rest));

    j=lines + 1 + 2 # 1 p/ q lea la prox linea 2 por el header
    for i=1:nmodes
        eval[i] = modes_text[j, 2]
        temp = permutedims(modes_text[(j+1):(lines+j), :], [2, 1])
        temp2 = reshape(temp, ncoords+(7-rest))
        for k=(rest+1):7
            pop!(temp2)
        end
        mode[:, i] = temp2
        j = j + lines + 1
    end

    if norma == true
        for i=1:nmodes
            mode[: ,i] = mode[:, i] / norm(mode[:, i])
        end
    end

    return mode, eval
end
#########
function displaceAA(in_frm, aa, aa_3, in_vec)
    # Preparo variables
    const in_top = Topology(in_frm)
    natoms = convert(Int64, size(in_top))
    const in_xyz = positions(in_frm)

    # Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(aa)
    ids = Array{Int64}(aa)
    [ ids[i+1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa-1 ]
    const idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido
    [ tmp[i+1] = size(Residue(in_top, i)) for i = 0:aa-1 ]
    const natom_aa = tmp[idx]

    # Paso el vector columna de tamaño 1xaa_3 a 3xaa
    const vector = reshape(in_vec, 3, aa)
    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    sum_mat = Array{Float64}(3, natoms)
    cursor = 0
    for i = 1:aa
        if i == 1
            sum_mat[:, 1:natom_aa[i]] = repmat(vector[:, 1], 1, natom_aa[i])
            cursor = natom_aa[i]
            continue
        end
        const rango = collect(cursor+1:cursor + natom_aa[i])
        sum_mat[:, rango] = repmat(vector[:, i], 1, natom_aa[i])
        cursor += natom_aa[i]
    end

    # Listo, ahora puedo mover el pdb
    out_frm = deepcopy(in_frm)
    out_xyz = positions(out_frm)

    # Tengo q hacer esto por ahora. Hasta q arreglemos Chemfiles.
    for i = 1:size(in_xyz)[1]
        for j = 1:size(in_xyz)[2]
            out_xyz[i, j]  = round(in_xyz[i, j] + sum_mat[i, j], 3)
        end
    end
    return out_frm
end
#########
function displaceAtoms(mod_pdb, vector1, multiplier)
  # Preparo variables
    pdb = copy(mod_pdb)
    struct_xyz = coordinatesmatrix(pdb)
#    new_struct_xyz = copy(struct_xyz)
    vector = Array{Float64}(1, 3)

    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    for i=1:3:length(vector1)
        if i== 1
            vector = reshape(vector1[i:i+2], 1, 3)
            continue
        end
        vector = vcat(vector, reshape(vector1[i:i+2], 1, 3))
    end

    # Listo, ahora puedo mover el pdb
    new_struct_xyz  = struct_xyz + vector .* multiplier
    pdb = change_coordinates(pdb, new_struct_xyz);
   return pdb
end
#########
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
        arg_type = Int
        default = 1
    "--index", "-i"
        help = "Mode number when reading cpptraj PCA modes"
        arg_type = Int
        default = 0
    "--script"
        help = "Only write output PDB"
        action = :store_true
end

##########
# main program
##########

# Read arguments from console
parsed_args = parse_args(ARGS, s)
args = Array{Any, 1}(0)
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
in_vec = Array{Float64, 1}
# Read PDB
const in_trj = Trajectory(inpdb)
const in_frm = read(in_trj)
const in_top = Topology(in_frm)
const aa = convert(Int64, count_residues(in_top))
const aa3 = aa * 3


in_vec = MVector{aa3, Float64}()
if (index == 0)
    try
        in_vec = convert(MVector{aa3, Float64}, readdlm(vector))
    catch e
        println(vector, " error:")
        error(e)
    end
else
    try
        in_modes, in_evals = read_ptraj_modes(vector, index)
        in_vec = convert(MVector{aa3, Float64}, in_modes[:, index])
    catch e
        println(vector, " error:")
        error(e)
    end
end

#
in_vec = in_vec / norm(in_vec) .* multiplier
const out_frm = displaceAA(in_frm, aa, aa3, in_vec);
# Y guardo
out_trj = Trajectory(outpdb, 'w')
write(out_trj, out_frm)

const script_filename = string("script_porky_", splitext(outpdb)[1] ,".py")
# Finalmente, hago el script
if script == true
    load = "cmd.load(\""
    f = open(script_filename, "w")
    write(f, "from pymol.cgo import *\n")
    write(f, "from pymol import cmd\n\n")
    write(f, load, inpdb,"\")\n")
    write(f, load, outpdb,"\")\n")
    write(f, load,"modevectors.py\")\n")
    write(f, "modevectors(\"", inpdb[1:end-4], "\", \"", outpdb[1:end-4], "\", ")
    write(f, "outname=\"modevectors\", head=1.0, tail=0.3, headrgb = \"1.0, 1.0, 0.0\", tailrgb = \"1.0, 1.0, 0.0\") ")
    close(f)
end
