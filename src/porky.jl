#!/home/german/julia-ae26b25d43/bin/julia
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
using DataFrames
using MIToS.PDB
using Distributions
##########
# functions
##########
function read_ptraj_modes(file, modes_elements, norma::Bool=true)
    modes_file=open(file, "r")
    modes_text = readdlm(modes_file, skipstart=0, skipblanks=true,
    ignore_invalid_chars=true, comments=true, comment_char='\*')
    close(modes_file)

    nmodes = modes_text[1, 5]
    ncoords = convert(Int64, modes_elements)
    lines = ceil(Int64, ncoords/7)
    rest = convert(Int64, ncoords % 7)

    eval=Array{Float64}(nmodes);
    mode = Array{Float64}(ncoords, nmodes);
    temp1=Array{Float64}(ncoords, 1);
    temp2 = Array{Float64}(ncoords+(7-rest));

    j=lines + 1 + 2 # 1 p/ q lea la prox linea 2 por el header

    for i=1:nmodes
        eval[i] = modes_text[j, 2]
        temp = transpose(modes_text[(j+1):(lines+j), :])
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
function displaceAA(mod_pdb, vector1, multiplier)
  # Preparo variables
    pdb = copy(mod_pdb)
    struct_xyz = coordinatesmatrix(pdb)
    new_struct_xyz = copy(struct_xyz)
    natom = Array{Int64}(1)
    vector = Array{Float64}(1, 3)
   aa = length(pdb)
   # Determino el nro de atomos de c/ aminoácido
   for i=1:aa
       push!(natom, length(pdb[i]))
   end
   shift!(natom)
   temp1 = Array{Int64}(natom[1],3)

   # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    for i=1:3:length(vector1)
        if i== 1
            vector = reshape(vector1[i:i+2], 1, 3)
            continue
        end
        vector = vcat(vector, reshape(vector1[i:i+2], 1, 3))
    end

   for i=1:aa
       if i == 1
           temp1 = repmat(vector[i, :], natom[i], 1)
           continue
       end
       temp2 = repmat(vector[i, :], natom[i], 1)
       temp1 = vcat(temp1, temp2)
   end
   sum_mat = temp1

   # Listo, ahora puedo mover el pdb
   new_struct_xyz  = struct_xyz + sum_mat .* multiplier
   pdb = change_coordinates(pdb, new_struct_xyz);
   return pdb
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

##########
# main program
##########

# Read arguments from console
message = string("\n\nUsage:\n",  "./porky.jl <input PDB> ",
"<input vector or AMBER PCA> <multiplier> <output PDB> \"AMBER mode number\"", "\n\n")
message_vec_not_found = string(" vector file not found \n\n")

# Check for bad input
if length(ARGS) < 3 || length(ARGS) > 5
    throw(ArgumentError(message))
end

# Get ready
in_vec = Array{Float64, 1}
in_pdb_filename = ARGS[1]
in_vec_filename = ARGS[2]
multiplier = ARGS[3]
multiplier = parse(Int64, multiplier)
out_pdb_filename = ARGS[4]

# Read PDB
in_pdb = read(string(in_pdb_filename), PDBFile, group="ATOM");
nres_xyz = 3*length(in_pdb)
natom_xyz = size(coordinatesmatrix(in_pdb))[1] * 3

if length(ARGS) == 5
# Vector de PCA Amber
    amber_idx = ARGS[5]
    amber_idx = parse(Int64, amber_idx)
    try
        in_vec = read_ptraj_modes(in_vec_filename, nres_xyz, true)[1][:, amber_idx]
    catch
        try
            in_vec = read_ptraj_modes(in_vec_filename, natom_xyz, true)[1][:, amber_idx]
        end
    end
else
# Vector puro
    in_vec = convert(Array{Float64}, readtable(in_vec_filename)[:, 1]);
end

# In case input vector file is note found
if in_vec == Array{Float64, 1}
    throw(ArgumentError(string("\n\n", in_vec_filename, message_vec_not_found)))
end

# Ahora desplazo
if nres_xyz == length(in_vec)
# El modo es de Calpha
    out_pdb = displaceAA(in_pdb, in_vec, multiplier);

elseif natom_xyz == length(test_vec)
# El modo es all-atom
    out_pdb = displaceAtoms(in_pdb, in_vec, multiplier);

else
# El modo no tiene el tamano adecuado
    println("PDB and input vector don't match.")
    println("PDB has ", length(in_pdb) , " amino acids and ", size(coordinatesmatrix(in_pdb))[1] * 3, " atoms.")
    println("Vector has ", length(in_vec), " elements, which should correspond to ", length(in_vec) / 3, " particles.")
end

# Y guardo
write(out_pdb_filename, out_pdb, PDBFile)

# Finalmente, hago el script
load = "cmd.load(\""
f = open("script_porky.py", "w")
write(f, "from pymol.cgo import *\n")
write(f, "from pymol import cmd\n\n")
write(f, load, in_pdb_filename,"\")\n")
write(f, load, out_pdb_filename,"\")\n")
write(f, load,"modevectors.py\")\n")
write(f, "modevectors(\"", in_pdb_filename[1:end-4], "\", \"", out_pdb_filename[1:end-4], "\", ")
write(f, "outname=\"modevectors\", head=1.0, tail=0.3, headrgb = \"1.0, 1.0, 0.0\", tailrgb = \"1.0, 1.0, 0.0\") ")
close(f)
