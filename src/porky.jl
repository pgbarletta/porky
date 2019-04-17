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

function write_porcu_script(script::String, inpdb::String, outpdb::String, red::Float64, green::Float64, blue::Float64)
    rgb = string(red, ", ", green, ", ", blue)
    load = "cmd.load(\""
    f = open(script, "w")
    write(f, "from pymol.cgo import *\n")
    write(f, "from pymol import cmd\n\n")
    write(f, load, inpdb, "\")\n")
    write(f, load, outpdb, "\")\n")
    write(f, load, "modevectors.py\")\n")
    write(f, "rgb=\"", rgb, "\"\n")
    write(f, "modevectors(\"", inpdb[1:end - 4], "\", \"", outpdb[1:end - 4], "\", ")
    write(f, "outname=\"", string(splitext(outpdb)[1], "_porky"), "\", head=0.5, tail=0.3, headrgb = rgb, tailrgb = rgb, cutoff=3.0)\n")
    write(f, "cmd.delete(\"", outpdb[1:end - 4], "\")\n")
    close(f)

end

#Arg Parse settings
s = ArgParseSettings()
@add_arg_table s begin
    "--inpdb", "-p"
    help = "Input PDB"
    arg_type = String
    required = true
    "--vector", "-v"
    help = "Input vector"
    arg_type = String
    default = "none"
    "--outpdb", "-o"
    help = "Output PDB"
    arg_type = String
    required = true
    "--matrix", "-M"
    help = "Input Amber modes"
    arg_type = String
    default = "none"
    "--multiplier", "-m"
    help = "Constant to scale the normalized vector. Defaults to 1"
    arg_type = Float64
    default = 1.
    "--index_min", "-i"
    help = "Mode number range when reading many modes"
    arg_type = Int
    default = 2
    "--index_max", "-I"
    help = "Mode number range when reading many modes"
    arg_type = Int
    default = 1
    "--amber_format", "-a"
    help = "Reading Amber format PCA modes?"
    action = :store_true
    "--script", "-s"
    help = "Only write output PDB"
    action = :store_true
    "--red", "-r"
    help = "Red color for modevectors arrows. If using matrices, a random number will be assigned. Defaults to 1."
    arg_type = Float64
    default = 1.
    "--green", "-g"
    help = "Green color for modevectors arrows. If using matrices, a random number will be assigned. Defaults to 1."
    arg_type = Float64
    default = 1.
    "--blue", "-b"
    help = "Blue color for modevectors arrows. If using matrices, a random number will be assigned. Defaults to 1."
    arg_type = Float64
    default = 1.
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
println("INPDB", "\t",  inpdb)
println("OUTPDB", "\t", outpdb)
println("MULTIPLIER", "\t", multiplier)
println("SCRIPT", "\t", script)

if (vector != "none")
    println("VECTOR", "\t", vector)
    println("RED", "\t", red)
    println("GREEN", "\t", green)
    println("BLUE", "\t", blue)
elseif (matrix != "none")
    println("MATRIX", "\t", matrix)
    println("INDEX_MIN", "\t", index_min)
    println("INDEX_MAX", "\t", index_max)
    println("AMBER_FORMAT", "\t", amber_format)
end

######### For testing purposes. #########
# inpdb = "ef_3fm7.pdb"
# outpdb = "test"
# matrix = "ef_3fm7.mods"
# index_min = 1
# index_max = 2
# script = true
# amber_format = false
# multiplier = 40
# red = .9
# green = .5
# blue = .2
######### For testing purposes. #########


# Get ready
in_vec = Array{Float64,1}
in_vecs = Array{Float64,2}
# Read PDB
const in_trj = Trajectory(inpdb)
const in_frm = read(in_trj)
const in_top = Topology(in_frm)
const aa = convert(Int64, count_residues(in_top))
const aa3 = aa * 3

if (matrix == "none")
    try
        global in_vec = convert(Array{Float64,1}, readdlm(vector)[:, 1])
    catch e
        println(vector, " error:")
        error(e)
    end
else
    if (index_min <= index_max)
        if (amber_format)
            try
                in_modes, in_evals = JUMD.readPtrajModes(matrix, index_max)
                global in_vecs = convert(Array{Float64,2}, in_modes[:, index_min:index_max])
            catch e
                println(matrix, " error:")
                error(e)
            end
        else
            try
                in_modes = convert(Array{Float64,2}, readdlm(matrix))
                global in_vecs = in_modes[:, index_min:index_max]
            catch e
                println(matrix, " error:")
                error(e)
            end
        end
    else
        error("index_min or index_max were not correctly set.")
    end
end

if(size(in_vecs)[1] != aa3)
    error("Size of vector: ", size(in_vecs)[1], " does not match number of residues: ", aa3)
end

if (matrix == "none")
    in_vec = in_vec ./ norm(in_vec) .* multiplier
    out_frm = displaceAA(in_frm, aa, aa3, in_vec);
    # Y guardo
    out_trj = Trajectory(outpdb, 'w')
    write(out_trj, out_frm)
    close(out_trj)

    if script == true
        script_filename = string("script_porky_", splitext(outpdb)[1], ".py")
        write_porcu_script(script_filename, inpdb, outpdb, red, green, blue)
    end
else
    for i = 1:size(in_vecs)[2]
        in_vec = in_vecs[:, i] ./ norm(in_vecs[:, i]) .* multiplier
        out_frm = displaceAA(in_frm, aa, aa3, in_vec);
        # Y guardo
        out_pdb = string(i) * "_" * outpdb
        out_trj = Trajectory(out_pdb, 'w')
        write(out_trj, out_frm)
        close(out_trj)

        if script == true
            script_filename = string("script_porky_", splitext(out_pdb)[1], ".py")
            write_porcu_script(script_filename, inpdb, out_pdb, rand(), rand(), rand())
        end
    end
end


