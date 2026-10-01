"""
identify_carriers.jl: A script that prints every carrier for each variant.

Written by Gilles-Philippe Morin.
"""
identify_carriers.jl
#!/usr/bin/env julia

using CSV, DataFrames

"""
    column_type_mapper(i::Int64, name::Symbol)
"""
function column_type_mapper(i::Int64, name::Symbol)
    if i == 1
        return Int
    elseif i == 2
        return String
    elseif i == 3
        return Int
    elseif i == 4
        return Int32
    else
        return String
    end
end

"""
    load_variant_data(tped_file::String)

Load variant data from PLINK .tped file and return variants and genotype matrix.
"""
function load_variant_data(tped_file::String)
    variant_df = CSV.read(tped_file, DataFrame;
        delim=' ', header=false, types=column_type_mapper)
    variants = String.(variant_df[:, 2])
    @views genotype_matrix = Matrix{String}(variant_df[:, 5:end])
    return variants, genotype_matrix
end

"""
    identify_carriers_for_selected(variants::Vector{String}, genotype_matrix::Matrix{String},
                                  individual_df::DataFrame, selected_variants::Set{String})

Identify carriers for selected variants only.
"""
function identify_carriers_for_selected(variants::Vector{String}, genotype_matrix::Matrix{String},
                                       individual_df::DataFrame, selected_variants::Set{String})

    n_variants = length(variants)
    n_individuals = size(genotype_matrix, 2) ÷ 2
    variant_to_carriers = Dict{String, Vector{String}}()

    # Build individual IDs
    individual_ids = ["$(individual_df[i, 1])_$(individual_df[i, 2])" for i in 1:n_individuals]

    for (variant_idx, variant) in enumerate(variants)
        # Skip if variant not in selected list
        variant ∈ selected_variants || continue

        carriers = String[]

        for ind_idx in 1:n_individuals
            allele1 = genotype_matrix[variant_idx, 2 * ind_idx - 1]
            allele2 = genotype_matrix[variant_idx, 2 * ind_idx]

            # Hétérozygotes (porteurs d’un seul allèle ALT)
            if allele1 != allele2
                push!(carriers, individual_ids[ind_idx])
            end
        end

        variant_to_carriers[variant] = carriers
    end

    return variant_to_carriers
end

"""
    main()
"""
function main()
    if length(ARGS) != 5
        println("Usage: julia get_carriers.jl <chromosome> <tped_file> <tfam_file> <variants_list.txt> <output.txt>")
        exit(1)
    end

    chr = parse(Int, ARGS[1])
    tped_file = ARGS[2]
    tfam_file = ARGS[3]
    variant_list_file = ARGS[4]
    output_file = ARGS[5]

    println("Loading variant data...")
    variants, genotype_matrix = load_variant_data(tped_file)

    println("Loading individuals...")
    individual_df = CSV.read(tfam_file, DataFrame; header=false, delim=' ')

    println("Loading selected variants...")
    selected_variants = Set(String.(strip.(readlines(variant_list_file))))

    println("Identifying carriers...")
    variant_to_carriers = identify_carriers_for_selected(variants, genotype_matrix, individual_df, selected_variants)

   println("Writing output...")

# Vérifie si le fichier existe déjà
append_mode = isfile(output_file)

open(output_file, append_mode ? "a" : "w") do io
    # Si le fichier n’existe pas encore, écrire l'en-tête
    if !append_mode
        println(io, "variant\tcarrier")
    end
    for (variant, carriers) in variant_to_carriers
        for carrier in carriers
            println(io, "$(variant)\t$(carrier)")
        end
    end
end

println("Done! Output written to $(output_file)")
end

# Run if script executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end



