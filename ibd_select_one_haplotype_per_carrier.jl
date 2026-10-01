"""
ibd_select_one_haplotype_per_carrier.jl: A script that counts carrier pairs
for identity by descent (IBD) segments at the variant's chromosome position
excluding within-family sharing and selects the haplotype with most sharing.

This script runs in parallel to process carrier pairs of several variants.

Written by Gilles-Philippe Morin on 2026-02-12.
"""

using CSV, DataFrames

# 1. Helper to map types for CSV loading
function column_type_mapper(i::Int64, name::Symbol)
    if i == 1 return Int
    elseif i == 2 return String
    elseif i == 3 return Int
    elseif i == 4 return Int32
    else return String
    end
end

function load_variant_data(tped_file::String)
    variant_df = CSV.read(tped_file, DataFrame; delim=' ', header=false, types=column_type_mapper)
    variants = String.(variant_df[:, 2])
    genotype_matrix = Matrix{String}(variant_df[:, 5:end])
    return variants, genotype_matrix
end

function identify_carriers(variants::Vector{String}, genotype_matrix::Matrix{String}, individual_df::DataFrame)
    n_individuals = size(genotype_matrix, 2) ÷ 2
    variant_to_carriers = Dict{String, Set{String}}()
    
    individual_ids = ["$(individual_df[i, 1])_$(individual_df[i, 2])" for i in 1:n_individuals]
    
    for (variant_idx, variant) in enumerate(variants)
        carriers = Set{String}()
        for ind_idx in 1:n_individuals
            # PLINK genotypes: cols 5,6 are Ind1, 7,8 are Ind2...
            # Ignore individuals with missing genotypes
            genotype_matrix[variant_idx, 2*ind_idx-1] == "0" && continue
            genotype_matrix[variant_idx, 2*ind_idx] == "0" && continue
            if genotype_matrix[variant_idx, 2*ind_idx-1] != genotype_matrix[variant_idx, 2*ind_idx]
                push!(carriers, individual_ids[ind_idx])
            end
        end
        if length(carriers) >= 2
            variant_to_carriers[variant] = carriers
        end
    end
    return variant_to_carriers
end

function process_ibd_data(chr::Int, ibd_file::String, variants::Vector{String}, variant_to_carriers::Dict{String, Set{String}})
    output_dir = "haplotype_per_carrier_chr$(chr)"
    !isdir(output_dir) && mkdir(output_dir)
    
    # Pre-parse variant positions to avoid string work in the loop
    # Assumes format "chr1:123456"
    var_positions = Dict(v => parse(Int, split(v, ':')[2]) for v in variants)

    println("Reading IBD file...")
    ibd_df = CSV.read(ibd_file, DataFrame,
        header=["id1", "hap1", "id2", "hap2", "chr", "start", "end", "score", "length"],
        types=[String, Int, String, Int, Int, Int, Int, Float32, Float32])
    
    # Filter IBD data by chromosome if they are merged
    ibd_chr = filter(row -> row.chr == chr, ibd_df)

    # Parallelize at the variant level
    Threads.@threads for variant in variants
        carriers = variant_to_carriers[variant]
        pos = var_positions[variant]
        
        # Local mapping for this variant
        carrier_list = collect(carriers)
        carrier_to_idx = Dict(c => i for (i, c) in enumerate(carrier_list))
        pairs_matrix = zeros(Int32, length(carrier_list), 2)

        # Iterate through chromosome-specific IBD segments
        for row in eachrow(ibd_chr)
            if row.start <= pos <= row.end
                if haskey(carrier_to_idx, row.id1) && haskey(carrier_to_idx, row.id2)
                    # Ignore IBD segments between individuals from the same family
                    split(row.id1, '_')[2] == split(row.id2, '_')[2] && continue
                    # Use row.hap to index (ensure they are 1 or 2)
                    pairs_matrix[carrier_to_idx[row.id1], row.hap1] += 1
                    pairs_matrix[carrier_to_idx[row.id2], row.hap2] += 1
                end
            end
        end

        split_haplotypes(row::AbstractArray{Int32}) = (row[1], row[2])
        select_haplotype(h1, h2) = h1 == h2 ? (h1 == 0 ? 0 : 3) : (h1 > h2 ? 1 : 2)
        selected_haplotypes = [select_haplotype(r[1], r[2]) for r in eachrow(pairs_matrix)]
        
        # Reset the pairs matrix for a second run
        pairs_matrix = zeros(Int32, length(carriers), 2)

        # Run a second time with the preferred haplotype
        for row in eachrow(ibd_chr)
            if row.start <= pos <= row.end
                if haskey(carrier_to_idx, row.id1) && haskey(carrier_to_idx, row.id2)
                    # Ignore IBD segments between individuals from the same family
                    split(row.id1, '_')[2] == split(row.id2, '_')[2] && continue
                    # Use row.hap to index (ensure they are 1 or 2)
                    sel_hap1 = selected_haplotypes[carrier_to_idx[row.id1]]
                    (row.hap1 != sel_hap1 && sel_hap1 != 3) && continue
                    sel_hap2 = selected_haplotypes[carrier_to_idx[row.id2]]
                    (row.hap2 != sel_hap2 && sel_hap1 != 3) && continue
                    pairs_matrix[carrier_to_idx[row.id1], row.hap1] += 1
                    pairs_matrix[carrier_to_idx[row.id2], row.hap2] += 1
                end
            end
        end

        # Write result for this variant
        open(joinpath(output_dir, "$(replace(variant, ":"=>"_")).tsv"), "w") do io
            for (i, carrier) in enumerate(carrier_list)
                h1, h2 = pairs_matrix[i, 1], pairs_matrix[i, 2]
                selected = h1 == h2 ? (h1 == 0 ? 0 : 3) : (h1 > h2 ? 1 : 2)
                println(io, "$carrier\t$h1\t$h2\t$selected")
            end
        end
    end
end

function analyze_variants(chr::Int, tped_file::String, tfam_file::String, ibd_file::String)
    variants, genotype_matrix = load_variant_data(tped_file)
    individual_df = CSV.read(tfam_file, DataFrame; header=false)
    
    variant_to_carriers = identify_carriers(variants, genotype_matrix, individual_df)
    
    # Only process variants that actually have carriers
    active_variants = intersect(variants, keys(variant_to_carriers))
    
    process_ibd_data(chr, ibd_file, active_variants, variant_to_carriers)
    println("Analysis complete!")
end

function main()
    if length(ARGS) != 4
        println("Usage: julia script.jl <chr> <tped> <tfam> <ibd.gz>")
        exit(1)
    end
    analyze_variants(parse(Int, ARGS[1]), ARGS[2], ARGS[3], ARGS[4])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
