# Biophysical.jl - Load biophysical constants and build gene info

"""
    load_biophysical_constants(path::String)::BiophysicalConstants

Load cell-free TX/TL biophysical constants from a JSON file.
"""
function load_biophysical_constants(path::String)::BiophysicalConstants
    raw = JSON.parsefile(path)
    bc = raw["biophysical_constants"]
    return BiophysicalConstants(
        parse(Float64, bc["RNAPII_concentration"]["value"]),
        parse(Float64, bc["ribosome_concentration"]["value"]),
        parse(Float64, bc["transcription_elongation_rate"]["value"]),
        parse(Float64, bc["translation_elongation_rate"]["value"]),
        parse(Float64, bc["transcription_saturation_constant"]["value"]),
        parse(Float64, bc["translation_saturation_constant"]["value"]),
        parse(Float64, bc["characteristic_initiation_time_transcription"]["value"]),
        parse(Float64, bc["characteristic_initiation_time_translation"]["value"]),
        parse(Float64, bc["mRNA_half_life_in_hr"]["value"]),
        parse(Float64, bc["protein_half_life_in_hr"]["value"]),
    )
end

"""
    build_gene_info()::GeneInfo

Build the gene information for the gluconate sensor circuit.
Order: [GntR, Venus, sigma_70].
Gene concentrations are initial DNA concentrations in μM.
"""
function build_gene_info()::GeneInfo
    return GeneInfo(
        [:GntR, :Venus, :sigma_70],
        [996.0, 720.0, 720.0],          # coding length (nt)
        [329.0, 238.0, 238.0],          # protein length (aa) ≈ nt/3.03
        [0.010, 0.007, 0.0],            # gene concentration (μM) - sigma_70 expressed exogenously
    )
end
