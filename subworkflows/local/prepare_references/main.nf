include { BWAMEM2_INDEX                   } from '../../../modules/nf-core/bwamem2/index/main'
include { UNTAR as UNTAR_VEP_CACHE        } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_GRIDSS_INDEX     } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HMF_ENSEMBL_DATA } from '../../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {

    take:
    val_genome_fasta            // string: [mandatory] path to reference genome FASTA
    val_bwamem2_index          // string: [optional] path to  pre-built BWA-MEM2 index directory
    val_ensembl_vep_cache      // string: [optional] path to  pre-extracted Ensembl VEP cache directory
    val_ensembl_vep_cache_tar   // string: [optional] path to  compressed Ensembl VEP cache archive
    val_gridss_index            // string: [optional] path to  pre-built GRIDSS index directory
    val_gridss_index_tar       // string: [optional] path to  compressed GRIDSS index archive
    val_hmf_ensembl_data        // string: [optional] path to  pre-extracted HMF Ensembl data directory
    val_hmf_ensembl_data_tar    // string: [optional] path to  compressed HMF Ensembl data archive

    main:
    def ch_genome_fasta     = channel.empty()
    def ch_vep_cache        = channel.empty()
    def ch_gridss_index     = channel.empty()
    def ch_hmf_ensembl_data = channel.empty()
    def ch_versions         = channel.empty()

    ch_genome_fasta = channel.fromPath(val_genome_fasta).map { it -> [[id: it.simpleName], it] }.collect()

    // Genome Indexing

    if (!val_bwamem2_index) {
        BWAMEM2_INDEX ( ch_genome_fasta )

        ch_bwamem2_index = BWAMEM2_INDEX.out.index
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    } else {
        ch_bwamem2_index = channel.fromPath(val_bwamem2_index).map { it -> [[id: it.simpleName], it] }.collect()
    }

    // Ensembl VEP Cache download and untar

    if (!val_ensembl_vep_cache) {
        ch_vep_cache_tar = channel.fromPath(val_ensembl_vep_cache_tar).map { it -> [[id: it.simpleName], it] }.collect()
        UNTAR_VEP_CACHE ( ch_vep_cache_tar )

        ch_vep_cache = UNTAR_VEP_CACHE.out.untar
    } else {
        ch_vep_cache = channel.fromPath(val_ensembl_vep_cache).map { it -> [[id: it.simpleName], it] }.collect()
    }


    // GRIDSS index download and untar
    if (!val_gridss_index) {
        ch_gridss_index_tar = channel.fromPath(val_gridss_index_tar).map { it -> [[id: it.simpleName], it] }.collect()
        UNTAR_GRIDSS_INDEX ( ch_gridss_index_tar )

        ch_gridss_index = UNTAR_GRIDSS_INDEX.out.untar
    } else {
        ch_gridss_index = channel.fromPath(val_gridss_index).map { it -> [[id: it.simpleName], it] }.collect()
    }


    // HMF Ensembl Data download and untar
    if (!val_hmf_ensembl_data) {
        ch_hmf_ensembl_data_tar = channel.fromPath(val_hmf_ensembl_data_tar).map { it -> [[id: it.simpleName], it] }.collect()
        UNTAR_HMF_ENSEMBL_DATA ( ch_hmf_ensembl_data_tar )

        ch_hmf_ensembl_data = UNTAR_HMF_ENSEMBL_DATA.out.untar
    } else {
        ch_hmf_ensembl_data = channel.fromPath(val_hmf_ensembl_data).map { it -> [[id: it.simpleName], it] }.collect()
    }


    emit:
    genome_fasta     = ch_genome_fasta              // channel: [ val(meta), path(fasta) ]
    bwamem2_index    = ch_bwamem2_index             // channel: [ val(meta), path(index) ]
    vep_cache        = ch_vep_cache                 // channel: [ val(meta), path(dir) ]
    gridss_index     = ch_gridss_index              // channel: [ val(meta), path(dir) ]
    hmf_ensembl_data = ch_hmf_ensembl_data          // channel: [ val(meta), path(dir) ]
    versions         = ch_versions                  // channel: [ path(versions.yml) ]
}
