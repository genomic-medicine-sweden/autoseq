include { BWAMEM2_INDEX                   } from '../../../modules/nf-core/bwamem2/index/main'
include { UNTAR as UNTAR_GRIDSS_INDEX     } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HMF_ENSEMBL_DATA } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_VEP_CACHE        } from '../../../modules/nf-core/untar/main'
include { channelFromPathWithMeta         } from '../utils_nfcore_autoseq_pipeline'

workflow PREPARE_REFERENCES {

    take:
    val_genome_fasta      // string: [mandatory] path to reference genome FASTA
    val_bwamem2_index     // string: [optional]  path to pre-built BWA-MEM2 index directory
    val_ensembl_vep_cache // string: [mandatory] path to Ensembl VEP cache directory or .tar.gz archive
    val_gridss_index      // string: [mandatory] path to GRIDSS index directory or .tar.gz archive
    val_hmf_ensembl_data  // string: [mandatory] path to HMF Ensembl data directory or .tar.gz archive

    main:
    def ch_bwamem2_index    = channel.empty()
    def ch_genome_fasta     = channel.empty()
    def ch_vep_cache        = channel.empty()
    def ch_gridss_index     = channel.empty()
    def ch_hmf_ensembl_data = channel.empty()

    ch_genome_fasta = channelFromPathWithMeta(val_genome_fasta)

    // Genome Indexing

    if (!val_bwamem2_index) {
        BWAMEM2_INDEX ( ch_genome_fasta )

        ch_bwamem2_index = BWAMEM2_INDEX.out.index
    } else {
        ch_bwamem2_index = channelFromPathWithMeta(val_bwamem2_index)
    }

    // Ensembl VEP cache: untar when given a .tar.gz archive, otherwise use the directory as is
    if (val_ensembl_vep_cache.endsWith('.tar.gz')) {
        ch_vep_cache_tar = channelFromPathWithMeta(val_ensembl_vep_cache)
        UNTAR_VEP_CACHE ( ch_vep_cache_tar )

        ch_vep_cache = UNTAR_VEP_CACHE.out.untar.collect()
    } else {
        ch_vep_cache = channelFromPathWithMeta(val_ensembl_vep_cache)
    }

    // GRIDSS index: untar when given a .tar.gz archive, otherwise use the directory as is
    if (val_gridss_index.endsWith('.tar.gz')) {
        ch_gridss_index_tar = channelFromPathWithMeta(val_gridss_index)
        UNTAR_GRIDSS_INDEX ( ch_gridss_index_tar )

        ch_gridss_index = UNTAR_GRIDSS_INDEX.out.untar.collect()
    } else {
        ch_gridss_index = channelFromPathWithMeta(val_gridss_index)
    }

    // HMF Ensembl data: untar when given a .tar.gz archive, otherwise use the directory as is
    if (val_hmf_ensembl_data.endsWith('.tar.gz')) {
        ch_hmf_ensembl_data_tar = channelFromPathWithMeta(val_hmf_ensembl_data)
        UNTAR_HMF_ENSEMBL_DATA ( ch_hmf_ensembl_data_tar )

        ch_hmf_ensembl_data = UNTAR_HMF_ENSEMBL_DATA.out.untar.collect()
    } else {
        ch_hmf_ensembl_data = channelFromPathWithMeta(val_hmf_ensembl_data)
    }


    emit:
    genome_fasta     = ch_genome_fasta     // channel: [ val(meta), path(fasta) ]
    bwamem2_index    = ch_bwamem2_index    // channel: [ val(meta), path(index) ]
    vep_cache        = ch_vep_cache        // channel: [ val(meta), path(dir) ]
    gridss_index     = ch_gridss_index     // channel: [ val(meta), path(dir) ]
    hmf_ensembl_data = ch_hmf_ensembl_data // channel: [ val(meta), path(dir) ]
}
