include { BWAMEM2_INDEX                   } from '../../../modules/nf-core/bwamem2/index/main'
include { UNTAR as UNTAR_VEP_CACHE        } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_GRIDSS_INDEX     } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HMF_ENSEMBL_DATA } from '../../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {

    take:
    ch_genome_fasta          // channel: [ val(meta), path(fasta) ]        genome FASTA to build the BWA-MEM2 index from
    ch_vep_cache_tar         // channel: [ val(meta), path(archive) ]      VEP cache tar archive
    ch_gridss_index_tar      // channel: [ val(meta), path(archive) ]      GRIDSS index tar archive
    ch_hmf_ensembl_data_tar  // channel: [ val(meta), path(archive) ]      HMF ensembl_data tar archive

    main:

    def ch_versions = channel.empty()

    //
    // MODULE: Build the BWA-MEM2 index on the fly from the genome FASTA
    //
    BWAMEM2_INDEX ( ch_genome_fasta )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    //
    // MODULE: Untar directory-type references
    //
    UNTAR_VEP_CACHE ( ch_vep_cache_tar )

    UNTAR_GRIDSS_INDEX ( ch_gridss_index_tar )

    UNTAR_HMF_ENSEMBL_DATA ( ch_hmf_ensembl_data_tar )

    emit:
    bwamem2_index    = BWAMEM2_INDEX.out.index             // channel: [ val(meta), path(index) ]
    vep_cache        = UNTAR_VEP_CACHE.out.untar           // channel: [ val(meta), path(dir) ]
    gridss_index     = UNTAR_GRIDSS_INDEX.out.untar        // channel: [ val(meta), path(dir) ]
    hmf_ensembl_data = UNTAR_HMF_ENSEMBL_DATA.out.untar    // channel: [ val(meta), path(dir) ]
    versions         = ch_versions                         // channel: [ path(versions.yml) ]
}
