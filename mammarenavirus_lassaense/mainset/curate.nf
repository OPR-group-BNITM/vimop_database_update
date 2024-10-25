nextflow.enable.dsl = 2


include minimap, get_orientations, orient_reads from '../shared/processes.nf'


process filter_short_seqs {
    input:
        tuple path('seqs.fasta'), val(minlen)
    output:
        tuple path('passed.fasta'), path('failed.txt')
    """
    seqtk comp data/refseqs/oriented_L.fasta | awk '{print \$1 " " \$2}' > lengths.txt
    awk '\$2 < ${minlen} {print \$1}' lengths.txt > failed.txt  
    awk '\$2 >= ${minlen} {print \$1}' lengths.txt > passed.txt
    seqtk subseq seqs.fasta passed.txt > passed.fasta
    """
}

process filter_n_share {
    input:
        tuple path('seqs.fasta'), val(max_n_share)
    output:
        tuple path('passed.fasta'), path('failed.txt')
    """
    seqtk comp data/refseqs/oriented_L.fasta | awk '{print \$1 " " \$9/\$2}' > n_share.txt
    awk '\$2 > ${max_n_share} {print \$1}' n_share.txt > failed.txt
    awk '\$2 =< ${max_n_share} {print \$1}' n_share.txt > passed.txt
    seqtk subseq seqs.fasta passed.txt > passed.fasta
    """
}

process minimap {
    // same as orient_references.nf
    input:
        tuple val(alias), path('refs.fasta'), path('queries.fasta')
    output:
        tuple val(alias), path('mapped.sam')
    """
    minimap2 --secondary=no --eqx -x map-ont -a refs.fasta queries.fasta > mapped.sam
    """
}

process get_orientations {
    // same as orient_references.nf
    input:
        tuple val(alias), path('mapped.sam')
    output:
        tuple val(alias), path('forward_mapped_ids.txt'), path('reverse_mapped_ids.txt'), path('unmapped_ids.txt') 
    """
    samtools view -F 16 mapped.sam | awk '{print \$1}' > forward_mapped_ids.txt
    samtools view -f 16 mapped.sam | awk '{print \$1}' > reverse_mapped_ids.txt
    samtools view -f 4 mapped.sam | awk '{print \$1}' > unmapped_ids.txt
    """
}

process filter_aligned {
    input:
        tuple
    output:
        tuple
    """ TODO (three filtering steps?)
    - get unaligned
    - get badly aligned
      - fixed threshold (e.g. 60 or 70 %) or/AND dynamic (e.g. 2 times std-dev of error rate)
      - make a plot and show it!
    - get too short compared to ref (in %?)
    """
}

// TODO: mapping quality filters
// - on whole distribution (e.g. 2 times stddev)
// - against hard cutoff (e.g. 60% min similarity)

