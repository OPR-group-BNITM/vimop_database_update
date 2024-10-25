nextflow.enable.dsl = 2


include { minimap; get_orientations; orient_reads; output } from '../shared/processes.nf'


process extract_sequences {
    input:
        tuple val(segment_info), path('sequences.fasta')
    output:
        tuple val(segment_info.alias), path("refs.fasta"), path("queries.fasta")
    """
    for id in ${segment_info['seqs'].join(' ')}
    do
        echo \$id >> seq_ids_queries.txt
    done
    seqtk subseq sequences.fasta seq_ids_queries.txt > queries.fasta
    for id in ${segment_info['refs'].join(' ')}
    do
        echo \$id >> seq_ids_refs.txt
    done
    seqtk subseq sequences.fasta seq_ids_refs.txt > refs.fasta
    """
}


workflow {
    segment_info = Channel.from(params.segments).combine(Channel.of(params.fasta))

    sequences = segment_info | extract_sequences
    mapped = sequences | minimap
    orientations = mapped | get_orientations
    oriented_reads = orientations | join(sequences | map{alias, refs, queries -> [alias, queries]}, by: 0) | orient_reads

    Channel.empty()
    | mix(
        orientations | map {alias, forward_ids, reverse_ids, unmapped_ids -> [unmapped_ids, '.', "unmapped_${alias}.txt"]},
        orientations | map {alias, forward_ids, reverse_ids, unmapped_ids -> [forward_ids, '.', "forward_ids_${alias}.txt"]},
        orientations | map {alias, forward_ids, reverse_ids, unmapped_ids -> [reverse_ids, '.', "reverse_ids_${alias}.txt"]},
        oriented_reads | map {alias, oriented -> [oriented, '.', "oriented_${alias}.fasta"]}
    )
    | output
}
