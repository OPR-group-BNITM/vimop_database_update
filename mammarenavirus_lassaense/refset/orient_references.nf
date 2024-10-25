nextflow.enable.dsl = 2


include {
    get_id_lists;
    extract_sequences;
    minimap as minimap_refseqs;
    get_orientations;
    orient_reads;
    output
} from '../shared/processes.nf'


workflow orient_references {
    take:
        segment_info
        sequences
    main:
        id_lists = segment_info | get_id_lists
        sequences = id_lists.combine(sequences) | extract_sequences
        orientations = sequences | minimap_refseqs | get_orientations
        oriented = orientations | join(sequences | map{segment, refs, queries -> [segment, queries]}, by: 0) | orient_reads
        // some output formatting
        query_ids = id_lists | map {segment, ref_ids, query_ids -> [segment, query_ids]}
        write_this = Channel.empty()
        | mix(
            orientations | map {segment, forward_ids, reverse_ids, unmapped_ids -> [unmapped_ids, "ref", "unampped_${segment}.txt"]},
            query_ids | map {segment, query_ids -> [query_ids, "ref", "ids_${segment}.txt"]},
            oriented | map {segment, oriented -> [oriented, "ref", "oriented_${segment}.fasta"]}
        )
    emit:
        id_lists = query_ids
        oriented_seqs = oriented
        write_this = write_this
}


workflow {

    segment_info = Channel.from(params.segments)
    sequences = Channel.of(params.fasta)

    // part one: organize the references
    oriented_refs = orient_references(segment_info, sequences)

    // part two: organize all reads
    

    Channel.empty()
    | mix(
        oriented_refs.write_this
    )
    | output
}
