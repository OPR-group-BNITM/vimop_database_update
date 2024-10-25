nextflow.enable.dsl = 2


include {
    // reference orientation
    minimap as minimap_refseqs;
    orient_reads as orient_refs;
    extract_sequences;
    get_id_lists;
    get_orientations as get_orientations_refs;    
    // in between
    concat_fasta;
    // second part
    minimap as minimap_all;
    get_orientations as get_orientations_all;
    orient_reads;
    seq_lengths as seq_lengths_references;
    seq_lengths as seq_lengths_queries;
    n_share;
    edit_distances;
    merge_orientations;
    collect_seq_and_map_stats;
    // finally
    output
} from './processes.nf'


workflow orient_references {
    take:
        segment_info
        sequences
    main:
        id_lists = segment_info | get_id_lists
        sequences = id_lists.combine(sequences) | extract_sequences
        orientations = sequences | minimap_refseqs | get_orientations_refs
        oriented = orientations | join(sequences | map{segment, refs, queries -> [segment, queries]}, by: 0) | orient_refs
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


workflow compare_genomes_to_references {
    take:
        refseqs
        sequences
    main:
        lengths_references = seq_lengths_references(refseqs)
        lengths_queries = seq_lengths_queries(sequences)
        n_shares = n_share(sequences)

        alignments = refseqs | map {refs -> ['all', refs, sequences]} | minimap_all
        edists = alignments | map {alias, align -> align} | edit_distances
        
        orientations = alignments | get_orientations_all
        
        merged_orientations = orientations
        | map {alias, fwd_id, rev_ids, unmapped -> [fwd_id, rev_ids, unmapped]}
        | merge_orientations

        collected_stats = merged_orientations
        | combine(n_shares)
        | combine(lengths_queries)
        | combine(lengths_references)
        | combine(edists)
        | collect_seq_and_map_stats

        write_this = Channel.empty()
        | mix(
            collected_stats | map {stats -> [stats, 'all', null]}
        )
    emit:
        write_this = write_this
}


workflow {

    segment_info = Channel.from(params.segments)
    sequences = Channel.of(params.fasta)

    // part one: organize the references
    oriented_refs = orient_references(segment_info, sequences)

    // part two: organize all reads and get the comparisons and statistics
    // TODO: export to sub-workflow!
    refs_concat = oriented_refs.oriented_seqs
    | map {segment, oriented_fasta -> oriented_fasta}
    | collect
    | concat_fasta

    // Todo: merge id lists and add segment info

    oriented_seqs = compare_genomes_to_references(refs_concat, params.fasta)

    Channel.empty()
    | mix(
        oriented_refs.write_this,
        oriented_seqs.write_this
    )
    | output
}
