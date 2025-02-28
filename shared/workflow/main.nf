nextflow.enable.dsl = 2


include {
    get_sequence_sets_for_curation;
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
    sam_info;
    seq_lengths as seq_lengths_references;
    seq_lengths as seq_lengths_queries;
    n_share;
    collect_seq_and_map_stats;
    add_col;
    concat_tsv;
    // finally
    output
} from './lib/processes.nf'


workflow orient_references {
    take:
        segment_info
        // contains label (e.g. LASV), segment (e.g. L), main_ref and other_refs,
        // sequence_file, sequence_index_file
    main:
        id_lists = segment_info
        | get_id_lists

        sequences_segment = id_lists.combine([segment_info.fasta])
        | extract_sequences

        orientations = sequences_segment
        | minimap_refseqs
        | get_orientations_refs

        oriented = orientations
        | join(sequences_segment | map{segment, refs, queries -> [segment, queries]}, by: 0)
        | orient_refs

        // some output formatting
        query_ids = id_lists
        | map {segment, ref_ids, query_ids -> [segment, query_ids]}

        write_this = Channel.empty()
        | mix(
            orientations | map {meta, forward_ids, reverse_ids, unmapped_ids -> [unmapped_ids, "${meta.label}/ref", "unmapped_${meta.segment}.txt"]},
            query_ids | map {meta, query_ids -> [query_ids, "${meta.label}/ref", "ids_${meta.segment}.txt"]},
            oriented | map {meta, oriented -> [oriented, "${meta.label}/ref", "oriented_${meta.segment}.fasta"]}
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
        segment_table
    main:
        lengths_references = seq_lengths_references(refseqs)
        lengths_queries = seq_lengths_queries(sequences)
        n_shares = n_share(sequences)

        alignments = refseqs | map {refs -> ['all', refs, sequences]} | minimap_all        
        align_stats = alignments | map {alias, align -> [align]} | sam_info

        collected_stats = align_stats
        | combine(n_shares)
        | combine(lengths_queries)
        | combine(lengths_references)
        | combine(segment_table)
        | collect_seq_and_map_stats

        write_this = Channel.empty()
        | mix(
            alignments | map {alias, align -> [align, 'all/alignments', "${alias}.sam"]},
            collected_stats | map {stats -> [stats, 'all', null]}
        )
    emit:
        write_this = write_this
}


workflow {

    sequences_per_group = Channel.of(params.taxa_config)
    | map {fasta -> [fasta, params.fasta_sequences]}
    | get_sequence_sets_for_curation
    | flatten()  // Ensure individual files are passed
    | map { file -> 
        def label = file.baseName.replace('.fasta', '') // Extract identifier
        return tuple(label, file)
    }

    sequences_dict = sequences_per_group
    | map {label, fasta -> fasta}
    | flatten()  // Ensure individual files are passed
    | reduce([:]) { acc, file ->
        def label = file.baseName.replace('.fasta', '') // Extract identifier
        acc[label] = file
        return acc
    }

    // TODO
    // - get the info for the groups
    // - merge the sequences per group

    // map sequences per group -> config[label].ref, config[label].

    // LASV:
    // segments:
    // S:
    //   refs:
    //     - NC_004296.1
    //   seqs:

    //segment_info = Channel.from(params.segments)

    // // part one: organize the references
    // oriented_refs = orient_references(segment_info, params.fasta_refs)

    // // part two: organize all reads and get the comparisons and statistics
    // refs_concat = oriented_refs.oriented_seqs
    // | map {segment, oriented_fasta -> oriented_fasta}
    // | collect  // groupby!
    // | concat_fasta

    // segment_table = oriented_refs.id_lists
    // | add_col
    // | collect
    // | concat_tsv

    // oriented_seqs = compare_genomes_to_references(refs_concat, params.fasta_seqs, segment_table)

    Channel.empty()
    | mix(
        sequences_per_group | map {label, fasta -> [fasta, "sequences", null]}
        // oriented_refs.write_this,
        // oriented_seqs.write_this
    )
    | output
}
