nextflow.enable.dsl = 2
import ReferenceConfigs

include {
    get_curation_sequences;
    get_filter_sequences;
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
    segment_table_add_columns;
    concat_tsv;
    concat_fasta_nolabel;
    mark_sequences_to_filter;
    orient_and_filter_fasta;
    add_unknown_segment_info;
    cdhit;
    build_blast_db;
    write_output_config;
    // finally
    output
} from './lib/processes.nf'


workflow orient_references {
    take:
        segment_info
        // contains label (e.g. LASV), segment (e.g. L), main_ref and other_refs,
        // sequence_file
    main:

        id_lists = segment_info
        | get_id_lists

        sequences_segment = id_lists
        | map {meta, refs, queries -> [meta, refs, queries, meta.fasta]}
        | extract_sequences

        orientations = sequences_segment
        | minimap_refseqs
        | get_orientations_refs

        query_sequences = sequences_segment
        | map { meta, refs, queries -> [meta.label, meta.segment, queries] }

        oriented = orientations
        | map { meta, fwd, rev, unmapped -> [meta.label, meta.segment, meta, fwd, rev, unmapped] }
        | join(query_sequences, by: [0, 1])
        | map {label, segment, meta, fwd, rev, unmapped, queries -> [meta, fwd, rev, unmapped, queries]}        
        | orient_refs

        // some output formatting
        query_ids = id_lists
        | map {meta, ref_ids, query_ids -> [meta, query_ids]}

        write_this = Channel.empty()
        | mix(
            orientations | map {meta, forward_ids, reverse_ids, unmapped_ids -> [unmapped_ids, "ref/${meta.label}", "unmapped_${meta.segment}.txt"]},
            query_ids | map {meta, query_ids -> [query_ids, "ref/${meta.label}", "ids_${meta.segment}.txt"]},
            oriented | map {meta, oriented -> [oriented, "ref/${meta.label}", "oriented_${meta.segment}.fasta"]}
        )
    emit:
        id_lists = query_ids
        oriented_seqs = oriented
        write_this = write_this
}


workflow compare_genomes_to_references {
    take:
        input  // holds label, fasta_ref, segment_table, fasta_all 
    main:
        lengths_references = input
        | map { meta -> [meta.label, meta.fasta_ref] }
        | seq_lengths_references
        
        lengths_queries = input
        | map { meta -> [meta.label, meta.fasta_all] }
        | seq_lengths_queries

        n_shares = input
        | map { meta -> [meta.label, meta.fasta_all] }
        | n_share

        alignments = input
        | map { meta -> [meta.label, meta.fasta_ref, meta.fasta_all] }
        | minimap_all

        align_stats = alignments
        | sam_info

        segment_table = input
        | map {meta -> [meta.label, meta.segment_table]} 

        collected_stats = align_stats
        | join(n_shares, by: 0)
        | join(lengths_queries, by: 0)
        | join(lengths_references, by: 0)
        | join(segment_table, by: 0)
        | collect_seq_and_map_stats

        filtering_stats = collected_stats
        | mark_sequences_to_filter

        sequences = input
        | map { meta -> [meta.label, meta.fasta_all] }

        filtered_sequences = filtering_stats
        | join(sequences, by: 0)
        | orient_and_filter_fasta

        write_this = Channel.empty()
        | mix(
            alignments | map {label, align -> [align, "alignments", "${label}.sam"]},
            filtering_stats | map {label, stats -> [stats, "stats", "${label}.stats.tsv"]},
            filtered_sequences | map {label, fasta -> [fasta, "sequences", "${label}.filtered.fasta"]}
        )
    emit:
        filtered_sequences = filtered_sequences
        write_this = write_this
}


def createFastaDict(List<String> filePaths) {
    def result = [:] // Initialize empty map
    filePaths.each { path ->
        def file = new File(path) // Create File object
        def label = file.name.replace('.fasta', '') // Remove ".fasta" extension
        result[label] = path // Store filename as string
    }
    return result
}


// TODO: add download of organism names to the tax ids


workflow {

    config_channel = Channel.of(params.taxa_config)

    sequences_per_group = config_channel
    | map { fasta -> [fasta, params.fasta_sequences] }
    | get_curation_sequences
    | flatten  // Ensure individual files are passed
    | map { file -> 
        def label = file.baseName.replace('.fasta', '') // Extract identifier
        return tuple(label, file)
    }

    ref_conf = new ReferenceConfigs(params.taxa_config)

    oriented_refs = Channel.from(ref_conf.getRefs())
    | combine(sequences_per_group)
    | filter { meta, label, fasta -> meta.label == label }
    | map { meta, label, fasta -> meta + ['fasta': fasta] }
    | orient_references

    // part two: organize all reads and get the comparisons and statistics
    refs_concat = oriented_refs.oriented_seqs
    | map { meta, oriented_fasta -> [meta.label, oriented_fasta] }
    | groupTuple(by: 0)
    | concat_fasta

    segment_table = oriented_refs.id_lists
    | map { meta, query_ids -> [meta.label, meta.segment, query_ids]}
    | segment_table_add_columns
    | groupTuple(by: 0)
    | concat_tsv

    oriented_seqs = refs_concat
    | join(segment_table, by: 0)
    | join(sequences_per_group, by: 0)
    | map { label, fasta_ref, segment_table, fasta_all ->
        ['label': label,
         'fasta_ref': fasta_ref,
         'segment_table': segment_table,
         'fasta_all': fasta_all]
    }
    | compare_genomes_to_references

    // clustering
    clustered = oriented_seqs.filtered_sequences
    | cdhit

    // merging
    nogroup = sequences_per_group
    | filter { label, fasta -> label == "NOGROUP"}
    | map { label, fasta -> fasta }
    | add_unknown_segment_info

    fasta_all = clustered
    | map { label, fasta -> fasta }
    | mix(nogroup)
    | collect
    | concat_fasta_nolabel

    // extract the family-filters
    family_filters = config_channel
    | combine(fasta_all)
    | get_filter_sequences
    | flatten

    blast_db = fasta_all
    | build_blast_db

    output_config = config_channel
    | write_output_config

    Channel.empty()
    | mix(
        sequences_per_group | map {label, fasta -> [fasta, "sequences", "${label}.all.fasta"]},
        oriented_refs.write_this,
        oriented_seqs.write_this,
        clustered | map { label, fasta -> [fasta, "virus", "${label}.fasta"] },
        fasta_all | map { fasta -> [fasta, "virus", "ALL.fasta"] },
        family_filters | map { fasta -> [fasta, "virus", null] },
        blast_db  | map { blast_dir -> [blast_dir, "virus", "blast_db"] },
        output_config | map { config -> [config, "virus", "virus.yaml"] }
    )
    | output
}
