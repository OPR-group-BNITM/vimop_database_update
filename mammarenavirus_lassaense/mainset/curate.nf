nextflow.enable.dsl = 2


// TODO: delete this !!!

include { 
    minimap;
    get_orientations;
    orient_reads;
    output;
    seq_lengths as seq_lengths_references;
    seq_lengths as seq_lengths_queries;
    n_share;
    edit_distances;
    merge_orientations;
    concat_fasta
} from '../shared/processes.nf'


workflow {

    input_genomes = Channel.of(params.raw_seqs)
    references_concatenated = Channel.from(params.refseqs) | map {entry -> entry.fasta} | collect | concat_fasta

    lengths_references = seq_lengths_references(references_concatenated)
    lengths_queries = seq_lengths_queries(input_genomes)
    n_shares = n_share(input_genomes)

    // TODO: move this to the new combined script!

    alignments = references_concatenated | map {refs -> ['all', refs, params.raw_seqs]} | minimap
    edists = alignments | map {alias, align -> align} | edit_distances
    
    orientations = alignments | get_orientations
    
    merged_orientations = orientations
    | map {alias, fwd_id, rev_ids, unmapped -> [fwd_id, rev_ids, unmapped]}
    | merge_orientations

    collected_stats = merged_orientations
    | combine(n_shares)
    | combine(lengths_queries)
    | combine(lengths_references)
    | combine(edists)
    | collect_seq_and_map_stats

    Channel.empty()
    | mix(
        lengths_references,
        lengths_queries,
        n_shares,
        alignments | map {alias, align -> align},
        edists,
        merged_orientations,
        collected_stats,
        orientations | map {alias, fwd_id, rev_ids, unmapped -> unmapped},
        orientations | map {alias, fwd_id, rev_ids, unmapped -> rev_ids},
        orientations | map {alias, fwd_id, rev_ids, unmapped -> fwd_id}
    )
    | map {stats -> [stats, '.', null]}
    | output
}
