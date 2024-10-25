nextflow.enable.dsl = 2


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


process collect_seq_and_map_stats {
    input:
        tuple path('orientations.tsv'),
            path('n_share.tsv'),
            path('seq_lengths.tsv'),
            path('seq_lengths_targets.tsv'),
            path('edit_distances.tsv')
    output:
        path('collected_stats.tsv')
    """
    #!/usr/bin/env python

    import pandas as pd

    def tsv(fname, *col_names):
        return pd.read_csv(fname, sep='\\t', header=None, names=col_names)

    n_share = tsv('n_share.tsv', 'Sequence', 'N_share')
    seq_lengths = tsv('seq_lengths.tsv', 'Sequence', 'Length')
    seq_lengths_targets = tsv('seq_lengths_targets.tsv', 'Reference', 'ReferenceLength')
    orientations = tsv('orientations.tsv', 'Sequence', 'Orientation')
    targets_and_distances = tsv('edit_distances.tsv', 'Sequence', 'Reference', 'EditDistance')

    merged_df = (
        n_share
        .merge(seq_lengths, on='Sequence', how='outer')
        .merge(orientations, on='Sequence', how='outer')
        .merge(targets_and_distances, on='Sequence', how='outer')
        .merge(seq_lengths_targets, on='Reference', how='outer')
    )
    merged_df.to_csv('collected_stats.tsv', sep='\\t', index=False)
    """
}

// TODO Filters
// - absolute read length
// - N share
// - un-aligned
//
// - length with respect to reference
// - absolute edit distance (normalized by length)
//
// - edit distance based on error distribution 
//
// - on whole distribution (e.g. 2 times stddev)
// - against hard cutoff (e.g. 60% min similarity)

workflow {

    input_genomes = Channel.of(params.raw_seqs)
    references_concatenated = Channel.from(params.refseqs) | map {entry -> entry.fasta} | collect | concat_fasta

    // TODO: get refence segments to later sort by segment

    lengths_references = seq_lengths_references(references_concatenated)
    lengths_queries = seq_lengths_queries(input_genomes)
    n_shares = n_share(input_genomes)

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
