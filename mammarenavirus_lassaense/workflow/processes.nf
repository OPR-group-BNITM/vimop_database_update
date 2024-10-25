nextflow.enable.dsl = 2


process get_id_lists {
    input:
        val(segment_info)
    output:
        tuple val(segment_info.alias), path('ref_ids.txt'), path('query_ids.txt')
    """
    for id in ${segment_info['seqs'].join(' ')}
    do
        echo \$id >> query_ids.txt
    done
    for id in ${segment_info['refs'].join(' ')}
    do
        echo \$id >> ref_ids.txt
    done
    """
}

process extract_sequences {
    input:
        tuple val(segment), path('ref_ids.txt'), path('query_ids.txt'), path("sequences.fasta")
    output:
        tuple val(segment), path("refs.fasta"), path("queries.fasta")
    """
    seqtk subseq sequences.fasta query_ids.txt > queries.fasta
    seqtk subseq sequences.fasta ref_ids.txt > refs.fasta
    """
}

process minimap {
    input:
        tuple val(alias), path('refs.fasta'), path('queries.fasta')
    output:
        tuple val(alias), path('mapped.sam')
    """
    minimap2 -k ${params.minimap.kmer_length} -w ${params.minimap.window_size} -p ${params.minimap.p} --secondary=no --eqx -x map-ont -a refs.fasta queries.fasta > mapped.sam
    """
}

process get_orientations {
    input:
        tuple val(alias), path('mapped.sam')
    output:
        tuple val(alias), path('forward_mapped_ids.txt'), path('reverse_mapped_ids.txt'), path('unmapped_ids.txt')
    """
    samtools view -F 16 -F 4 mapped.sam | awk '{print \$1}' > forward_mapped_ids.txt
    samtools view -f 16 -F 4 mapped.sam | awk '{print \$1}' > reverse_mapped_ids.txt
    samtools view -f 4 mapped.sam | awk '{print \$1}' > unmapped_ids.txt
    """
}

process orient_reads {
    input:
        tuple val(alias),
            path('forward_mapped_ids.txt'),
            path('reverse_mapped_ids.txt'),
            path('unmapped_ids.txt'),
            path('seqs.fasta')
    output:
        tuple val(alias), path('oriented.fasta')
    """
    #!/usr/bin/env python

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    def read_ids(fname):
        with open(fname) as f:
            return {id.strip() for id in f}

    fwd_ids = read_ids('forward_mapped_ids.txt')
    rev_ids = read_ids('reverse_mapped_ids.txt')
    skip_ids = read_ids('unmapped_ids.txt')

    with open('oriented.fasta', 'w') as out_fasta:
        for record in SeqIO.parse('seqs.fasta', "fasta"):
            if record.id in skip_ids:
                continue
            forward = True
            if record.id in fwd_ids:
                flag = "|forward"
            elif record.id in rev_ids:
                forward = False
                flag = "|reverse"
            else:
                flag = "|ERROR"
            new_record = SeqRecord(
                seq=record.seq if forward else record.seq.reverse_complement(),
                id=record.id,
                name=record.name,
                description=record.description + flag
            )
            SeqIO.write(new_record, out_fasta, "fasta")
    """
}

process seq_lengths {
    input:
        path('seqs.fasta')
    output:
        path('seq_lengths.tsv')
    """
    seqtk comp seqs.fasta | awk '{print \$1"\t"\$2}' > seq_lengths.tsv 
    """
}

process n_share {
    input:
        path('seqs.fasta')
    output:
        path('n_share.tsv')
    """
    seqtk comp seqs.fasta | awk '{print \$1"\t"\$9/\$2}' > n_share.tsv
    """
}

process edit_distances {
    input:
        path('mapped.sam')
    output:
        path('edit_distances.tsv')
    """
    samtools view -F 4 mapped.sam | awk '{for(i=12;i<=NF;i++) if(\$i ~ /^NM:i:/) print \$1"\t"\$3"\t"substr(\$i, 6)}' > edit_distances.tsv
    """
}

process merge_orientations {
    input:
        tuple path('forward_mapped_ids.txt'), path('reverse_mapped_ids.txt'), path('unmapped_ids.txt')
    output:
        path('orientations.tsv')
    """
    awk '{print \$1"\t""forward"}' forward_mapped_ids.txt > orientations.tsv
    awk '{print \$1"\t""reverse"}' reverse_mapped_ids.txt >> orientations.tsv
    awk '{print \$1"\t""unmapped"}' unmapped_ids.txt >> orientations.tsv
    """
}

process concat_fasta {
    input:
        path('in_*.fasta')
    output:
        path('concat.fasta')
    """
    cat in_*.fasta > concat.fasta
    """
}

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

process output {
    // publish inputs to output directory
    cpus 1
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: {
            dirname ? (fname_out ? "$dirname/$fname_out" : "$dirname/$fname_in") : fname_in
        }
    )
    input:
        tuple path(fname_in), val(dirname), val(fname_out)
    output:
        path(fname_in)
    """
    """
}