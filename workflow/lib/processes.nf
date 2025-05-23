nextflow.enable.dsl = 2


process get_curation_sequences {
    label 'general'
    input:
        tuple path('groups.yaml'), path('sequences.fasta')
    output:
        path("*.fasta")
    """
    sort_ids_to_groups.py \\
        --taxgroups groups.yaml \\
        --sequences sequences.fasta \\
        --outdir out \\
        --remove_set ${params.remove_set} \\
        --category curated

    for fname_ids in out/*.txt
    do
        bname=\$(basename \$fname_ids)
        fname_out="\${bname%.txt}.fasta" # Remove .txt extension
        seqtk subseq sequences.fasta \$fname_ids > \$fname_out
    done
    """
}


process get_filter_sequences {
    label 'general'
    input:
        tuple path('groups.yaml'), path('sequences.fasta')
    output:
        path("out/*.fasta")
    """
    sort_ids_to_groups.py \\
        --taxgroups groups.yaml \\
        --sequences sequences.fasta \\
        --outdir tmp \\
        --remove_set ${params.remove_set} \\
        --category filters \\
        --organism_label_position 3

    if [ -f tmp/NOGROUP.txt ]
    then
        rm tmp/NOGROUP.txt
    fi

    for fname_ids in tmp/*.txt
    do
        bname=\$(basename \$fname_ids)
        fname_out="\${bname%.txt}.fasta" # Remove .txt extension
        seqtk subseq sequences.fasta \$fname_ids > tmp/\$fname_out
    done
    mv tmp out
    """
}


process get_id_lists {
    label 'general'
    input:
        val(segment_info)
    output:
        tuple val(segment_info), path('ref_ids.txt'), path('query_ids.txt')
    """
    for id in ${segment_info.refs.join(' ')}
    do
        echo \$id >> ref_ids.txt
    done
    for id in ${segment_info.seqs.join(' ')}
    do
        echo \$id >> query_ids.txt
    done
    """
}


process extract_sequences {
    label 'general'
    input:
        tuple val(segment_info), path('ref_ids.txt'), path('query_ids.txt'), path("sequences.fasta")
    output:
        tuple val(segment_info), path("refs.fasta"), path("queries.fasta")
    """
    seqtk subseq sequences.fasta query_ids.txt > queries.fasta
    seqtk subseq sequences.fasta ref_ids.txt > refs.fasta
    """
}


process minimap {
    memory '32 GB'
    label 'general'
    input:
        tuple val(meta), path('refs.fasta'), path('queries.fasta')
    output:
        tuple val(meta), path('mapped.sam')
    """
    set -o pipefail
    minimap2 \\
        -k ${params.minimap.kmer_length} \\
        -w ${params.minimap.window_size} \\
        -p ${params.minimap.p} \\
        --secondary=no \\
        --eqx \\
        -x map-ont \\
        -a refs.fasta queries.fasta -o mapped.sam
    """
}


process get_orientations {
    label 'general'
    input:
        tuple val(meta), path('mapped.sam')
    output:
        tuple val(meta), path('forward_mapped_ids.txt'), path('reverse_mapped_ids.txt'), path('unmapped_ids.txt')
    """
    # -F 2308 removes unmapped reads (4), secondary alignments (256) and supplementary alignments (2048) 
    samtools view -F 16 -F 2308 mapped.sam | awk '{print \$1}' > forward_mapped_ids.txt
    samtools view -f 16 -F 2308 mapped.sam | awk '{print \$1}' > reverse_mapped_ids.txt
    samtools view -f 4 mapped.sam | awk '{print \$1}' > unmapped_ids.txt
    """
}

process orient_reads {
    label 'general'
    input:
        tuple val(meta),
            path('forward_mapped_ids.txt'),
            path('reverse_mapped_ids.txt'),
            path('unmapped_ids.txt'),
            path('seqs.fasta')
    output:
        tuple val(meta), path('oriented.fasta')
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


process segment_table_add_columns {
    label 'general'
    input:
        tuple val(label), val(segment), path('input_table_to_extend.tsv')
    output:
        tuple val(label), path('table_out.tsv')
    """
    awk -v label="${label}" -v segment="${segment}" '{print \$0"\\t"label"\\t"segment}' input_table_to_extend.tsv > table_out.tsv
    """
}

process concat_tsv {
    label 'general'
    input:
        tuple val(label), path('in_*.tsv')
    output:
        tuple val(label), path('out.tsv')
    """
    cat in_*.tsv > out.tsv
    """
}

process concat_fasta_nolabel {
    label 'general'
    input:
        path('*.fasta')
    output:
        path('out.fasta')
    """
    cat *.fasta > out.fasta
    """
} 

process seq_lengths {
    label 'general'
    input:
        tuple val(label), path('seqs.fasta')
    output:
        tuple val(label), path('seq_lengths.tsv')
    """
    seqtk comp seqs.fasta | awk '{print \$1"\t"\$2}' > seq_lengths.tsv 
    """
}

process n_share {
    label 'general'
    input:
        tuple val(label), path('seqs.fasta')
    output:
        tuple val(label), path('n_share.tsv')
    """
    seqtk comp seqs.fasta | awk '{print \$1"\t"\$9/\$2}' > n_share.tsv
    """
}

process sam_info {
    label 'general'
    input:
        tuple val(meta), path('mapped.sam')
    output:
        tuple val(meta), path('stats.tsv')
    """
    #!/usr/bin/env python
    import pysam
    import pandas as pd

    def cast_int(val):
        try:
            return int(val)
        except TypeError:
            return -1

    with pysam.AlignmentFile("mapped.sam", "r") as samfile:
        data = [
            {
                'Sequence': read.query_name,
                'Reference': samfile.get_reference_name(read.reference_id),
                'IsForward': not read.is_reverse, 
                'ReferenceStart': cast_int(read.reference_start),
                'ReferenceEnd': cast_int(read.reference_end),
                'QueryStart': cast_int(read.query_alignment_start),
                'QueryEnd': cast_int(read.query_alignment_end),
                'EditDistance': cast_int(dict(read.tags).get('NM')),
                'IsSupplementaryAlignment': read.is_supplementary,  
            }
            for read in samfile.fetch()
            if not read.is_supplementary and not read.is_secondary
        ]
    columns = [
        "Sequence",
        "Reference",
        "IsForward", 
        "ReferenceStart",
        "ReferenceEnd", 
        "QueryStart",
        "QueryEnd", 
        "EditDistance",
        "IsSupplementaryAlignment"
    ]
    pd.DataFrame(data, columns=columns).to_csv('stats.tsv', sep='\\t')
    """
}

process concat_fasta {
    label 'general'
    input:
        tuple val(label), path('in_*.fasta')
    output:
        tuple val(label), path('concat.fasta')
    """
    cat in_*.fasta > concat.fasta
    """
}

process collect_seq_and_map_stats {
    label 'general'
    input:
        tuple val(label),
            path('samstats.tsv'),
            path('n_share.tsv'),
            path('seq_lengths.tsv'),
            path('seq_lengths_targets.tsv'),
            path('segment_table.tsv')
    output:
        tuple val(label), path('collected_stats.tsv')
    """
    #!/usr/bin/env python

    import pandas as pd

    def tsv(fname, **cols_dtypes):
        return pd.read_csv(fname, sep='\\t', header=None, names=list(cols_dtypes), dtype=cols_dtypes)

    dtypes_samstats = {
        'Sequence': 'str',
        'Reference': 'str',
        'IsForward': 'bool',
        'ReferenceStart': 'int',
        'ReferenceEnd': 'int',
        'QueryStart': 'int',
        'QueryEnd': 'int',
        'EditDistance': 'int',
        'IsSupplementaryAlignment': 'bool',
    }
    samstats = pd.read_csv('samstats.tsv', sep='\\t', dtype=dtypes_samstats)
    n_share = tsv('n_share.tsv', Sequence='str', N_share='float')
    seq_lengths = tsv('seq_lengths.tsv', Sequence='str', Length='int')
    seq_lengths_targets = tsv('seq_lengths_targets.tsv', Reference='str', ReferenceLength='int')
    segments = tsv('segment_table.tsv', Reference='str', Label='str', Segment='str')

    merged_df = (
        samstats
        .merge(n_share, on='Sequence', how='outer')
        .merge(seq_lengths, on='Sequence', how='outer')
        .merge(seq_lengths_targets, on='Reference', how='outer')
        .merge(segments, on='Reference', how='outer')
    )
    merged_df.to_csv('collected_stats.tsv', sep='\\t', index=False)
    """
}


process mark_sequences_to_filter {
    label 'general'
    input:
        tuple val(label), path('collected_stats.tsv')
    output:
        tuple val(label), path('collected_stats_filtered.tsv')
    """
    #!/usr/bin/env python3
    import pandas as pd
    seqstats = pd.read_csv('collected_stats.tsv', sep='\\t')
    seqstats['RelativeLength'] = seqstats['Length'] / seqstats['ReferenceLength']

    unmapped = seqstats['ReferenceStart'] == -1
    too_many_N = seqstats['N_share'] > ${params.filter_max_n_share}
    too_short = seqstats['RelativeLength'] < ${params.filter_min_relative_length}

    seqstats['FilteringStatus'] = 'Ok'

    seqstats['FilteringStatus'][too_short] = 'TooShort'
    seqstats['FilteringStatus'][unmapped] = 'Unmapped'
    seqstats['FilteringStatus'][too_many_N] = 'TooManyN'

    seqstats.to_csv('collected_stats_filtered.tsv', sep='\\t')
    """
}


process orient_and_filter_fasta {
    label 'general'
    input:
        tuple val(label), path('collected_stats_filtered.tsv'), path('sequences.fasta')
    output:
        tuple val(label), path("${label}.filtered.fasta")
    """
    #!/usr/bin/env python3

    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    df = pd.read_csv('collected_stats_filtered.tsv', sep='\\t')
    df_filtered_final = df[df['FilteringStatus'] == 'Ok']

    segment_and_orientation = {
        row['Sequence']: (row['IsForward'], row['Segment'])
        for _, row in df_filtered_final.iterrows()
    }

    with open('${label}.filtered.fasta', 'w') as f_out:
        for record in SeqIO.parse('sequences.fasta', 'fasta'):
            if record.id not in segment_and_orientation:
                continue
            is_forward, segment = segment_and_orientation[record.id]
            if is_forward is None:
                raise RuntimeError(f"Invalid orientation {is_forward}")
            orientation = 'forward' if is_forward else 'reverse'
            new_record = SeqRecord(
                seq=record.seq if is_forward else record.seq.reverse_complement(),
                id=record.id,
                name=record.name,
                description=f'{record.description}|{orientation}|{segment}'
            )
            SeqIO.write(new_record, f_out, 'fasta')
    """
}


process cdhit {
    label 'general'
    input:
        tuple val(label), path(fasta)
    output:
        tuple val(label), path("${label}.fasta")
    cpus 4
    """
    cd-hit-est -i ${fasta} -o ${label}.fasta -c ${params.cdhit_threshold} -T ${task.cpus}
    """
}

process add_unknown_segment_info {
    label 'general'
    input:
        path('in.nogroup.fasta')
    output:
        path('out.nogroup.fasta')
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    with open('out.nogroup.fasta', 'w') as f_out:
        for seq in SeqIO.parse('in.nogroup.fasta', 'fasta'):
            new_record = SeqRecord(
                id=seq.id,
                name=seq.name,
                description=seq.description + '|Unknown|Unknown',
                seq=seq.seq,
            )
            SeqIO.write(new_record, f_out, 'fasta')
    """
}

process build_blast_db {
    label 'general'
    input:
        path(all_fasta)
    output:
        path('blast_db')
    """
    mkdir -p blast_db

    makeblastdb \
        -dbtype nucl \
        -in ${all_fasta} \
        -out blast_db/ALL \
        -parse_seqids \
        -blastdb_version 5
    """
}


process write_output_config {
    label 'general'
    input:
        path('input.yaml')
    output:
        path('output.yaml')
    """
    #!/usr/bin/env python
    import yaml
    with open('input.yaml') as f_in:
        config = yaml.safe_load(f_in)
    out = {
        'all': {
            'blast_db': 'blast_db',
            'blast_prefix': 'ALL',
            'fasta': 'ALL.fasta',
        },
        'curated': {
            label: {
                'fasta': f'{label}.fasta',
                'name': conf['name'],
                'organisms': sorted(set([
                    org
                    for tax in conf['taxa']
                    for org in tax['organisms']
                ])),
                'segments': list(conf['segments'].keys()),
            }
            for label, conf in config['curated'].items()
        },
        'filters': {
            label: f'{label}.fasta'
            for label in config['filters']
        },
        'version': '${params.output_version}',
        'description': '${params.output_description}',
        'params.fasta_sequences': '${params.fasta_sequences}',
        'params.taxa_config': '${params.taxa_config}',
        'params.filter_max_n_share': ${params.filter_max_n_share},
        'params.filter_min_relative_length': ${params.filter_min_relative_length},
        'params.cdhit_threshold': ${params.cdhit_threshold},
    }
    with open('output.yaml', 'w') as f_out:
        yaml.dump(dict(out), f_out, default_flow_style=False, sort_keys=False)
    """
}


process output {
    label 'general'
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
