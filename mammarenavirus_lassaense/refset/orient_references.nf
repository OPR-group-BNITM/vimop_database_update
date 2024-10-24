nextflow.enable.dsl = 2


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


process minimap {
    input:
        tuple val(alias), path('refs.fasta'), path('queries.fasta')
    output:
        tuple val(alias), path('mapped.sam')
    """
    minimap2 --secondary=no --eqx -x map-ont -a refs.fasta queries.fasta > mapped.sam
    """
}


process get_orientations {
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

workflow {

    // WARNING: at the moment we can only have one msa remapping parameter set,
    // as the msa remapping parameters are not written to the results

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