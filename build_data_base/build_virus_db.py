"""
Merge curated and non-curated content to build a complete database
"""

# TODO: gzip fasta files?

import yaml
import argparse
import pathlib
import subprocess
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def run_sub(cmd):
    print(' '.join(cmd))
    try:
        subprocess.run(cmd, shell=False, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Command '{e.cmd}' failed with return code {e.returncode}")
        print("Standard Output:\n", e.stdout)
        print("Standard Error:\n", e.stderr)


def cat(fnames, outfile):
    with open(outfile, 'w') as f_out:
        for fn in fnames:
            with open(fn) as f_in:
                for line in f_in:
                    f_out.write(line)


def pstr(path):
    return str(path.resolve())


def minimap2_index(path_fasta, outpath):
    path_index = outpath / f'{path_fasta.stem}.mmi'
    run_sub(['minimap2', '-d', pstr(path_index), pstr(path_fasta)])
    return path_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'config',
        help='Config file in yaml format.'
    )
    parser.add_argument(
        'output_version',
        help='Version to write to the output config.'
    )
    args = parser.parse_args()

    with open(args.config) as f_config:
        config = yaml.safe_load(f_config)

    outpath = pathlib.Path(config['output_directory'])
    outpath.mkdir(exist_ok=True, parents=True)

    print(list(config['curated'].values())[0])

    # Merge curated and non-curated species.
    curated_organisms = set(
        organism_name.lower().strip()
        for curated_entry in config['curated'].values()
        for organism_name in curated_entry['organisms']
    )
    all_fasta = outpath / 'ALL.fasta'
    with all_fasta.open(mode='w') as f_all_fasta:
        # first write all curated sequences to the collection
        for curated_entry in config['curated'].values():
            for segment_info in curated_entry['segments']:
                with open(segment_info['file']) as f_segment:
                    for line in f_segment:
                        f_all_fasta.write(line)
        # now add sequences for all species, that are not included in the
        # curated data sets
        for record in SeqIO.parse(config['rest']['fasta'], 'fasta'):
            organism = record.description.rsplit('|', 1)[1].lower().strip()
            if organism not in curated_organisms:
                SeqIO.write(
                    SeqRecord(
                        seq=record.seq,
                        id=record.id,
                        name=record.name,
                        description=record.description.strip() + '|forward|Unknown'
                    ),
                    f_all_fasta,
                    'fasta'
                )

    out_config = {}

    # create the blast index
    subdir_blast = 'blast_db'
    blast_db_path = outpath / subdir_blast
    blast_db_path.mkdir(exist_ok=True, parents=True)
    blast_db_prefix = pstr(blast_db_path / 'ALL')
    run_sub([
        'makeblastdb',
        '-dbtype', 'nucl',
        '-in', pstr(all_fasta),
        '-out', blast_db_prefix,
        '-parse_seqids',
        '-blastdb_version', '5'
    ])

    #fname_all_filter = minimap2_index(all_fasta, outpath)
    out_config['all'] = {
        'fasta': all_fasta.name,
        # 'minimap2_index': pstr(fname_all_filter),
        'blast_db': subdir_blast,
    }

    out_config['curated'] = {}
    for dataset_name, dataset_info in config['curated'].items():
        fasta_files = [segment_info['file'] for segment_info in dataset_info['segments']]
        path_fasta_concat = outpath / f'{dataset_name}.fasta'
        cat(fasta_files, path_fasta_concat)
        # fname_mmi = minimap2_index(path_fasta_concat, outpath)
        msa_paths = {}
        for segment_info in dataset_info['segments']:
            if segment_info['msa'] is not None:
                segment_name = segment_info['segment']
                msa_path = outpath / f'{dataset_name}.{segment_name}.msa.fasta'
                msa_paths[segment_info['segment']] = msa_path.name
                shutil.copyfile(segment_info['msa'], msa_path)
            else:
                msa_paths[segment_info['segment']] = None
        out_config['curated'][dataset_name] = {
            'name': dataset_info['name'],
            'organisms': dataset_info['organisms'],
            'segments': sorted(msa_paths.keys()),
            'msa': msa_paths,
            'fasta': path_fasta_concat.name,
            # 'minimap2_index': pstr(fname_mmi),
        }

    out_config['filters'] = {}
    for filter_name, fname_filter_in in config['filters'].items():
        filter_path = outpath / f'{filter_name}.fasta'
        out_config['filters'][filter_name] = filter_path.name
        shutil.copyfile(fname_filter_in, filter_path)

    out_config['version'] = args.output_version

    path_config_out = outpath / 'db.yaml'
    with path_config_out.open(mode='w') as f_config_out:
        yaml.dump(out_config, f_config_out, default_flow_style=False)


if __name__ == '__main__':
    main()
