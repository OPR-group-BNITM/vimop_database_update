#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from typing import Any, Dict
import yaml


def load_yaml(path: str) -> Any:
    with open(path, 'r', encoding='utf-8') as fh:
        return yaml.safe_load(fh)


def build_manifest(version: str,
                   description: str,
                   virus_data: Any,
                   contaminants_data: Any,
                   centrifuge_data: Any) -> Dict[str, Any]:
    return {
        'version': version,
        'description': description,
        'sub_databases': {
            'virus': virus_data,
            'contaminants': contaminants_data,
            'centrifuge': centrifuge_data,
        }
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog='build_vimop_manifest',
        description='Combine virus, contaminants, and centrifuge YAMLs into a ViMOP bucket manifest.'
    )
    p.add_argument('--virus', required=True, help='Path to virus YAML file.')
    p.add_argument('--contaminants', required=True, help='Path to contaminants YAML file.')
    p.add_argument('--centrifuge', required=True, help='Path to centrifuge YAML file.')
    p.add_argument('--version', required=True, help='Manifest version (e.g., 1.1).')
    p.add_argument('--description', default='ViMOP_Bucket',
                   help='Manifest description (default: ViMOP_Bucket).')
    p.add_argument('--output', '-o', default='-',
                   help='Output path for combined YAML (default: stdout).')
    return p.parse_args()


def main() -> None:
    args = parse_args()

    try:
        virus_data = load_yaml(args.virus)
    except Exception as e:
        sys.stderr.write(f'Failed to read --virus YAML: {e}\n')
        sys.exit(2)

    try:
        contaminants_data = load_yaml(args.contaminants)
    except Exception as e:
        sys.stderr.write(f'Failed to read --contaminants YAML: {e}\n')
        sys.exit(3)

    try:
        centrifuge_data = load_yaml(args.centrifuge)
    except Exception as e:
        sys.stderr.write(f'Failed to read --centrifuge YAML: {e}\n')
        sys.exit(4)

    manifest = build_manifest(
        version=args.version,
        description=args.description,
        virus_data=virus_data,
        contaminants_data=contaminants_data,
        centrifuge_data=centrifuge_data
    )

    try:
        dumped = yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True)
    except Exception as e:
        sys.stderr.write(f'Failed to serialize combined YAML: {e}\n')
        sys.exit(5)

    if args.output == '-' or args.output.lower() == 'stdout':
        sys.stdout.write(dumped)
    else:
        try:
            with open(args.output, 'w', encoding='utf-8') as fh:
                fh.write(dumped)
        except Exception as e:
            sys.stderr.write(f'Failed to write output file: {e}\n')
            sys.exit(6)


if __name__ == '__main__':
    main()
