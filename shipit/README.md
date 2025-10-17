# Assembly of the downloadable content

Each part of the data base (virus, centrifuge and conaminants) is packaged separately and has its own version.
Finally, a common config is created to define one data base version.

## Compress the data bases

Run `./shipit/archive_and_split.sh /path/to/db virus|centrifuge|contaminants output_dir`.
Replace `path/to/db` with the directory of you new data base, select the type of data base (virus, centrifuge or contaminants) and give an output directory.

## Assemble configs for a new data base version

To create a whole data base, we need a download config combining the information about the three sub data bases.
To merge them run `python merge_config.py --virus path/to/download/configs/virus.v1.0.files.yaml --contaminants path/to/download/configs/virus.v1.0.files.yaml --centrifuge path/to/download/configs/centrifuge.v1.0.files.yaml --output vimop_db.1.0.yaml`.
Replace the version numbers with the correct version numbers.
