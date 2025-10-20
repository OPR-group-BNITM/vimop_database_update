# ViMOP DB centrifuge index creation

This repository contains scripts and documentation to update the vimop database centrifuge index.

To build a new index build the virus data base first (see `virus` directory of this repository).
Then run

```
./slurmrun_rs_setup.sh
```

to get refseq sequences for bacteria, archaea and eukaryotes (only human and mouse).
Subsequently run

```
./slurmrun_assign_virus_taxids.sh
```

to assign taxIDs to all virus sequences.
The virus sequences come from the previously built virus data base.
Hence, make sure to update the version and path in `env.sh`.
Finally run

```
./slurmrun_build_centrifuge_index.sh
```