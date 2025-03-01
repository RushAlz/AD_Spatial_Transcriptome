# cell2location_smk

Snakemake pipeline for running cell2location

First copy the input files to the bucket specified at `default-remote-prefix`

These files are required (path can be changed in the Snakefile):

```
input_path = "input/split_by_section/"
nuclei_file = "input/split_by_section/avgNuclei_perSample.txt"
reference_file = "resources/references/Tsai_reference/cellsubtypeSignatures.csv"
```

Run snakemake. Remember to use **screen**

```
snakemake -s Snakefile \
--google-lifesciences --default-remote-prefix rvialle/projects/ST/cell2location \
--use-conda --google-lifesciences-region us-central1 -j8 \
--container-image ricardovialle/cell2location-docker -npr
```

Issue regarding the *gpu* with *google lifesciences*:

https://github.com/snakemake/snakemake/issues/457
