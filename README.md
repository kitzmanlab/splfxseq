# splfxseq

Scripts and python modules to design and process data from massively parallel splicing assay experiments.

### installation

It is recommended to install this package in a dedicated miniconda or miniforge environment.  After checking out from github, navigate to the spflxseq folder and run:

```
pip install .
```

If you want an editable install (ie, where source files can be edited so that changes take effect 'in place', add the -e flag)

```
pip install -e .
```

Requirements:
 - commonly used python modules including: numpy, scipy, pysam, pandas, cutadapt
 - snakemake is required to run the included pipelines
 - these pipelines require certain other software be installed/available in your enviroment: STAR, samtools, starcode, multiqc, fastqc

