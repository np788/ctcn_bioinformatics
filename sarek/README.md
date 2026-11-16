# Sarek joint germline workflow example

This has an example of a WGS workflow using **nf-core/sarek v3.10.0** under **Nextflow** with
**Singularity** and **SLURM** to perform GATK HaplotypeCaller germline variant
calling and joint genotyping from CRAM/CRAI files.

This created gVCFs, performed joint germline calling, and produced a recalibrated VCF used in downstream donor matching.

Launcher:

```bash
nextflow run nf-core/sarek \
    -r 3.10.0 \
    -profile singularity \
    -c config/slurm.config \
    -params-file config/sarek_joint_params.yaml \
    -work-dir work/joint_full
```