#!/usr/bin/python3
from pybedtools import BedTool
import subprocess as sp

# Create IGV batch script to create screenshots at given position in bed

# For snoscan predictions in S. pombe
snoscan_bed = BedTool(snakemake.input.snoscan_bed)
sp.call('mkdir -p results/figures/screenshots/snoscan/S_pombe/', shell=True)
snoscan_bed.igv(name=True, path="./results/figures/screenshots/snoscan/S_pombe/", img="png")
sp.call(f'mv {snoscan_bed.igv_script} {snakemake.output.batch_script_snoscan}', shell=True)

# For snoreport predictions in S. pombe
snoreport_bed = BedTool(snakemake.input.snoreport_bed)
sp.call('mkdir -p results/figures/screenshots/snoreport/S_pombe/', shell=True)
snoreport_bed.igv(name=True, path="./results/figures/screenshots/snoreport/S_pombe/", img="png")
sp.call(f'mv {snoreport_bed.igv_script} {snakemake.output.batch_script_snoreport}', shell=True)

# For infernal predictions in S. pombe
infernal_bed = BedTool(snakemake.input.infernal_bed)
sp.call('mkdir -p results/figures/screenshots/infernal_rfam/S_pombe/', shell=True)
infernal_bed.igv(name=True, path="./results/figures/screenshots/infernal_rfam/S_pombe/", img="png")
sp.call(f'mv {infernal_bed.igv_script} {snakemake.output.batch_script_infernal}', shell=True)

# For snoBIRD predictions in S. pombe
snoBIRD_bed = BedTool(snakemake.input.snoBIRD_bed)
sp.call('mkdir -p results/figures/screenshots/snoBIRD/S_pombe/', shell=True)
snoBIRD_bed.igv(name=True, path="./results/figures/screenshots/snoBIRD/S_pombe/", img="png")
sp.call(f'mv {snoBIRD_bed.igv_script} {snakemake.output.batch_script_snoBIRD}', shell=True)


