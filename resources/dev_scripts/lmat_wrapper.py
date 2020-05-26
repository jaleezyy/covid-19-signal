# This script is no longer used!

# It runs LMAT-1.2.6 from Fin's Docker container, and was previously invoked
# by Snakemake, in the 'run_lmat' rule.  However, we've now switched to running
# LMAT via Fin's conda recipe, which is cleaner than using Docker.

# TODO: write comment explaining why this klunky script is necessary,
# instead of just doing 'docker run...' in the Snakefile!

import os
import sys
import shutil

print('lmat_wrapper.py: start')

# The following parameters are passed in from 'snakemake'
lmat_basedir = snakemake.params.lmat_basedir
lmat_db = snakemake.params.lmat_db
logfile = snakemake.params.logfile
outdir = snakemake.params.outdir

if os.path.exists(logfile):
    print(f'lmat_wrapper.py: deleting logfile {logfile}')
    os.remove(logfile)

docker_cmd = (
    f'docker run --rm'
    f' -v {lmat_basedir}/data:/data'
    f' -v {lmat_basedir}/runtime_inputs:/runtime_inputs'
    f' -v {os.getcwd()}/{outdir}:/pipeline'
    f' finlaymaguire/lmat:1.2.6'
    f' bash /bin/run_rl.sh'
    f' --db_file=/data/{lmat_db}'
    f' --query_file=/pipeline/assembly.tiled.fasta'
    f' --odir=/pipeline'
    f' --overwrite --verbose --threads={snakemake.threads}'
)

print(f'lmat_wrapper.py: {docker_cmd}')
ret = os.system(docker_cmd)

for f in os.listdir(outdir):
    src = os.path.join(f'{outdir}/{f}')
    dst = os.path.join(f'{outdir}/tmp-{f}')
    
    if os.stat(src).st_uid == 0:
        print(f'lmat_wrapper.py: copy {src} -> {dst}')
        shutil.copyfile(src, dst)

        print(f'lmat_wrapper.py: rename {dst} -> {src}')
        os.rename(dst, src)

if ret != 0:
    print(f"lmat_wrapper.py: dockerized LMAT run failed (nonzero exit status)")
    sys.exit(1)
    
if not os.path.exists(logfile):
    print(f"lmat_wrapper.py: dockerized LMAT run failed (didn't write logfile {logfile})")
    sys.exit(1)

with open(logfile) as f:
    for line in f:
        if line.find("Exit status: 0"):
            print('lmat_wrapper.py: success')
            sys.exit(0)

print(f"lmat_wrapper.py: dockerized LMAT run failed (logfile {logfile} didn't contain 'Exit status: 0'")
sys.exit(1)
            
