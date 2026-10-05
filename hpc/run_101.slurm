
cd ~/u-dales-run/101
nano namoptions.101
cp ~/u-dales-run/101/namoptions.101 $SCRATCH_FLASH/u-dales-run/101/ 
diff -q ~/u-dales-run/101/namoptions.101 $SCRATCH_FLASH/u-dales-run/101/namoptions.101  % check differents
OUTDIR=$SCRATCH_FLASH/output/101_3600s
mkdir -p $SCRATCH_FLASH/output
cd ~/u-dales-run/101
nano run_101.slurm
===========================================================================================================
# ---- Software environment ----
module load spack/v0.23
module load gcc/12.4.0
module load openmpi/5.0.7_gcc12
module load cmake/3.26
module load fftw/3.3.10
spack load netcdf-fortran@4.6.1

# ---- OpenMP ----
export OMP_NUM_THREADS=1

# ---- POLITO OpenMPI/hwloc workaround ----
export HWLOC_LIBXML=0

# ---- Paths ----
INPUTDIR=$HOME/u-dales-run/101
OUTDIR=$SCRATCH_FLASH/output/101_3600s
DA_BUILD=/home/jili/u-dales/build/release/u-dales

# ---- Prepare output directory ----
mkdir -p "$OUTDIR"

# ---- Copy case inputs, as official hpc_execute.sh does ----
cp -r "$INPUTDIR"/* "$OUTDIR"/

# ---- Enter simulation directory ----
cd "$OUTDIR"

echo "============================================================"
echo "uDALES Case 101"
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Working directory: $PWD"
echo "Executable: $DA_BUILD"
echo "MPI ranks: 4"
echo "============================================================"

# ---- Run uDALES ----
mpirun --bind-to none --map-by slot -np 4 \
    "$DA_BUILD" namoptions.101

echo "============================================================"
echo "uDALES finished"
echo "End time: $(date)"
echo "============================================================"
============================================================================================================
sbatch run_101.slurm %run
squeue -j （1977745）
squeue -j 1975002 % check process
ls -lh $SCRATCH_FLASH/output/101_3600s
