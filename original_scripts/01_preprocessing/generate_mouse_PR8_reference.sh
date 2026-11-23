module purge
cd /data/scratch/projects/punim1236/Reference/pr8_genome
module load GCCcore/11.3.0 gffread/0.12.7
echo "modules loaded successfully"
#gffread GCF_000865725.1_ViralMultiSegProj15521_genomic.gff --force-exons --gene2exon -T -o pr8_genes.gtf
# commented out as this worked on the command line but not inside the script for some reason

conda_load
conda activate star_2.6.1

cd /data/scratch/projects/punim1236/Reference
mkdir mm10_pr8
cat mouse_mm10*/fasta/genome.fa pr8_genome/*genomic.fna >> mm10_pr8/genome.fa
cat mouse_mm10*/genes/genes.gtf pr8_genome/pr8_genes.gtf >> mm10_pr8/genes.gtf

cd /data/scratch/projects/punim1236/curioseeker-2.0.0/generate_seeker_reference
bash generate_seeker_reference.sh \
  /data/scratch/projects/punim1236/Reference/mm10_pr8/genome.fa \
  /data/scratch/projects/punim1236/Reference/mm10_pr8/genes.gtf \
  chrM \
  mm10_opt_PR8
