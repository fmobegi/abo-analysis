source activate minionvironment

#python src/AnalyzeAbo_Main.py  \
# --reference="/home/ben/Github/abo-analysis/input/A1_reference_Exon7.fasta" \
# --alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
# --output="/home/ben/Github/abo-analysis/results/Exon7PolymorphismResults" \
# --analysis-type="ALLELES"
 
#python src/AnalyzeAbo_Main.py  \
# --reference="/home/ben/Github/abo-analysis/input/A1_reference_Exon6.fasta" \
#--alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
# --output="/home/ben/Github/abo-analysis/results/Exon6PolymorphismResults" \
# --analysis-type="ALLELES"

cd /home/ben/Github/abo-analysis

python src/AnalyzeAbo_Main.py  \
 --reference="/home/ben/ben_share/2017.Nov.8.ABO.PatientStudy/reads_bc51/A1_01_01_1_reference_Exon6.fasta" \
 --alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
 --output="/home/ben/ben_share/2017.Nov.8.ABO.PatientStudy/reads_bc51/bc51_exon6" \
 --analysis-type="READS" \
 --reads="/home/ben/ben_share/2017.Nov.8.ABO.PatientStudy/reads_bc51/reads_BC51_Pass.fastq" \
 
 python src/AnalyzeAbo_Main.py  \
 --reference="/home/ben/ben_share/2017.Nov.8.ABO.PatientStudy/reads_bc51/A1_01_01_1_reference_Exon7.fasta" \
 --alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
 --output="/home/ben/ben_share/2017.Nov.8.ABO.PatientStudy/reads_bc51/bc51_exon7" \
 --analysis-type="READS" \
 --reads="/home/ben/ben_share/2017.Nov.8.ABO.PatientStudy/reads_bc51/reads_BC51_Pass.fastq"
 
source deactivate

