# The directory with raw data and genome files
cd data_drop

# Generate genome indexes
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir HS38_STAR --genomeFastaFiles GRCh38.p13.genome.fa --sjdbGTFfile gencode.v40.annotation.gtf --sjdbOverhang 149
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir MM10_STAR --genomeFastaFiles GRCm39.genome.fa --sjdbGTFfile gencode.vM29.annotation.gtf --sjdbOverhang 149

# Filter and align
fastp --detect_adapter_for_pe -i scWT_02_R1.fastq.gz -I scWT_02_R2.fastq.gz -o scWT_02_R1.filter.fastq.gz -O scWT_02_R2.filter.fastq.gz
python3 sub-script/BCFilter.py
STAR --runThreadN 16 --genomeDir HS38_STAR --readFilesIn scWT_02_R1.filter.fastq.gz scWT_02_R2.filter.fastq.gz scWT_02_I2.filter.fastq.gz --soloType CB_samTagOut --soloCBwhitelist 737K-cratac-v1_rev.txt --soloCBmatchWLtype 1MM --outFileNamePrefix scWT_02_hs --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes CR CB XS --outSAMstrandField intronMotif --soloBarcodeReadLength 0
java -jar picard.jar MarkDuplicates I=scWT_02_hsAligned.sortedByCoord.out.bam O=scWT_02_hs.bam CREATE_INDEX=false REMOVE_DUPLICATES=false M=scWT_02_hs.dup.txt BARCODE_TAG=CB
rm scWT_02_hsAligned.sortedByCoord.out.bam
samtools index scWT_02_hs.bam
qualimap rnaseq -outdir scWT_02_hs_qc -bam scWT_02_hs.bam -gtf gencode.v40.annotation.gtf --java-mem-size=16G
STAR --runThreadN 16 --genomeDir MM10_STAR --readFilesIn scWT_02_R1.filter.fastq.gz scWT_02_R2.filter.fastq.gz scWT_02_I2.filter.fastq.gz --soloType CB_samTagOut --soloCBwhitelist 737K-cratac-v1_rev.txt --soloCBmatchWLtype 1MM --outFileNamePrefix scWT_02_mm --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes CR CB XS --outSAMstrandField intronMotif --soloBarcodeReadLength 0
java -jar picard.jar MarkDuplicates I=scWT_02_mmAligned.sortedByCoord.out.bam O=scWT_02_mm.bam CREATE_INDEX=false REMOVE_DUPLICATES=false M=scWT_02_mm.dup.txt BARCODE_TAG=CB
rm scWT_02_mmAligned.sortedByCoord.out.bam
samtools index scWT_02_mm.bam
qualimap rnaseq -outdir scWT_02_mm_qc -bam scWT_02_mm.bam -gtf gencode.vM29.annotation.gtf --java-mem-size=16G
fastp --detect_adapter_for_pe -i SRR19619979_1.fastq -I SRR19619979_2.fastq -o SRR19619979_1.filter.fastq -O SRR19619979_2.filter.fastq
STAR --runThreadN 16 --genomeDir /md01/xuyc7/library/HS38_150 --readFilesIn SRR19619979_1.filter.fastq SRR19619979_2.filter.fastq --outFileNamePrefix SRR19619979 --outSAMtype BAM SortedByCoordinate
java -jar ~/tools/bin/picard.jar MarkDuplicates I=SRR19619979Aligned.sortedByCoord.out.bam O=SRR19619979.bam CREATE_INDEX=false REMOVE_DUPLICATES=false M=SRR19619979.dup.txt
samtools index SRR19619979.bam
qualimap rnaseq -outdir SRR19619979_qc -bam SRR19619979.bam -gtf ~/library/gencode.v40.annotation.gtf --java-mem-size=16G

# Annotate and count
python3 sub-script/calcHybrid.py
python3 sub-script/calcReads.py
Rscript sub-script/calcCells.R
Rscript sub-script/calcCells_other.R

# Calculate coverage and transcripts
python3 sub-script/calcCov.py
stringtie scWT_02_hs_filter.bam -o scWT_02_hs_filter.gtf -p 16 -G gencode.v40.annotation.gtf -l HSS
stringtie SRR19619979.bam -o 293T_bulk.gtf -p 16 -G gencode.v40.annotation.gtf -l HSB
Rscript sub-script/calcTrans.R
python3 sub-script/calcTrans.py

# Collect the results
cp scWT_02.rds ../
cp scWT_feature_info.tsv ../
cp cell_info.csv ../
cp scWT_02_hs_filter_coverage.csv ../
cp scWT_02_mm_filter_coverage.csv ../
cp or_02.rds ../
cp trans_info_all.txt ../
cp *_trans_info.txt ../
cp res_info_density.txt ../
cp res_info_coord.txt ../
cp res_info_junction.txt ../
cp scWT_02_trans_info_filter.h5 ../

