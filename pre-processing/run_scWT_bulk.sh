# The directory with raw data and genome files
cd data_bulk

# Generate genome indexes
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir HS38_STAR --genomeFastaFiles GRCh38.p13.genome.fa --sjdbGTFfile gencode.v40.annotation.gtf --sjdbOverhang 149

# scWT(40%Met)
fastp --detect_adapter_for_pe -i met40.R1.fastq -I met40.R2.fastq -o met40.R1.filter.fastq -O met40.R2.filter.fastq
STAR --runThreadN 16 --genomeDir HS38_STAR --readFilesIn met40.R1.filter.fastq met40.R2.filter.fastq --outFileNamePrefix met40 --outSAMtype BAM SortedByCoordinate
java -jar ~/tools/bin/picard.jar MarkDuplicates I=met40Aligned.sortedByCoord.out.bam O=met40.bam CREATE_INDEX=false REMOVE_DUPLICATES=false M=met40.dup.txt
samtools index met40.bam
qualimap rnaseq -outdir met40_qc -bam met40.bam -gtf gencode.v40.annotation.gtf --java-mem-size=16G
htseq-count -r pos -f bam met40.bam gencode.v40.annotation.gtf > met40.count.txt

# scWT(60%Met)
fastp --detect_adapter_for_pe -i met60.R1.fastq -I met60.R2.fastq -o met60.R1.filter.fastq -O met60.R2.filter.fastq
STAR --runThreadN 16 --genomeDir HS38_STAR --readFilesIn met60.R1.filter.fastq met60.R2.filter.fastq --outFileNamePrefix met60 --outSAMtype BAM SortedByCoordinate
java -jar ~/tools/bin/picard.jar MarkDuplicates I=met60Aligned.sortedByCoord.out.bam O=met60.bam CREATE_INDEX=false REMOVE_DUPLICATES=false M=met60.dup.txt
samtools index met60.bam
qualimap rnaseq -outdir met60_qc -bam met60.bam -gtf gencode.v40.annotation.gtf --java-mem-size=16G
htseq-count -r pos -f bam met60.bam gencode.v40.annotation.gtf > met60.count.txt

# scWT(80%Met)
fastp --detect_adapter_for_pe -i met80.R1.fastq -I met80.R2.fastq -o met80.R1.filter.fastq -O met80.R2.filter.fastq
STAR --runThreadN 16 --genomeDir HS38_STAR --readFilesIn met80.R1.filter.fastq met80.R2.filter.fastq --outFileNamePrefix met80 --outSAMtype BAM SortedByCoordinate
java -jar ~/tools/bin/picard.jar MarkDuplicates I=met80Aligned.sortedByCoord.out.bam O=met80.bam CREATE_INDEX=false REMOVE_DUPLICATES=false M=met80.dup.txt
samtools index met80.bam
qualimap rnaseq -outdir met80_qc -bam met80.bam -gtf gencode.v40.annotation.gtf --java-mem-size=16G
htseq-count -r pos -f bam met80.bam gencode.v40.annotation.gtf > met80.count.txt

# RNA-seq(std.)
fastp --detect_adapter_for_pe -i SRR307005_1.fastq -I SRR307005_2.fastq -o SRR307005.R1.filter.fastq -O SRR307005.R2.filter.fastq
STAR --runThreadN 16 --genomeDir HS38_STAR --readFilesIn SRR307005.R1.filter.fastq SRR307005.R2.filter.fastq --outFileNamePrefix SRR307005 --outSAMtype BAM SortedByCoordinate
java -jar ~/tools/bin/picard.jar MarkDuplicates I=met80Aligned.sortedByCoord.out.bam O=SRR307005.bam CREATE_INDEX=false REMOVE_DUPLICATES=false M=SRR307005.dup.txt
samtools index SRR307005.bam
qualimap rnaseq -outdir SRR307005_qc -bam SRR307005.bam -gtf gencode.v40.annotation.gtf --java-mem-size=16G
htseq-count -r pos -f bam SRR307005.bam gencode.v40.annotation.gtf > SRR307005.count.txt

# Collect the results
mkdir ../bulk_met
cp met40.count.txt  ../bulk_met/met40.cov.txt
cp "met40_qc/raw_data_qualimapReport/rnaseq_qc_results.txt" ../bulk_met/met40.stat.txt
cp "met40_qc/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" ../bulk_met/met40.cov.txt
cp met60.count.txt  ../bulk_met/met60.cov.txt
cp "met60_qc/raw_data_qualimapReport/rnaseq_qc_results.txt" ../bulk_met/met60.stat.txt
cp "met60_qc/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" ../bulk_met/met60.cov.txt
cp met80.count.txt  ../bulk_met/met80.cov.txt
cp "met80_qc/raw_data_qualimapReport/rnaseq_qc_results.txt" ../bulk_met/met80.stat.txt
cp "met80_qc/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" ../bulk_met/met80.cov.txt
cp SRR307005.count.txt  ../bulk_met/norm.cov.txt
cp "SRR307005_qc/raw_data_qualimapReport/rnaseq_qc_results.txt" ../bulk_met/norm.stat.txt
cp "SRR307005_qc/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" ../bulk_met/norm.cov.txt

