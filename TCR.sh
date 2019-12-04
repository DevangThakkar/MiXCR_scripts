# get the areas of interest from bam file into multiple bams based on conditions
# get mates somehow and make into paired end fastq

sample=$1

bam_dir=$2

result_dir="$3"/"$1"
mkdir -p $result_dir

threads=$4

# get index
#samtools index -@ $threads "$sample".bam

# get unmapped reads
#samtools view -b -h -f 4 -@ $threads "$bam_dir"/"$sample".bam > "$result_dir"/"$sample".unmapped.bam

# get TRA reads
samtools view -b -h -@ $threads "$bam_dir"/"$sample".bam "chr14:21621904-22552132" > "$result_dir"/"$sample".TRA.bam

# get TRB reads
samtools view -b -h -@ $threads "$bam_dir"/"$sample".bam "chr7:142299011-142813287" > "$result_dir"/"$sample".TRB.bam

# get TRG reads
samtools view -b -h -@ $threads "$bam_dir"/"$sample".bam "chr7:38240024-38368055" > "$result_dir"/"$sample".TRG.bam

# get TRD reads
samtools view -b -h -@ $threads "$bam_dir"/"$sample".bam "chr14:22422546-22466577" > "$result_dir"/"$sample".TRD.bam

# merge TR reads
samtools merge -f -@ $threads "$result_dir"/"$sample".TR.bam "$result_dir"/"$sample".TRA.bam "$result_dir"/"$sample".TRB.bam "$result_dir"/"$sample".TRG.bam "$result_dir"/"$sample".TRD.bam

# sort TR reads
samtools sort -@ $threads -n -o "$result_dir"/"$sample".TR.sorted.bam "$result_dir"/"$sample".TR.bam

# sort unmapped reads
#samtools sort -n -o "$result_dir"/"$sample".unmapped.sorted.bam "$result_dir"/"$sample".unmapped.bam

# get unmapped fastq
#bedtools bamtofastq -i "$result_dir"/"$sample".unmapped.sorted.bam -fq "$result_dir"/"$sample".unmapped.R1.fastq -fq2 "$result_dir"/"$sample".unmapped.R2.fastq

# get TR fastq
bedtools bamtofastq -i "$result_dir"/"$sample".TR.sorted.bam -fq "$result_dir"/"$sample".TR.R1.fastq -fq2 "$result_dir"/"$sample".TR.R2.fastq

# get combined fastq
cat "$result_dir"/"$sample".unmapped.R1.fastq "$result_dir"/"$sample".TR.R1.fastq > "$result_dir"/"$sample".combined.T.R1.fastq
cat "$result_dir"/"$sample".unmapped.R2.fastq "$result_dir"/"$sample".TR.R2.fastq > "$result_dir"/"$sample".combined.T.R2.fastq

# mixcr align
mixcr align -f -r "$result_dir"/"$sample".report -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$result_dir"/"$sample".combined.T.R1.fastq "$result_dir"/"$sample".combined.T.R2.fastq "$result_dir"/"$sample".T.vdjca

# mixcr assemblePartial
mixcr assemblePartial -f "$result_dir"/"$sample".T.vdjca "$result_dir"/"$sample".T.rescued.vdjca
mixcr assemblePartial -f "$result_dir"/"$sample".T.rescued.vdjca "$result_dir"/"$sample".T.rescued2.vdjca

# mixcr assemble
mixcr assemble -f "$result_dir"/"$sample".T.rescued2.vdjca "$result_dir"/"$sample".T.clns

# mixcr exportClones
mixcr exportClones -f -o -t "$result_dir"/"$sample".T.clns "$result_dir"/"$sample".T.clones.tsv

