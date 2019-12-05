# get the areas of interest from bam file into multiple bams based on conditions
# get mates somehow and make into paired end fastq

sample=$1
threads=$2

# get index
#samtools index -@ $threads "$sample".bam

# get unmapped reads
samtools view -b -h -f 4 -@ $threads "$sample".bam > "$sample".unmapped.bam

# get TRA reads
samtools view -b -h -@ $threads "$sample".bam "chr14:21621904-22552132" > "$sample".TRA.bam

# get TRB reads
samtools view -b -h -@ $threads "$sample".bam "chr7:142299011-142813287" > "$sample".TRB.bam

# get TRG reads
samtools view -b -h -@ $threads "$sample".bam "chr7:38240024-38368055" > "$sample".TRG.bam

# get TRD reads
samtools view -b -h -@ $threads "$sample".bam "chr14:22422546-22466577" > "$sample".TRD.bam

# merge TR reads
samtools merge -f -@ $threads "$sample".TR.bam "$sample".TRA.bam "$sample".TRB.bam "$sample".TRG.bam "$sample".TRD.bam

# sort TR reads
samtools sort -@ $threads -n -o "$sample".TR.sorted.bam "$sample".TR.bam

# sort unmapped reads
samtools sort -n -o "$sample".unmapped.sorted.bam "$sample".unmapped.bam

# get unmapped fastq
bedtools bamtofastq -i "$sample".unmapped.sorted.bam -fq "$sample".unmapped.R1.fastq -fq2 "$sample".unmapped.R2.fastq

# get TR fastq
bedtools bamtofastq -i "$sample".TR.sorted.bam -fq "$sample".TR.R1.fastq -fq2 "$sample".TR.R2.fastq

# get combined fastq
cat "$sample".unmapped.R1.fastq "$sample".TR.R1.fastq > "$sample".combined.T.R1.fastq
cat "$sample".unmapped.R2.fastq "$sample".TR.R2.fastq > "$sample".combined.T.R2.fastq

# mixcr align
mixcr align -f -r "$sample".T.report -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$sample".combined.T.R1.fastq "$sample".combined.T.R2.fastq "$sample".T.vdjca

# mixcr assemblePartial
mixcr assemblePartial -f "$sample".T.vdjca "$sample".T.rescued.vdjca
mixcr assemblePartial -f "$sample".T.rescued.vdjca "$sample".T.rescued2.vdjca

# mixcr assemble
mixcr assemble -f "$sample".T.rescued2.vdjca "$sample".T.clns

# mixcr exportClones
mixcr exportClones -f -o -t "$sample".T.clns "$sample".T.clones.tsv
