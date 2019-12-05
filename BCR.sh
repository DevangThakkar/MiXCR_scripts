# get the areas of interest from bam file into multiple bams based on conditions
# get mates somehow and make into paired end fastq

sample=$1
threads=$2

report=$3
output=$4

# get index
#samtools index -@ $threads "$sample".bam

# get unmapped reads
samtools view -b -h -f 4 -@ $threads "$sample".bam > "$sample".unmapped.bam

# get IGH reads
samtools view -b -h -@ $threads "$sample".bam "chr14:105586437-106879844" > "$sample".IGH.bam

# get IGK reads
samtools view -b -h -@ $threads "$sample".bam "chr2:88857361-90235368" > "$sample".IGK.bam

# get IGL reads
samtools view -b -h -@ $threads "$sample".bam "chr22:22026076-22922913" > "$sample".IGL.bam

# merge IG reads
samtools merge -f -@ $threads "$sample".IG.bam "$sample".IGH.bam "$sample".IGK.bam "$sample".IGL.bam

# sort IG reads
samtools sort -@ $threads -n -o "$sample".IG.sorted.bam "$sample".IG.bam

# sort unmapped reads
samtools sort -@ $threads -n -o "$sample".unmapped.sorted.bam "$sample".unmapped.bam

# get unmapped fastq
bedtools bamtofastq -i "$sample".unmapped.sorted.bam -fq "$sample".unmapped.R1.fastq -fq2 "$sample".unmapped.R2.fastq

# get IG fastq
bedtools bamtofastq -i "$sample".IG.sorted.bam -fq "$sample".IG.R1.fastq -fq2 "$sample".IG.R2.fastq

# get combined fastq
cat "$sample".unmapped.R1.fastq "$sample".IG.R1.fastq > "$sample".combined.IG.R1.fastq
cat "$sample".unmapped.R2.fastq "$sample".IG.R2.fastq > "$sample".combined.IG.R2.fastq

# mixcr align
mixcr align -f -r "$report" -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$sample".combined.IG.R1.fastq "$sample".combined.IG.R2.fastq "$sample".IG.vdjca

# mixcr assemblePartial
mixcr assemblePartial -f "$sample".IG.vdjca "$sample".IG.rescued.vdjca
mixcr assemblePartial -f "$sample".IG.rescued.vdjca "$sample".IG.rescued2.vdjca

# mixcr assemble
mixcr assemble -f "$sample".IG.rescued2.vdjca "$sample".IG.clns

# mixcr exportClones
mixcr exportClones -f -o -t "$sample".IG.clns "$output"
