# get the areas of interest from bam file into multiple bams based on conditions
# get mates somehow and make into paired end fastq

sample=$1
threads=$2

report=$3
output=$4

# get index
#samtools index -@ $threads "$sample"

# get unmapped reads
samtools view -b -h -f 4 -@ $threads "$sample" > "$sample".unmapped.bam

# get TRA reads
samtools view -b -h -@ $threads "$sample" "chr14:21621904-22552132" > "$sample".TRA.bam

# get TRB reads
samtools view -b -h -@ $threads "$sample" "chr7:142299011-142813287" > "$sample".TRB.bam

# get TRG reads
samtools view -b -h -@ $threads "$sample" "chr7:38240024-38368055" > "$sample".TRG.bam

# get TRD reads
samtools view -b -h -@ $threads "$sample" "chr14:22422546-22466577" > "$sample".TRD.bam

# merge TR reads
samtools merge -f -@ $threads "$sample".TR.bam "$sample".TRA.bam "$sample".TRB.bam "$sample".TRG.bam "$sample".TRD.bam

# sort TR reads
samtools sort -@ $threads -n "$sample".TR.bam "$sample".TR.sorted

# sort unmapped reads
samtools sort -n "$sample".unmapped.bam "$sample".unmapped.sorted

# get unmapped fastq
bedtools bamtofastq -i "$sample".unmapped.sorted.bam -fq "$sample".unmapped.R1.fastq -fq2 "$sample".unmapped.R2.fastq

# get TR fastq
bedtools bamtofastq -i "$sample".TR.sorted.bam -fq "$sample".TR.R1.fastq -fq2 "$sample".TR.R2.fastq

# get combined fastq
cat "$sample".unmapped.R1.fastq "$sample".TR.R1.fastq > "$sample".combined.T.R1.fastq
cat "$sample".unmapped.R2.fastq "$sample".TR.R2.fastq > "$sample".combined.T.R2.fastq

# mixcr align
mixcr align -f -r "$report" -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$sample".combined.T.R1.fastq "$sample".combined.T.R2.fastq "$sample".T.vdjca

# mixcr assemblePartial
mixcr assemblePartial -f "$sample".T.vdjca "$sample".T.rescued.vdjca
mixcr assemblePartial -f "$sample".T.rescued.vdjca "$sample".T.rescued2.vdjca

# mixcr assemble
mixcr assemble -f "$sample".T.rescued2.vdjca "$sample".T.clns

# mixcr exportClones
mixcr exportClones -f -o -t "$sample".T.clns "$output"

# parse output
cat "$sample".T.clones.tsv | grep -v 'IGH' | grep -v 'IGK' | grep -v 'IGL' > "$sample".T.clones.tsv.filter1
head -1 "$sample".T.clones.tsv.filter1 > "$sample".T.clones.tsv.head1
grep -m 1 "TRA" "$sample".T.clones.tsv.filter1 > "$sample".T.clones.tsv.TRA
grep -m 1 "TRB" "$sample".T.clones.tsv.filter1 > "$sample".T.clones.tsv.TRB
grep -m 1 "TRG" "$sample".T.clones.tsv.filter1 > "$sample".T.clones.tsv.TRG
grep -m 1 "TRD" "$sample".T.clones.tsv.filter1 > "$sample".T.clones.tsv.TRD
cat "$sample".T.clones.tsv.head1 "$sample".T.clones.tsv.TRA "$sample".T.clones.tsv.TRB "$sample".T.clones.tsv.TRG "$sample".T.clones.tsv.TRD > "$sample".T.clones.tsv.filter2
cut -f 2,3,4,6,7,8,9 "$sample".T.clones.tsv.filter2 > "$sample".T.clones.tsv