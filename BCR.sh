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
samtools view -b -h -f 4 -@ $threads "$bam_dir"/"$sample".bam > "$result_dir"/"$sample".unmapped.bam

# get IGH reads
samtools view -b -h -@ $threads "$bam_dir"/"$sample".bam "chr14:105586437-106879844" > "$result_dir"/"$sample".IGH.bam

# get IGK reads
samtools view -b -h -@ $threads "$bam_dir"/"$sample".bam "chr2:88857361-90235368" > "$result_dir"/"$sample".IGK.bam

# get IGL reads
samtools view -b -h -@ $threads "$bam_dir"/"$sample".bam "chr22:22026076-22922913" > "$result_dir"/"$sample".IGL.bam

# merge IG reads
samtools merge -f -@ $threads "$result_dir"/"$sample".IG.bam "$result_dir"/"$sample".IGH.bam "$result_dir"/"$sample".IGK.bam "$result_dir"/"$sample".IGL.bam

# sort IG reads
samtools sort -@ $threads -n -o "$result_dir"/"$sample".IG.sorted.bam "$result_dir"/"$sample".IG.bam

# sort unmapped reads
samtools sort -@ $threads -n -o "$result_dir"/"$sample".unmapped.sorted.bam "$result_dir"/"$sample".unmapped.bam

# get unmapped fastq
bedtools bamtofastq -i "$result_dir"/"$sample".unmapped.sorted.bam -fq "$result_dir"/"$sample".unmapped.R1.fastq -fq2 "$result_dir"/"$sample".unmapped.R2.fastq

# get IG fastq
bedtools bamtofastq -i "$result_dir"/"$sample".IG.sorted.bam -fq "$result_dir"/"$sample".IG.R1.fastq -fq2 "$result_dir"/"$sample".IG.R2.fastq

# get combined fastq
cat "$result_dir"/"$sample".unmapped.R1.fastq "$result_dir"/"$sample".IG.R1.fastq > "$result_dir"/"$sample".combined.IG.R1.fastq
cat "$result_dir"/"$sample".unmapped.R2.fastq "$result_dir"/"$sample".IG.R2.fastq > "$result_dir"/"$sample".combined.IG.R2.fastq

# mixcr align
mixcr align -f -r "$result_dir"/"$sample".report -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$result_dir"/"$sample".combined.IG.R1.fastq "$result_dir"/"$sample".combined.IG.R2.fastq "$result_dir"/"$sample".IG.vdjca

# mixcr assemblePartial
mixcr assemblePartial -f "$result_dir"/"$sample".IG.vdjca "$result_dir"/"$sample".IG.rescued.vdjca
mixcr assemblePartial -f "$result_dir"/"$sample".IG.rescued.vdjca "$result_dir"/"$sample".IG.rescued2.vdjca

# mixcr assemble
mixcr assemble -f "$result_dir"/"$sample".IG.rescued2.vdjca "$result_dir"/"$sample".IG.clns

# mixcr exportClones
mixcr exportClones -f -o -t "$result_dir"/"$sample".IG.clns "$result_dir"/"$sample".IG.clones.tsv

