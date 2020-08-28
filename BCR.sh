# get the areas of interest from bam file into multiple bams based on conditions
# get mates somehow and make into paired end fastq

sample=$1
threads=$2

report=$3
output=$4

echo "ARGUMENTS"
echo "sample:" "$sample"
echo "threads:" "$threads"
echo "report:" "$report"
echo "output:" "$output"
echo "========="
echo

# get index
# echo "indexing started"
# samtools index -@ $threads "$sample"
# sleep 5
# echo "indexing done"

# get unmapped reads
echo "unmapped started"
time samtools view -b -h -f 4 -@ $threads "$sample" > "$sample".unmapped.bam
sleep 5
echo "unmapped done"

# get IGH reads
echo "IGH started"
time samtools view -b -h -@ $threads "$sample" "chr14:105586437-106879844" > "$sample".IGH.bam
sleep 5
echo "IGH done"

# get IGK reads
echo " IGK started"
time samtools view -b -h -@ $threads "$sample" "chr2:88857361-90235368" > "$sample".IGK.bam
sleep 5
echo "IGK done"

# get IGL reads
echo "IGL started"
time samtools view -b -h -@ $threads "$sample" "chr22:22026076-22922913" > "$sample".IGL.bam
sleep 5
echo "IGL done"

# merge IG reads
echo "merge started"
time samtools merge -f -@ $threads "$sample".IG.bam "$sample".IGH.bam "$sample".IGK.bam "$sample".IGL.bam
sleep 5
echo "merge done"

# sort IG reads
echo "sort started"
time samtools sort -@ $threads -n "$sample".IG.bam "$sample".IG.sorted
sleep 5
echo "sort done"

echo "more processing started"

# sort unmapped reads
samtools sort -@ $threads -n "$sample".unmapped.bam "$sample".unmapped.sorted

# get unmapped fastq
bedtools bamtofastq -i "$sample".unmapped.sorted.bam -fq "$sample".unmapped.R1.fastq -fq2 "$sample".unmapped.R2.fastq 2> tmp

# get IG fastq
bedtools bamtofastq -i "$sample".IG.sorted.bam -fq "$sample".IG.R1.fastq -fq2 "$sample".IG.R2.fastq 2> tmp

# get combined fastq
cat "$sample".unmapped.R1.fastq "$sample".IG.R1.fastq > "$sample".combined.IG.R1.fastq
cat "$sample".unmapped.R2.fastq "$sample".IG.R2.fastq > "$sample".combined.IG.R2.fastq

sleep 5
echo "more processing done"

echo "mixcr started"

# mixcr align
mixcr align -f -r "$report" -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$sample".combined.IG.R1.fastq "$sample".combined.IG.R2.fastq "$sample".IG.vdjca

# mixcr assemblePartial
mixcr assemblePartial -f "$sample".IG.vdjca "$sample".IG.rescued.vdjca
mixcr assemblePartial -f "$sample".IG.rescued.vdjca "$sample".IG.rescued2.vdjca

# mixcr assemble
mixcr assemble -f "$sample".IG.rescued2.vdjca "$sample".IG.clns

# mixcr exportClones
mixcr exportClones -count -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -aaFeature CDR3 "$sample".IG.clns "$sample".IG.clones.tsv
# mixcr exportClones -f -o -t "$sample".IG.clns "$sample".IG.clones.tsv

sleep 5
echo "mixcr done"

echo "parse output started"

# parse output
cat "$sample".IG.clones.tsv | grep -v 'TRA' | grep -v 'TRB' | grep -v 'TRG' | grep -v 'TRD' > "$sample".IG.clones.tsv.filter1
head -1 "$sample".IG.clones.tsv.filter1 > "$sample".IG.clones.tsv.head1
grep -m 1 "IGH" "$sample".IG.clones.tsv.filter1 > "$sample".IG.clones.tsv.IGH
grep -m 1 "IGK" "$sample".IG.clones.tsv.filter1 > "$sample".IG.clones.tsv.IGK
grep -m 1 "IGL" "$sample".IG.clones.tsv.filter1 > "$sample".IG.clones.tsv.IGL
cat "$sample".IG.clones.tsv.head1 "$sample".IG.clones.tsv.IGH "$sample".IG.clones.tsv.IGK "$sample".IG.clones.tsv.IGL > "$output"
# "$sample".IG.clones.tsv.filter2
# cut -f 2,3,4,6,7,8,9 "$sample".IG.clones.tsv.filter2 > "$sample".IG.clones.tsv.filter3
# sed 's/[(][^)]*[)]//g' "$sample".IG.clones.tsv.filter3 > "$output"

sleep 5
echo "parse output done"