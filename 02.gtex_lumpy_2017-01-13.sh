# 2017-01-13

# ========================================
# Run SVTYPER on a batch of GTEx samples
# ========================================

pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/lumpy_2016-06-29

# ----------------------------------------
# 1. Create sample map
# ----------------------------------------

cat notes/gtex_504.txt \
	| cut -f 1 \
	| awk '{OFS="\t"; print $1,"/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/"$1"/"$1".bam","/gscmnt/gc2802/halllab/gtex_2016-06-01/lumpy_2016-06-29/"$1}' \
	> notes/sample.map

# make output directories
mkdir -p /gscmnt/gc2802/halllab/gtex_2016-06-01/lumpy_2016-06-29/gt
mkdir -p /gscmnt/gc2802/halllab/gtex_2016-06-01/lumpy_2016-06-29/gt/logs


# ----------------------------------------
# 3. Run svtyper
# ----------------------------------------

#using svtyper v0.1.1

while read SAMPLE BAM LUMPYDIR
do
 echo $SAMPLE
 base=`basename $BAM .bam`
 bsub -M 30000000 -q long -R 'select[mem>30000] rusage[mem=30000]' -J $SAMPLE.gt -g /ajscott/svtyper -o /gscmnt/gc2802/halllab/gtex_2016-06-01/lumpy_2016-06-29/gt/logs/$SAMPLE.gt.%J.log -e /gscmnt/gc2802/halllab/gtex_2016-06-01/lumpy_2016-06-29/gt/logs/$SAMPLE.gt.%J.log \
 "zcat $LUMPYDIR/$SAMPLE.sv.vcf.gz \
  | /gscmnt/gc2802/halllab/sv_aggregate/svtyper_repo/svtyper/svtyper \
    -B $BAM \
    > /gscmnt/gc2802/halllab/gtex_2016-06-01/lumpy_2016-06-29/gt/$SAMPLE.vcf "
done < notes/sample.map


