#!/bin/bash

pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17


WORKDIR=`pwd`
REF=/gscmnt/gc2802/halllab/sv_aggregate/refs/all_sequences.fa
EXCLUDE=/gscmnt/gc2802/halllab/sv_aggregate/exclusion/exclude.cnvnator_100bp.112015.bed
SVTOOLS=/gscmnt/gc2802/halllab/ccdg_resources/bin/svtools-0.3.1
SAMPLE_MAP=sample.map
HALLLAB_BIN=/gscmnt/gc2719/halllab/bin

mkdir -p lumpy_non_ref/logs


# GTEX-WHWD-0002-SM-5SOEC has a weird library and failed. omitted from  /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/merged_2015-06-26/merge_sample_list.txt

# the following samples have RNA-seq contamination, so omitting from /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/merged_2015-06-26/merge_sample_list.txt
# GTEX-QCQG-0003
# GTEX-OOBJ-0003
# GTEX-PX3G-0004
# GTEX-S7SF-0003
# GTEX-XUYS-0001
# GTEX-SNMC-0004
# GTEX-XUZC-0002
# GTEX-WRHU-0003
# GTEX-XQ8I-0001

# the following samples have weird read-depth profiles, omitting them (MAD of CN (every 100 bp) is greater than 0.5)
# GTEX-11EI6-0004
# GTEX-12584-0004
# GTEX-11DZ1-0002
# GTEX-11DXZ-0003
# GTEX-1211K-0003
# GTEX-11DYG-0004
# GTEX-12BJ1-0001
# GTEX-11EQ9-0003
# GTEX-11GSP-0004
# GTEX-14PJ2-0002
# GTEX-11DXY-0004
# GTEX-11GS4-0004
# GTEX-14PJM-0003
# GTEX-11EQ8-0001
# GTEX-WHSE-1726
# GTEX-11DXX-0004
# GTEX-11P81-0002
# GTEX-12WSG-0001
# GTEX-1192W-0001
# GTEX-NPJ8-0004



# remove homoygous reference SVs from VCF
while read TRUNC SAMPLE BAM LUMPYDIR 
do
  mkdir -p lumpy_non_ref/$TRUNC
  bsub -g /ajscott/sv_aggregate -e lumpy_non_ref/logs/$SAMPLE.log -o lumpy_non_ref/logs/$SAMPLE.log -q long  "cat $LUMPYDIR/gt/$SAMPLE.vcf \
  | vawk --header '{if(S\$*\$GT!=\"0/0\" && S\$*\$GT!=\"./.\") print \$0}' \
  > lumpy_non_ref/$TRUNC/$TRUNC.sv.non_ref.vcf "
done < notes/sample.map.485

while read TRUNC SAMPLE BAM LUMPYDIR 
do
  mkdir -p lumpy_non_ref/$TRUNC
  bsub -g /ajscott/sv_aggregate -e lumpy_non_ref/logs/$SAMPLE.log -o lumpy_non_ref/logs/$SAMPLE.log -q long  "cat /gscmnt/gc2802/halllab/sv_aggregate/mega_round3/cohorts/GTEX/$SAMPLE/$SAMPLE.sv.vcf \
  | vawk --header '{if(S\$*\$GT!=\"0/0\" && S\$*\$GT!=\"./.\") print \$0}' \
  > lumpy_non_ref/$TRUNC/$TRUNC.sv.non_ref.vcf "
done < notes/sample.map.146

while read TRUNC SAMPLE BAM LUMPYDIR
do
ls -cltr lumpy_non_ref/$TRUNC/$TRUNC.sv.non_ref.vcf
done < $SAMPLE_MAP > logs/lumpy_non_ref.txt



# concatenate and sort the variants
echo -n "$SVTOOLS lsort " > sort_cmd.sh
while read TRUNC SAMPLE BAM LUMPYDIR
do
 echo -ne " \\\\\n\tlumpy_non_ref/$TRUNC/$TRUNC.sv.non_ref.vcf"
done >> sort_cmd.sh < $SAMPLE_MAP
echo "" >> sort_cmd.sh
chmod 755 sort_cmd.sh

HNUM=`bsub -e logs/sort.%J.log -o logs/sort.%J.log -M 32000000 -R 'select[mem>32000] rusage[mem=32000]' -q long -g /ajscott/sv_aggregate "bash sort_cmd.sh | bgzip -c > sv.sorted.vcf.gz " | perl -ne 'print "$1\n" if /<([\d]+)>/' `

bsub -M 48000000 -q long  -R 'select[mem>48000] rusage[mem=48000]' -e logs/merge.%J.log -o logs/merge.%J.log  -w "done($HNUM)" \
"zcat sv.sorted.vcf.gz \
 | $SVTOOLS lmerge -i /dev/stdin -f 20 \
 | bgzip -c > sv.merged.vcf.gz "



mkdir -p gt/logs

while read TRUNC SAMPLE BAM LUMPYDIR
do
	echo $TRUNC
	base=`basename $BAM .bam`
	ROOT=$LUMPYDIR/$SAMPLE/temp/cnvnator-temp/$base.bam.hist.root
	SM=`sambamba view -H $BAM | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^SM:") SM=$i; gsub("^SM:","",SM); } print SM }'`
	bsub -M 30000000 -q long -R 'select[mem>30000] rusage[mem=30000]' -J $SAMPLE.gt -g /ajscott/sv_aggregate -o gt/logs/$SAMPLE.gt.%J.log -e gt/logs/$SAMPLE.gt.%J.log \
	"zcat sv.merged.vcf.gz \
	| vawk --header '{  \$6==\".\"; print }' \
	| $SVTOOLS genotype \
		-B $BAM \
		-l $BAM.json \
	| sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
	> gt/$SAMPLE.vcf"
	sleep 1
done < $SAMPLE_MAP
exit 0

mkdir -p cn/logs
mkdir -p status

# generate coordinate bed file required for cn annotations
zcat sv.merged.vcf.gz | /gscmnt/gc2802/halllab/ccdg_resources/bin/create_coordinates-0.3.1 -o coordinates.txt

#ROOT libraries path
source /gsc/pkg/root/root/bin/thisroot.sh

while read TRUNC SAMPLE BAM LUMPYDIR
do
	echo $TRUNC
	base=`basename $BAM .bam`
	ROOT=$LUMPYDIR/$SAMPLE/temp/cnvnator-temp/$base.bam.hist.root
	SM=`sambamba view -H $BAM | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^SM:") SM=$i; gsub("^SM:","",SM); } print SM }'`
	bsub -M 30000000 -q long -R 'select[mem>30000] rusage[mem=30000]' -J $SAMPLE.gt -g /ajscott/sv_aggregate -o cn/logs/$SAMPLE.gt.%J.log -e cn/logs/$SAMPLE.gt.%J.log \
	"$SVTOOLS copynumber \
	    --cnvnator cnvnator-multi \
	    -s $SM \
	    -w 100 \
	    -r $ROOT \
	    -c coordinates.txt \
	    -i gt/$SAMPLE.vcf \
	 > cn/$SAMPLE.vcf"
	sleep 1
done < $SAMPLE_MAP
exit 0

# check to make sure all of the files are the same length
while read TRUNC SAMPLE BAM LUMPYDIR
do
 bsub -q long \
 " cat cn/$SAMPLE.vcf | wc -l > status/$SAMPLE.txt  "
done < $SAMPLE_MAP
exit 0

# paste individual vcfs together
ls cn/*.vcf > paste.txt

bsub -q long -M 48000000 -R 'select[mem>48000] rusage[mem=48000]'  -J paste.cn -o logs/paste.%J.log -e logs/paste.%J.log \
   "$SVTOOLS vcfpaste \
       -m sv.merged.vcf.gz \
       -f paste.txt \
       -q \
       | bgzip -c \
       > sv.merged.pasted.new.vcf.gz "
exit 0


bsub -q long -e logs/prune.%J.log -o logs/prune.%J.log -M 64000000 -R 'select[mem>64000] rusage[mem=64000]'  \
	"zcat sv.merged.pasted.vcf.gz \
	| $SVTOOLS afreq \
	| $SVTOOLS vcftobedpe \
	| $SVTOOLS bedpesort \
	| $SVTOOLS prune -s -d 100 -e "AF" \
	| $SVTOOLS bedpetovcf \
	| bgzip -c \
	> sv.pruned.vcf.gz"


# filter batch effects
bsub -e logs/batch.%J.log -o logs/batch.%J.log -M 1000000 -R 'select[mem>1000] rusage[mem=1000]' -q long \
"cat /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/qc/flagstat_qc.plusOld.txt \
  | zjoin -a stdin -b sample.map -1 1 -2 2 -w a \
  | awk '{OFS=\"\t\"; split(\$1,a,\"-\"); print \"GTEX-\"a[2]\"-\"a[3],\$2,\$3,\$4,\$5,\$6}' \
  | awk '{ s[NR]=\$1; m[NR]=\$5 } END { printf \"ID\"; for (i=1;i<=length(s);++i) printf \"\t\"s[i] ; printf \"\n\"; printf \"PLATFORM\"; for (i=1;i<=length(m);++i) printf \"\t\"m[i] ; printf \"\n\" }' \
  > gtex_wgs_covar.long.txt"

bsub -e logs/filterBatch.%J.log -o logs/filterBatch.%J.log -M 1000000 -R 'select[mem>1000] rusage[mem=1000]' -q long \
"zcat sv.pruned.vcf.gz \
    | vawk --header 'I\$AF>0' \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/vcf_covar.py -c gtex_wgs_covar.long.txt \
    | bgzip -c \
    > sv.pruned.covar.vcf.gz"

# flag deletions that are supported only by paired-end reads and have size less than 423
# 423 is the minimum deletion size at which every library has distro overlap of less than 5% between concordant and discordant
bsub -e logs/flagDels.%J.log -o logs/flagDels.%J.log -M 4000000 -R 'select[mem>4000] rusage[mem=4000]' -q long \
"zcat sv.pruned.covar.vcf.gz \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i \"LOW_RES_PE\" -c FILTER -d \"Deletions supported only by paired-end evidence that are smaller than paired-end resolution (423 bp)\" \
    | vawk --header '{ if (I\$SVTYPE==\"DEL\" && I\$SR==0 && I\$SVLEN>-423) \$7=\"LOW_RES_PE\"; print }' \
    | $SVTOOLS vcfsort \
    | bgzip -c \
    > sv.pruned.covar.pe_filter.vcf.gz; \
tabix -p vcf -f sv.pruned.covar.pe_filter.vcf.gz"


############################################################################################################################################################
# Check for overlap w/ 1KGP
# 25% recip overlap and type matching validates any SV - go stricter for lumpy calls (check 50 & 90% overlap)

# convert vcf to bed format
zcat sv.pruned.covar.pe_filter.vcf.gz | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/vcfToBed.py | bgzip -c > 1KGP_comparison/sv.pruned.covar.pe_filter.bed.gz

zcat 1KGP_comparison/sv.pruned.covar.pe_filter.bed.gz \
  | vawk -c 9 '{print $1,$2,$3,$4,I$SVTYPE}' \
  | bedtools intersect -f 0.9 -r -wo -a stdin -b /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/1kg_sv_truthset_2015-12-15/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.bed.gz \
  | vawk -c 14 '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,I$SVTYPE}' \
  | awk '{if($5==$11){print $0} else if($5=="DEL" && $11~/DEL/){print $0} else if($5=="DEL" && $10 ~ /CN0/){print $0} else if($5=="DUP" && $11=="CNV"){print $0} else if($5=="BND"){print $0}}' \
  > 1KGP_comparison/gtex_merged.lumpy.1KGP_0.9.txt
wc -l 1KGP_comparison/gtex_merged.lumpy.1KGP_0.9.txt
# 12273

zcat 1KGP_comparison/sv.pruned.covar.pe_filter.bed.gz \
  | vawk -c 9 '{print $1,$2,$3,$4,I$SVTYPE}' \
  | bedtools intersect -f 0.5 -r -wo -a stdin -b /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/1kg_sv_truthset_2015-12-15/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.bed.gz \
  | vawk -c 14 '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,I$SVTYPE}' \
  | awk '{if($5==$11){print $0} else if($5=="DEL" && $11~/DEL/){print $0} else if($5=="DEL" && $10 ~ /CN0/){print $0} else if($5=="DUP" && $11=="CNV"){print $0} else if($5=="BND"){print $0}}' \
  > 1KGP_comparison/gtex_merged.lumpy.1KGP_0.5.txt
wc -l 1KGP_comparison/gtex_merged.lumpy.1KGP_0.5.txt
# 15862

# 50% overlap looks like the way to go after looking at overlaps in IGV - flag these variants in our dataset
zcat sv.pruned.covar.pe_filter.vcf.gz \
  | grep -v "#" \
  | zjoin -a stdin -b 1KGP_comparison/gtex_merged.lumpy.1KGP_0.5.txt -1 3 -2 4 -w a \
  | vawk '{$7="MATCH_1KGP"; print}' \
  | bgzip -c \
  > 1KGP_comparison/sv.pruned.covar.pe_filter.1KGP_0.5.vcf.gz

zcat sv.pruned.covar.pe_filter.vcf.gz \
  | grep -v "#" \
  | zjoin -a stdin -b 1KGP_comparison/gtex_merged.lumpy.1KGP_0.5.txt -1 3 -2 4 -v \
  | bgzip -c \
  > 1KGP_comparison/sv.pruned.covar.pe_filter.1KGP_0.5.unmatched.vcf.gz

# get vcf header
zcat sv.pruned.covar.pe_filter.vcf.gz | head -100 | grep "^#" | bgzip -c > 1KGP_comparison/vcf_header.txt.gz

bsub -e logs/1KGP.%J.log -o logs/1KGP.%J.log -M 4000000 -R 'select[mem>4000] rusage[mem=4000]' -q long \
"zcat 1KGP_comparison/vcf_header.txt.gz 1KGP_comparison/sv.pruned.covar.pe_filter.1KGP_0.5.vcf.gz 1KGP_comparison/sv.pruned.covar.pe_filter.1KGP_0.5.unmatched.vcf.gz \
  | $SVTOOLS vcfsort \
  | vawk --header '{if(\$7==\"MATCH_1KGP\"){\$8=\$8\";MATCH=1\"; print} else if(\$7==\"LOW_RES_PE\"){\$8=\$8\";LOWRES=1\"; print} else{print}}' \
  | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i \"MATCH_1KGP\" -c FILTER -d \"SVs validated by 1KGP callset\" \
  | uniq \
  | bgzip -c \
  > sv.pruned.covar.pe_filter.1KGP_0.5.vcf.gz"

############################################################################################################################################################


# reclassify (using large-sample method, here) and mark low_qual vars
bsub -e logs/reclass.%J.log -o logs/reclass.%J.log -M 32000000 -R 'select[mem>32000] rusage[mem=32000]' \
"zcat sv.pruned.covar.pe_filter.1KGP_0.5.vcf.gz  \
  | $SVTOOLS vcfsort \
  | $SVTOOLS classify \
   -m large_sample \
   -a /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz \
   -g /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/gtex.sex.txt \
  | vawk --header '{
     split(I\$STRANDS,x,\",\");
     split(x[1],y,\":\");
     split(x[2],z,\":\");
     if (I\$SVTYPE==\"DEL\" || I\$SVTYPE==\"DUP\" || I\$SVTYPE==\"MEI\"){
       \$7=\"PASS\"; print \$0;
     }  else if ( I\$SVTYPE==\"INV\" && \$8>=100 && (I\$SR/I\$SU)>=0.1 && (I\$PE/I\$SU)>=0.1 && (y[2]/I\$SU)>0.1 && (z[2]/I\$SU)>0.1){
       \$7=\"PASS\"; print \$0;
     } else if ( I\$SVTYPE==\"BND\" && \$8>=100 && (I\$SR/I\$SU)>=0.25 && (I\$PE/I\$SU)>=0.25){
       \$7=\"PASS\"; print \$0;
     } else {
       \$7=\"LOW\"; print \$0;}}' \
  | bgzip -c > sv.pruned.reclassed.vcf.gz"


# Convert SV type tags back to standard format, add appropriate flags
bsub -e logs/fixReclass.%J.log -o logs/fixReclass.%J.log -M 1000000 -R 'select[mem>1000] rusage[mem=1000]' -q long \
"zcat sv.pruned.reclassed.vcf.gz \
  | vawk --header '{if(I\$MATCH==1){\$7=\"MATCH_1KGP\"; \$8=substr(\$8,1,length(\$8)-8); print} else if(I\$LOWRES==1){\$7=\"LOW_RES_PE\";\$8=substr(\$8,1,length(\$8)-9); print} else{print}}' \
  | bgzip -c \
  > sv.pruned.pe_filter.reclassed.fixed.vcf.gz"

bsub -e logs/dosage.%J.log -o logs/dosage.%J.log -M 4000000 -R 'select[mem>4000] rusage[mem=4000]' -q long \
"zcat sv.pruned.pe_filter.reclassed.fixed.vcf.gz \
  | stripprob \
  | vawk --header '{ \$3=\"LUMPY_\"I\$SVTYPE\"_\"\$3; if (I\$SVTYPE==\"BND\") { gsub(\";MATEID=\",\";MATEID=LUMPY_\"I\$SVTYPE\"_\",\$8); gsub(\";EVENT=\",\";EVENT=LUMPY_\"I\$SVTYPE\"_\",\$8) } else gsub(\";EVENT[^;]\",\"\",\$8); print }' \
  | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i DS -c FORMAT -t Float -n A -d \"Dosage of alt allele\" \
  | awk '{ if (\$0~\"^#\") { print; next } \$9=\$9\":DS\"; split(\$9,fmt,\":\"); for (i=1;i<=length(fmt);++i) { if (fmt[i]==\"AB\") fmt_idx=i } for (i=10;i<=NF;++i) { split(\$i,gt,\":\"); \$i=\$i\":\"gt[fmt_idx]*2 } print }' OFS=\"\t\" \
  | $SVTOOLS vcfsort \
  | sed 's/GTEX-\(....\)-..../GTEX-\1/g' \
  | bgzip -c \
  > sv.pruned.pe_filter.reclassed.dosage.vcf.gz; \
tabix -p vcf -f sv.pruned.pe_filter.reclassed.dosage.vcf.gz"

# apply filters
# filter variants with > 0.1 explained variance by platform and BND or INV variants with qual < 20
# Remove SNAME
bsub -e logs/filter.%J.log -o logs/filter.%J.log -M 2000000 -R 'select[mem>2000] rusage[mem=2000]' -q long \
"zcat sv.pruned.pe_filter.reclassed.dosage.vcf.gz \
  | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i \"PLATFORM\" -c FILTER -d \"Variant genotypes correlate with sequencing platform with R^2 greater than 0.1\" \
  | vawk --header '{ if (I\$EXVAR>0.1) \$7=\"PLATFORM\"; print }' \
  | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i \"MSQ_20\" -c FILTER -d \"Variant without read-depth support with MSQ less than 20\" \
  | vawk --header '{ if ((I\$SVTYPE==\"BND\" || I\$SVTYPE==\"INV\") && I\$MSQ<20) \$7=\"MSQ_20\"; print }' \
  | vawk --header '{ if (I\$SVTYPE==\"BND\") gsub(\";SVLEN=[^;]*;\", \";\", \$8); print }' \
  | $SVTOOLS vcfsort \
  | uniq \
  | bgzip -c \
  > gtex_merged.filt_annot.vcf.gz; \
tabix -p vcf -f gtex_merged.filt_annot.vcf.gz"

# remove flagged variants
bsub -e logs/low.%J.log -o logs/low.%J.log -M 2000000 -R 'select[mem>2000] rusage[mem=2000]' -q long  \
"zcat gtex_merged.filt_annot.vcf.gz \
  | vawk --header '\$7==\"PASS\" || \$7==\"LOW\" || \$7==\"MATCH_1KGP\"' \
  | bgzip -c \
  > gtex_merged.filt_annot.low_pass.vcf.gz; \
tabix -p vcf -f gtex_merged.filt_annot.low_pass.vcf.gz"

# fix sample names in header
zcat gtex_merged.filt_annot.low_pass.vcf.gz \
	| awk '{OFS="\t"; if($1~/^#CHROM/){line=$1; for(x=2;x<10;x++){line=line"\t"$x}; for(x=10;x<=NF;x++){split($x,a,"-"); line=line"\tGTEX-"a[2]}; print line} else{print $0}}' \
	| bgzip -c \
	> gtex_merged.filt_annot.low_pass.fixed.vcf.gz
	
mv gtex_merged.filt_annot.low_pass.fixed.vcf.gz gtex_merged.filt_annot.low_pass.vcf.gz
tabix -p vcf -f gtex_merged.filt_annot.low_pass.vcf.gz






