# 2016-09-01

# generate a unified callset with LUMPY and Genome STRiP calls

# -----------------------------------
# strategy:
    # overlap
    #     Visualize in IGV
    #     keep in mind that GS variants are imprecise, perhaps 1kb confidence window around the breakpoints
    # prioritizing
    #     When GS and LUMPY overlap, choose the LUMPY variant because the coordinates are more precise
    #     if it's a BND, choose the GS variant because it will contain read depth information


pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/combine_lumpy_gs_2017-01-26/

GS=/gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gtex.genome_strip.filt_annot.merged.low_pass.vcf.gz
LUMPY=/gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/gtex_merged.filt_annot.low_pass.vcf.gz

# make the combined header
zcat $GS \
    | awk '{ if ($0~"^#") print; else exit }' \
    | cat - <(zcat $LUMPY | awk '{ if ($0~"^#") print; else exit }') \
    > combined.header.txt

# manually sort the fields in this file, and add GSMERGE and TOOL INFO fields

# make some BED files
zcat $GS | vawk '{ print $1,$2,I$END,$3,I$GSCNCATEGORY }' | bgzip -c > gs.bed.gz
zcat gs.bed.gz | awk '$5=="DEL"' | bgzip -c > gs.del.bed.gz
zcat gs.bed.gz | awk '$5=="DUP"' | bgzip -c > gs.dup.bed.gz
zcat gs.bed.gz | awk '$5=="MIXED"' | bgzip -c > gs.mixed.bed.gz
zcat $LUMPY | cut -f -8 | vcfToBedpe | bgzip -c > lumpy.bedpe.gz
zcat lumpy.bedpe.gz | bedpeToBed12 -i stdin -n lumpy | bgzip -c > lumpy.bed.gz

# merge GS and LUMPY variants with 0.5 reciprical overlap when variant type matches or MIXED.
FRAC=0.5
zcat lumpy.bedpe.gz \
    | vawk -c 13 '{ if ($1!=$4 || $9==$10 || $11=="BND") next; split(I$CIPOS,CIPOS,","); split(I$CIEND,CIEND,","); CHROM=$1; CSTART=$2-CIPOS[1]; CEND=$6-CIEND[2]; print $1,CSTART,CEND,$7,$8,$9,$10,$11,$12,$13 }' \
    | bedtools intersect -wo -r -f $FRAC -a stdin -b gs.bed.gz \
    | awk '($6=="+" && $7=="-" && $15=="DEL") || ($6=="-" && $7=="+" && $15=="DUP") || $15=="MIXED"' \
    | awk '{ print $0,$16/($13-$12) }' OFS="\t" \
    | sort -k14,14 \
    | groupBy -g 14 -c 17 -o max -full \
    | bgzip -c \
    > recip_f_$FRAC.bed.gz

# make a tab delim file of the variant pairs
zcat recip_f_$FRAC.bed.gz \
    | awk '{ print $4,$14 }' OFS="\t" \
    > recip_f_$FRAC.txt

# identify those which should be merged from the 0.5 reciprocal overlap calls
cat recip_f_$FRAC.txt | sort -k1,1 | groupBy -g 1 -c 2 -o collapse > gs.recip.f_0.5.collapse.txt

# 2.
# additionally, merge when 90% of the Genome strip call is contained within a LUMPY call, where the genotype concordance has R^2 of greater than 0.25
# when a GS call hits more than one LUMPY variant, choose the one that has the highest reciprocal overlap from LUMPY
FRAC=0.9
zcat lumpy.bedpe.gz \
    | vawk -c 13 '{ if ($1!=$4 || $9==$10 || $11=="BND") next; split(I$CIPOS,CIPOS,","); split(I$CIEND,CIEND,","); CHROM=$1; CSTART=$2-CIPOS[1]; CEND=$6-CIEND[2]; print $1,CSTART,CEND,$7,$8,$9,$10,$11,$12,$13 }' \
    | bedtools intersect -wo -f $FRAC -a gs.bed.gz -b stdin \
    | zjoin -v -a stdin -b recip_f_0.5.txt -w a -1 4 -2 2 \
    | awk '($11=="+" && $12=="-" && $5=="DEL") || ($11=="-" && $12=="+" && $5=="DUP") || $5=="MIXED"' \
    | awk '{ print $0,$16/($8-$7) }' OFS="\t" \
    | sort -k4,4 \
    | groupBy -g 4 -c 17 -o max -full \
    | bgzip -c \
    > gs.contained.f_$FRAC.bed.gz

# make a tab delim file of the variant pairs
zcat gs.contained.f_$FRAC.bed.gz \
    | awk '{ print $9,$4 }' OFS="\t" \
    > gs.contained.f_$FRAC.txt

# get the genotype correlations for the contained GS calls
/gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/src/var_gt_corr.py \
    -a $LUMPY \
    -b $GS \
    -v gs.contained.f_$FRAC.txt \
    -af AB -bf CN \
    > gs.contained.f_$FRAC.r_sq.txt

# identify those which should be merged from contained GS calls.
cat gs.contained.f_$FRAC.r_sq.txt | awk '$3**2>0.25' | sort -k1,1 | groupBy -g 1 -c 2 -o collapse > gs.contained.f_0.9.r_sq_0.25.collapse.txt

# ------------------------------------------
# Merge the files.

# first, quadruple check that the Genome STRiP and LUMPY vcf sample columns are in the same order
zcat $GS | grep -m 1 "^#CHROM" | tr '\t' '\n' | paste - <(zcat $LUMPY | grep -m 1 "^#CHROM" | tr '\t' '\n') | awk '$1!=$2'
zcat $GS | grep -m 1 "^#CHROM" | tr '\t' '\n' | paste - <(cat combined.header.txt | grep -m 1 "^#CHROM" | tr '\t' '\n') | awk '$1!=$2'


# Merge! and set low/high confidence
# For high confidence, GSCNQUAL >=4 (10%FDR)
NUM_COL=`zcat $LUMPY | grep -v "^#" | awk '{ print NF }' | head -n 1`
zcat $LUMPY $GS \
    | grep -v "#" \
    | zjoin -v -a stdin -b <(cat gs.recip.f_0.5.collapse.txt gs.contained.f_0.9.r_sq_0.25.collapse.txt | cut -f 2 | tr ',' '\n') -1 3 -2 1 \
    | zjoin -r -a stdin -b <(cat gs.recip.f_0.5.collapse.txt gs.contained.f_0.9.r_sq_0.25.collapse.txt | sort -k1,1 | groupBy -g 1 -c 2 -o collapse) -1 3 -2 1 -p "#" \
    | vawk -v NUM_COL=$NUM_COL '{ if ($(NUM_COL+2)!="NA") { $8=$8";GSMERGE="$(NUM_COL+2); } print }' \
    | cut -f -${NUM_COL} \
    | cat combined.header.txt - \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i LOW -c FILTER -d "Variants that do not meet strict quality thresholds (min GSCNQUAL cutoffs or GS merged, depth or annotation supported LUMPY variant, INV min QUAL 100 with 10% PE and SR contribution, BND min QUAL 100 with 25% PE and SR contribution" \
    | vawk --header '{split(I$STRANDS,x,","); split(x[1],y,":"); split(x[2],z,":"); \
        if($3~"^GS" && I$GSCNQUAL>=4){ $7="PASS"; print} \
        else if ($3~"^GS_" && I$GSMERGE!="") { $7="PASS"; print } \
        else if ($3~"^GS_"){ $7="LOW"; print } \
        else { print} }' \
    | /gscmnt/gc2802/halllab/ccdg_resources/bin/svtools-0.3.1 vcfsort \
    | bgzip -c \
    > ../gtex.lumpy.gs.low_pass.vcf.gz
tabix -f -p vcf ../gtex.lumpy.gs.low_pass.vcf.gz

# make full callset
zcat ../gtex.lumpy.gs.low_pass.vcf.gz \
    | cat - <(zcat /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/gtex_merged.filt_annot.vcf.gz /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gtex.genome_strip.filt_annot.merged.vcf.gz | grep -v "^#" | awk '$7!="." && $7!="PASS" && $7!="LOW" && $7!="MATCH_1KGP"') \
    | /gscmnt/gc2802/halllab/ccdg_resources/bin/svtools-0.3.1 vcfsort \
    | bgzip -c \
    > ../gtex.lumpy.gs.vcf.gz
tabix -p vcf -f ../gtex.lumpy.gs.vcf.gz

############################################################################################################################################################
# Check for GS overlap w/ 1KGP
# 25% recip overlap and type matching validates any SV (looked at 50% overlap in IGV & lose a lot of variants that seem to coincide w/ 1KGP)
cd ..
pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/

zcat gtex.lumpy.gs.vcf.gz | /gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/src/vcfToBed.py | awk '$4~/^GS/ || $1~"^#"' | awk '$8=="PASS" || $8=="LOW"' |  bgzip -c > 1KGP_comparison/gtex.gs.bed.gz

zcat 1KGP_comparison/gtex.gs.bed.gz \
  | vawk -c 9 '{print $1,$2,$3,$4,I$GSCNCATEGORY}' \
  | bedtools intersect -f 0.25 -r -wo -a stdin -b /gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/1kg_sv_truthset_2015-12-15/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.bed.gz \
  | vawk -c 14 '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,I$SVTYPE}' \
  | awk '$5==$11 || $11=="CNV" || ($5=="MIXED" && ($11 ~ "^DEL" || $11=="DUP")) || ($5=="DEL" && $11~"^DEL")' \
  > 1KGP_comparison/gtex.gs.1KGP_0.25.txt

# 25% overlap looks like the way to go after looking at overlaps in IGV - flag these variants in our dataset
zcat gtex.lumpy.gs.vcf.gz \
  | grep -v "#" \
  | zjoin -a stdin -b 1KGP_comparison/gtex.gs.1KGP_0.25.txt -1 3 -2 4 -w a \
  | vawk '{$7="MATCH_1KGP"; print}' \
  | /gscmnt/gc2802/halllab/ccdg_resources/bin/svtools-0.3.1 vcfsort \
  | uniq \
  | bgzip -c \
  > 1KGP_comparison/gtex.gs.1KGP_0.25.vcf.gz

zcat gtex.lumpy.gs.vcf.gz \
  | grep -v "#" \
  | zjoin -a stdin -b 1KGP_comparison/gtex.gs.1KGP_0.25.txt -1 3 -2 4 -v \
  | bgzip -c \
  > 1KGP_comparison/gtex.gs.1KGP_0.25.unmatched.vcf.gz

bomb -q long -m 4 -w 'done(777444)' \
"zcat combine_lumpy_gs_2017-01-26/combined.header.txt.gz 1KGP_comparison/gtex.gs.1KGP_0.25.vcf.gz 1KGP_comparison/gtex.gs.1KGP_0.25.unmatched.vcf.gz \
  | /gscmnt/gc2802/halllab/ccdg_resources/bin/svtools-0.3.1 vcfsort \
  | bgzip -c \
  > gtex.lumpy.gs.1KGP_0.5.vcf.gz; \
tabix -p vcf -f gtex.lumpy.gs.1KGP_0.5.vcf.gz"

############################################################################################################################################################

pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/

# Add allele frequency and NSAMP for genome strip calls
zcat gtex.lumpy.gs.1KGP_0.5.vcf.gz \
    | awk '$1 !~ /^GL/' \
    | vawk --header '$1 != "Y"' \
    | /gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/src/vcf_allele_freq_gs.py \
    | bgzip -c \
    >  gtex.lumpy.gs.1KGP_0.5.af.vcf.gz
tabix -p vcf gtex.lumpy.gs.1KGP_0.5.af.vcf.gz

# mark pseudoautosomal LUMPY DEL and DUPs as filtered

for VCF in gtex.lumpy.gs.1KGP_0.5.af.vcf.gz
do
    NCOL=`zcat $VCF | grep -v "^#" | head -n 1 | awk '{ print NF }'`
    tabix -h $VCF X Y | vawk '(I$SVTYPE=="DEL" || I$SVTYPE=="DUP") && I$NSAMP<=10' | vcfToBed | bedtools intersect -a stdin -b /gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/annotations/par.bed.gz \
        | cut -f 4 \
        | zjoin -r -a $VCF -b stdin -1 3 -2 1 -p "#" \
        | vawk --header '{ if ($NF!="NA") $7="PAR"; print }' \
        | cut -f -$NCOL \
        | grep -v "^#" \
        | cat combine_lumpy_gs_2017-01-26/combined.header.txt - \
        | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i "PAR" -c FILTER -d "Low frequency DEL or DUP in the pseudoautosomal region of X and Y" \
        | bgzip -c \
        > tmp.vcf.gz &&
    mv tmp.vcf.gz $VCF &&
    tabix -p vcf -f $VCF
done


# remove variants smaller than 50 bp
zcat gtex.lumpy.gs.1KGP_0.5.af.vcf.gz | awk '$3~/^LUMPY_BND/' | vcfToBedpe | awk '{OFS="\t"; if($1==$4){len=$6-$2; if(len<0){len=len*-1}; if(len<50){print $7}}}' > size_limit_50bp/bnd_under50.txt
zcat gtex.lumpy.gs.1KGP_0.5.af.vcf.gz | vawk --header '{if($3~/^LUMPY/ && $3!~/^LUMPY_BND/){len=I$SVLEN; if(len<0){len=len*-1}; if(len>=50){print $0}} else{print $0}}' | bgzip -c > size_limit_50bp/gtex.lumpy.gs.1KGP_0.5.size.withBND.vcf.gz
zcat size_limit_50bp/gtex.lumpy.gs.1KGP_0.5.size.withBND.vcf.gz | vawk --header '{if($3!~/^LUMPY_BND/){print $0}}' | bgzip -c > size_limit_50bp/gtex.lumpy.gs.1KGP_0.5.size.noBND.vcf.gz
zcat gtex.lumpy.gs.1KGP_0.5.af.vcf.gz | vawk '{if($3~/^LUMPY_BND/){print I$EVENT,$0}}' | zjoin -a stdin -b size_limit_50bp/bnd_under50.txt -1 1 -2 1 -v | cut -f 2- | bgzip -c > size_limit_50bp/bnd.size.vcf.gz
zcat size_limit_50bp/gtex.lumpy.gs.1KGP_0.5.size.noBND.vcf.gz size_limit_50bp/bnd.size.vcf.gz | /gscmnt/gc2802/halllab/ccdg_resources/bin/svtools-0.3.1 vcfsort | bgzip -c > gtex.lumpy.gs.1KGP_0.5.size.vcf.gz


# low conf
zcat gtex.lumpy.gs.1KGP_0.5.size.vcf.gz \
    | vawk --header '$7=="LOW" || $7=="PASS" || $7=="MATCH_1KGP"' \
    | bgzip -c \
    > gtex.lumpy.gs.1KGP_0.5.size.low_pass.vcf.gz
tabix -p vcf -f gtex.lumpy.gs.1KGP_0.5.size.low_pass.vcf.gz


# high conf
zcat gtex.lumpy.gs.1KGP_0.5.size.low_pass.vcf.gz \
    | vawk --header '$7=="PASS" || $7=="MATCH_1KGP"' \
    | bgzip -c \
    > gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz


# fix GS allele frequency to be based on GSNVARIANT instead of mode cn
zcat gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz \
    | vawk --header '{if($3~/^GS/){$8="AF="I$GSNVARIANT/631";MODE"$8; print $0} else{print $0}}' \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i MODEAF -c INFO -d "GS Allele frequency based on MODECN" \
    | bgzip -c \
    > gtex.lumpy.gs.1KGP_0.5.size.gscnqual.af.high_conf.vcf.gz
mv -f gtex.lumpy.gs.1KGP_0.5.size.gscnqual.af.high_conf.vcf.gz gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz
tabix -p vcf -f gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz


# remove redundant GS calls
# see gs_redundancy_2017-07-18.sh

# add MELT calls
# first check everything in same order
zcat gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz \
    | head -n 300 \
    | grep "#CHROM" \
    | cut -f 10- \
    | tr '\t' '\n' \
    | paste - <(zcat /gscmnt/gc2802/halllab/ajscott/gtex_v7/melt_mei/gtex.631.melt_mei.assess_5.platform.vcf.gz | head -n 100 | grep "#CHROM" | cut -f 10- | tr '\t' '\n') \
    | awk '$1!=$2'
zcat gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz /gscmnt/gc2802/halllab/ajscott/gtex_v7/melt_mei/gtex.631.melt_mei.assess_5.platform.vcf.gz \
        | grep -v "#" \
        | cat header.txt - \
        | vcfsort \
        | bgzip -c \
        > gtex.lumpy.gs.melt.1KGP_0.5.size.high_conf.vcf.gz
tabix -p vcf -f gtex.lumpy.gs.melt.1KGP_0.5.size.high_conf.vcf.gz

zcat gtex.lumpy.gs.1KGP_0.5.size.no_redundancy.high_conf.vcf.gz /gscmnt/gc2802/halllab/ajscott/gtex_v7/melt_mei/gtex.631.melt_mei.assess_5.platform.vcf.gz \
        | grep -v "#" \
        | cat header.txt - \
        | vcfsort \
        | bgzip -c \
        > gtex.lumpy.gs.melt.1KGP_0.5.size.no_redundancy.high_conf.vcf.gz
tabix -p vcf -f gtex.lumpy.gs.melt.1KGP_0.5.size.no_redundancy.high_conf.vcf.gz




# count number of variants per sample
mkdir -p sv_count
mkdir -p sv_count/logs


# lumpy counts
while read SAMPLE LONG BAM LUMPDIR
do
    bomb -g /ajscott/sv_aggregate -m 2 -J $SAMPLE.svcount -q hall-lab -o sv_count/logs/$SAMPLE.high.log -e sv_count/logs/$SAMPLE.high.log \
        "bash sv_counts.sh \
            gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz $SAMPLE \
            > sv_count/$SAMPLE.high.count.txt"
done < sample.map

# consolidate into a single file
while read SAMPLE LONG BAM LUMPDIR
do
    cat sv_count/$SAMPLE.high.count.txt
done < sample.map | awk '{ print $1,$2,$4+$5 }' OFS="\t" | groupBy -g 1 -c 3 -o collapse | tr ',' '\t' > sv_count/post-merged.sv_counts.high.txt



# gs counts
while read SAMPLE LONG BAM LUMPDIR
do
    bomb -g /ajscott/sv_aggregate -q hall-lab -m 1 -J $SAMPLE.svcount -e sv_count/logs/$SAMPLE.gs.redundant.high.log -o sv_count/logs/$SAMPLE.gs.redundant.high.log "zcat gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz | vawk '{if(I\$SVTYPE==\"CNV\"){if(S\$$SAMPLE\$CN != I\$MODECN){print I\$GSCNCATEGORY}}}' | sort | uniq -c > sv_count/$SAMPLE.gs.high.count.txt"
done < sample.map

while read SAMPLE LONG BAM LUMPDIR
do
    bomb -g /ajscott/sv_aggregate -q hall-lab -m 1 -J $SAMPLE.svcount -e sv_count/logs/$SAMPLE.gs.high.log -o sv_count/logs/$SAMPLE.gs.high.log "zcat gtex.lumpy.gs.1KGP_0.5.size.non_redundant.high_conf.vcf.gz | vawk '{if(I\$SVTYPE==\"CNV\"){if(S\$$SAMPLE\$CN != I\$MODECN){print I\$GSCNCATEGORY}}}' | sort | uniq -c > sv_count/$SAMPLE.gs.non_redundant.high.count.txt"
done < sample.map

while read SAMPLE LONG BAM LUMPDIR
do
    awk -v sample=$SAMPLE '{OFS="\t"; print sample,$2,$1}' sv_count/$SAMPLE.gs.high.count.txt
done < sample.map | groupBy -g 1 -c 3 -o collapse | tr ',' '\t' > sv_count/gs.high.no_redundancy.sv_counts.txt


# while read SAMPLE LONG BAM LUMPDIR
# do
#     bomb -g /ajscott/sv_aggregate -q hall-lab -m 1 -J $SAMPLE.svcount -e sv_count/logs/$SAMPLE.gs.low.log -o sv_count/logs/$SAMPLE.gs.low.log "zcat gtex.lumpy.gs.1KGP_0.5.size.female_af.low_pass.vcf.gz | vawk '{if(I\$SVTYPE==\"CNV\"){if(S\$$SAMPLE\$CN != I\$MODECN){print I\$GSCNCATEGORY}}}' | sort | uniq -c > sv_count/$SAMPLE.gs.low.count.txt"
# done < sample.map

# while read SAMPLE LONG BAM LUMPDIR
# do
#     awk -v sample=$SAMPLE '{OFS="\t"; print sample,$2,$1}' sv_count/$SAMPLE.gs.low.count.txt
# done < sample.map | groupBy -g 1 -c 3 -o collapse | tr ',' '\t' > sv_count/gs.low.sv_counts.txt

# get length and allele frequency data for qc plots
zcat gtex.lumpy.gs.melt.1KGP_0.5.size.non_redundant.high_conf.vcf.gz \
    | vcfToBed \
    | vawk -c 9 '{if($4~/^LUMPY/){print $1,$2,$3,$4,I$SVTYPE,I$AF,I$SVLEN} else{print $1,$2,$3,$4,"GS_"I$GSCNCATEGORY,I$AF,I$GCLENGTH}}' \
    | awk '{OFS="\t"; if($4~/^LUMPY_BND/){print $1,$2,$3,$4,$5,$6,"1"} else if($7<0){print $1,$2,$3,$4,$5,$6,$7*-1} else{print $0}}' \
    | bgzip -c \
    > gtex.lumpy.gs.1KGP_0.5.size.non_redundant.high_conf.bed.gz

zcat gtex.lumpy.gs.1KGP_0.5.size.high_conf.vcf.gz \
    | vcfToBed \
    | vawk -c 9 '{if($4~/^LUMPY/){print $1,$2,$3,$4,I$SVTYPE,I$AF,I$SVLEN} else{print $1,$2,$3,$4,"GS_"I$GSCNCATEGORY,I$AF,I$GCLENGTH}}' \
    | awk '{OFS="\t"; if($4~/^LUMPY_BND/){print $1,$2,$3,$4,$5,$6,"1"} else if($7<0){print $1,$2,$3,$4,$5,$6,$7*-1} else{print $0}}' \
    | bgzip -c \
    > gtex.lumpy.gs.1KGP_0.5.size.high_conf.bed.gz


# # qual filtered gs counts
# while read SAMPLE LONG BAM LUMPDIR
# do
#     bomb -g /ajscott/sv_aggregate -q hall-lab -m 1 -J $SAMPLE.svcount -e sv_count/logs/$SAMPLE.gs.high.log -o sv_count/logs/$SAMPLE.gs.high.log "zcat gtex.lumpy.gs.1KGP_0.5.size.gscnqual.high_conf.vcf.gz | vawk '{if(I\$SVTYPE==\"CNV\"){if(S\$$SAMPLE\$CN != I\$MODECN){print I\$GSCNCATEGORY}}}' | sort | uniq -c > sv_count/$SAMPLE.gs.qual_24.high.count.txt"
# done < sample.map

# while read SAMPLE LONG BAM LUMPDIR
# do
#     awk -v sample=$SAMPLE '{OFS="\t"; print sample,$2,$1}' sv_count/$SAMPLE.gs.qual_24.high.count.txt
# done < sample.map | groupBy -g 1 -c 3 -o collapse | tr ',' '\t' > sv_count/gs.high.qual_24.sv_counts.txt


# # get length and allele frequency data for qc plots (qual filtered)
# zcat gtex.lumpy.gs.1KGP_0.5.size.gscnqual.high_conf.vcf.gz \
#     | vcfToBed \
#     | vawk -c 9 '{if($4~/^LUMPY/){print $1,$2,$3,$4,I$SVTYPE,I$AF,I$SVLEN} else{print $1,$2,$3,$4,"GS_"I$GSCNCATEGORY,I$AF,I$GCLENGTH}}' \
#     | awk '{OFS="\t"; if($4~/^LUMPY_BND/){print $1,$2,$3,$4,$5,$6,"1"} else if($7<0){print $1,$2,$3,$4,$5,$6,$7*-1} else{print $0}}' \
#     | bgzip -c \
#     > gtex.lumpy.gs.1KGP_0.5.size.gscnqual.high_conf.bed.gz






# -------------------------------------------
# calculate some statistics on this data set

# 1. original Genome STRiP distr:
zcat $GS | vawk '{ print I$GSCNCATEGORY }' | sort | uniq -c | awk '{ print $2,$1 }' OFS="\t"
# DEL 27455
# DUP 6893
# MIXED   18398

# total: 52746

# 2. Genome strip calls that were merged in (don't count these, were merged in files below):
zcat gtex.lumpy.gs.low_pass.vcf.gz \
    | vawk '{ if (I$GSMERGE!="") print I$GSMERGE }' \
    | tr ',' '\n' \
    | zapdups -u \
    | wc -l
# 16042

# GS variants merged from 0.5 reciprocal overlaps
wc -l combine_lumpy_gs_2017-01-26/recip_f_0.5.txt
# 9829 combine_lumpy_gs_2017-01-26/recip_f_0.5.txt

# GS variants merged from containment
cat combine_lumpy_gs_2017-01-26/gs.contained.f_0.9.r_sq_0.25.collapse.txt | cut -f 2 | tr ',' '\n' | zapdups -u | wc -l
# 5317

# Get unmerged GS calls
# check that this adds up for GS variants. 37155+9780+5283 = 552746. check.
zcat gtex.lumpy.gs.low_pass.vcf.gz | vawk '$3~"^GS"' | wc -l
# 37600



# check that it adds up for LUMPY variants. check.
zcat gtex.lumpy.gs.low_pass.vcf.gz | vawk '$3~"^LUMPY" && !I$SECONDARY' | wc -l
# 74242

zcat $LUMPY | vawk '$3~"^LUMPY" && !I$SECONDARY' | wc -l
# 74242




# get their type distribution of merged variants (contained)
zjoin -a $GS -b <(cat combine_lumpy_gs_2017-01-26/gs.contained.f_0.9.r_sq_0.25.collapse.txt | cut -f 2 | tr ',' '\n') -w a -1 3 -2 1 -p "#" \
    | vawk '{ print I$GSCNCATEGORY }' | sort | uniq -c | awk '{ print $2,$1 }' OFS="\t"
# DEL 1711
# DUP 947
# MIXED   2659

# (recip 0.5)
zjoin -a $GS -b combine_lumpy_gs_2017-01-26/recip_f_0.5.txt -w a -1 3 -2 2 -p "#" \
    | vawk '{ print I$GSCNCATEGORY }' | sort | uniq -c | awk '{ print $2,$1 }' OFS="\t"
# DEL 6451
# DUP 772
# MIXED   2606

# Classes of both 0.5 recipr overlap and 0.9 containment combined.
zjoin -a $GS -b <(cat combine_lumpy_gs_2017-01-26/recip_f_0.5.txt combine_lumpy_gs_2017-01-26/gs.contained.f_0.9.r_sq_0.25.collapse.txt | cut -f 2 | tr ',' '\n') -w a -1 3 -2 1 -p "#" \
    | vawk '{ print I$GSCNCATEGORY }' | sort | uniq -c | awk '{ print $2,$1 }' OFS="\t"
# DEL 8162
# DUP 1719
# MIXED   5265
