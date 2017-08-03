# 2016-08-30

pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28

# fix header 
 bomb -q long -m 1 \
    "zcat CNVdiscovery/run/results/gs_cnv.genotypes.vcf.gz \
        | vawk '{print $0}' \
        | cat gs.header.trunc.txt - \
        | bgzip -c \
        > gs_cnv.genotypes.vcf.gz"

# covariance table for platform bias filtering: /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2016-08-11/gtex_wgs_covar.txt

# exclude platform bias variants, add MODECN 
# annotate the filtering of the variants
zcat /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gs_cnv.genotypes.vcf.gz \
    | awk '{ if ($0~"^#") { print ; next } delete ALLELE; NUM_SAMPLES=NF-9; for (i=10;i<=NF;++i) { split($i,gt,":"); ALLELE[gt[2]]++; } MODE=0; MODE_ALLELE=0; for (a in ALLELE) { if (ALLELE[a] > MODE) { MODE=ALLELE[a]; MODE_ALLELE=a } } if (ALLELE[MODE_ALLELE]<NUM_SAMPLES) { $8=$8";MODECN="MODE_ALLELE; print } }' OFS="\t" \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i MODECN -c INFO -n 1 -t Integer -d "Mode of the CN values across samples" \
    | awk '{ if ($0~"^#") { print; next; } $9=$9":DS"; split($9,fmt,":"); for (i=1;i<=length(fmt);++i) { if (fmt[i]=="CN") fmt_idx=i } for (i=10;i<=NF;++i) { split($i,gt,":"); $i=$i":"gt[fmt_idx]*2 } print }' OFS="\t" \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i DS -c FORMAT -n A -t Float -d "Dosage of alt allele" \
    | vawk --header '{ $3="GS_"I$GSCNCATEGORY"_"$3; $8=$8";CIPOS=-100,100;CIEND=-100,100"; print }' \
    | /gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/src/vcf_covar.py -c /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2016-08-11/gtex_wgs_covar.txt -f CN \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i "PLATFORM" -c FILTER -d "Variant genotypes correlate with sequencing platform with R^2 greater than 0.1" \
    | vawk --header '{ if (I$EXVAR>0.1) $7="PLATFORM"; print }' \
    | vawk --header '{ if ($7=="PASS") $7="."; print }' \
    | bgzip -c \
    > gtex.genome_strip.filt_annot.vcf.gz
tabix -f -p vcf gtex.genome_strip.filt_annot.vcf.gz


# get the private and doubleton variants that we want to merge
zcat gtex.genome_strip.filt_annot.vcf.gz \
    | vawk --header '$7=="."' \
    | vawk --header '{ if (I$SVTYPE=="CNV") { if (I$GSCNALLELES==2) print $1,$2,$3,$4,$5,$6,$7,$8,"CN",S$*$CN; } else print }' \
    | vawk --header '{ if (I$SVTYPE!="CNV") { print; next } delete ALLELE; for (i=10;i<=NF;++i) ALLELE[$i]++; MODE=0; MODE_ALLELE=0; for (a in ALLELE) { if (ALLELE[a] > MODE) { MODE=ALLELE[a]; MODE_ALLELE=a } } NSAMP=0; for (i=10;i<=NF;++i) { if ($i!=MODE_ALLELE) NSAMP++ } if (NSAMP<=2) { printf $1; for (i=2;i<=8;++i) printf "\t"$i; printf "\tGT"; for (i=10;i<=NF;++i) { if ($i==MODE_ALLELE) printf "\t0/0"; else printf "\t0/1" } printf "\n" } }' \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_allele_freq.py \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i POSSAMP -c INFO -t String -n "." -d "Positively genotyped samples" \
    | /gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/src/get_pos_samples.py \
    | bgzip -c \
    > gtex.genome_strip.private_doubleton.low_pass.vcf.gz
tabix -p vcf -f gtex.genome_strip.private_doubleton.low_pass.vcf.gz


pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/rare_merge/

# Generate file with long ID to make sure only calls in same samples and of same type are merged:
bomb -q hall-lab -m 1 \
    "zcat ../gtex.genome_strip.private_doubleton.low_pass.vcf.gz \
        | vcfToBed \
        | vawk -c 9 '{print \$1\"_\"I\$POSSAMP\"_\"I\$GSCNCATEGORY,\$2,\$3,\$4}' \
        | bgzip -c \
        > gtex.genome_strip.private_doubleton.low_pass.longID.bed.gz"

zcat gtex.genome_strip.private_doubleton.low_pass.longID.bed.gz  | bedtools sort > input.bed
cat input.bed |  bedtools merge -i stdin -d 10000000 -c 4 -o collapse > merge1.bed
bedtools coverage -a input.bed -b merge1.bed | awk '$8>=0.1'  > merge1.gt10.bed
cat merge1.gt10.bed | awk '$8>=0.1' | awk '{OFS="\t"; split($1,x,"_"); print x[1],$2,$3,$4}' > merge2.bed
cat input.bed | bedtools intersect -v -a stdin -b merge1.gt10.bed | awk '{OFS="\t"; split($1,x,"_"); print x[1],$2,$3,$4}' > unmerged.bed
cat merge2.bed unmerged.bed | bedtools sort > gtex.genome_strip.private_doubleton.low_pass.longID.merged.bed

# add the merged rare variants and remove their constituent fragments
# First, for each new variant, grab the first old variant and modify its position.
VCF_NCOL=`zcat ../gtex.genome_strip.filt_annot.vcf.gz | grep -v "^#" | head -n 1 | awk '{ print NF }'`

zcat ../gtex.genome_strip.filt_annot.vcf.gz  \
    | grep -v "^#" \
    | zjoin -a stdin -b <(cat gtex.genome_strip.private_doubleton.low_pass.longID.merged.bed | awk '$4~","' | awk '{ IDX=$4; gsub(",.*","",IDX); print $0,IDX }' OFS="\t") -1 3 -2 5 \
    | vawk '{ $2=$(NF-3); gsub(";END=[^;]*;",";END="$(NF-2)";",$8); split($3,ID_ARR,"_"); $3=ID_ARR[1]"_"ID_ARR[2]"_"ID_ARR[3]"_"ID_ARR[4]"_"$(NF-3)"_"$(NF-2); $8=$8";GSMERGE="$(NF-1); print }' \
    | cut -f -${VCF_NCOL} \
    > merged.gs.calls.txt

# remove the old calls from our VCF file and replace with new ones
cat gtex.genome_strip.private_doubleton.low_pass.longID.merged.bed | awk '$4~","' \
    | cut -f 4 | tr ',' '\n' \
    | zjoin -v -a ../gtex.genome_strip.filt_annot.vcf.gz -b stdin -1 3 -2 1 -p "#" \
    | cat - merged.gs.calls.txt \
    | vcfsort \
    | bgzip -c \
    > ../gtex.genome_strip.filt_annot.merged.vcf.gz
tabix -p vcf -f ../gtex.genome_strip.filt_annot.merged.vcf.gz

# mark pseudoautosomal variants as filtered
cd /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/

NCOL=`zcat gtex.genome_strip.filt_annot.merged.vcf.gz | grep -v "^#" | head -n 1 | awk '{ print NF }'`
tabix -h gtex.genome_strip.filt_annot.merged.vcf.gz X Y | vcfToBed | bedtools intersect -a stdin -b /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/par.bed.gz \
    | cut -f 4 \
    | zjoin -r -a gtex.genome_strip.filt_annot.merged.vcf.gz -b stdin -1 3 -2 1 -p "#" \
    | vawk --header '{ if ($NF!="NA") $7="PAR"; print }' \
    | cut -f -$NCOL \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i "PAR" -c FILTER -d "Low frequency DEL or DUP in the pseudoautosomal region of X and Y" \
    | bgzip -c \
    > tmp.vcf.gz &&
mv tmp.vcf.gz gtex.genome_strip.filt_annot.merged.vcf.gz &&
tabix -p vcf -f gtex.genome_strip.filt_annot.merged.vcf.gz


# get the low pass variants
zcat gtex.genome_strip.filt_annot.merged.vcf.gz \
    | vawk --header '$7=="."' \
    | bgzip -c \
    > gtex.genome_strip.filt_annot.merged.low_pass.vcf.gz
tabix -p vcf -f gtex.genome_strip.filt_annot.merged.low_pass.vcf.gz















