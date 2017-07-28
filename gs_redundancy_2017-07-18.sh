# 2017-07-18

# collapse redundant GS calls post-merge with lumpy, but before filtering based on GSCNQUAL

pwd	
# /gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/gs_redundancy/

 zcat /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gtex.genome_strip.filt_annot.merged.low_pass.vcf.gz \
 	| vawk '{if($3~/^GS/){print $1,$2,I$END,$3}}' \
 	| bgzip -c \
 	> gs.1KGP_0.5.bed.gz

WINDOW_SIZE=10000
OVERLAP=$((WINDOW_SIZE/2))
LINKAGE=0.9

# make genome-wide windows
bedtools makewindows -g /gscmnt/gc2719/halllab/genomes/human/GRCh37/1kg_phase1/human_g1k_v37.genome -w $WINDOW_SIZE -s $OVERLAP \
	| grep -v MT \
	| grep -v GL \
	> b37.windows.bed

# get GS calls in each window
bedtools window -a b37.windows.bed -b gs.1KGP_0.5.bed.gz -w 0 \
	| bedtools groupby -g 1-3 -c 7 -o distinct \
	| cut -f 4 \
	| grep "," \
	| uniq \
	> overlapping_window.clusters.txt
/gscmnt/gc2719/halllab/users/ajscott/src/collapse_groups.py -i overlapping_window.clusters.txt > window_groups.collapsed.txt

# check that each variant only in 1 cluster. check.
for VAR in `cat window_groups.collapsed.txt | tr ',' '\n'`
do
	COUNT=`grep -c $VAR window_groups.collapsed.txt`
	echo -e $VAR"\t"$COUNT
done | awk '$2!=1'

# make list of all pairwise comparisions in each cluster
cat window_groups.collapsed.txt \
	| awk '{OFS="\t"; split($1,a,","); if(length(a)>1){for(x=1;x<=length(a);x++){for(y=1;y<=length(a);y++){print a[x],a[y]}}}}' \
	> window_groups.pairs.txt

# get linkage between each pair
/gscmnt/gc2802/halllab/cchiang/projects.2719/gtex/src/var_gt_corr.py \
    -a /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gtex.genome_strip.filt_annot.merged.low_pass.vcf.gz \
    -b /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gtex.genome_strip.filt_annot.merged.low_pass.vcf.gz \
    -v window_groups.pairs.txt \
    -af CN -bf CN \
    | awk '{OFS="\t"; print $1,$2,$3**2}' \
    > window_groups.linkage.txt

/gscmnt/gc2719/halllab/users/ajscott/src/gs_cluster_linkage.py -i window_groups.collapsed.txt -l window_groups.linkage.txt -r2 $LINKAGE \
	| sort \
	| uniq \
	> window_clusters.$LINKAGE.txt




# make tracks to look at merged variants in IGV
cat window_clusters.$LINKAGE.txt | tr ',' '\n' | zjoin -a stdin -b gs.1KGP_0.5.size.high_conf.bed.gz -1 1 -2 4 -w b > gs.merged.$LINKAGE.bed

# use this python code to get the span of each cluster:
# #!/usr/bin/env python

# infile = open("/gscmnt/gc2802/halllab/gtex_2016-06-01/merge_2017-01-17/gs_redundancy/window_clusters.0.9.txt", 'r')

# for line in infile:
# 	cluster = line.rstrip('\n').split(',')
# 	chrom = cluster[0].split("_")[3]
# 	start = cluster[0].split("_")[4]
# 	end = cluster[0].split("_")[5]
# 	for x in range(1,len(cluster)):
# 		if cluster[x].split("_")[4] < start:
# 			start = cluster[x].split("_")[4]
# 		if cluster[x].split("_")[5] > end:
# 			end = cluster[x].split("_")[5]
# 	output = chrom + "\t" + start + "\t" + end + "\t" + line.rstrip('\n')
# 	print(output)

# infile.close()

python temp.py > cluster_span.0.9.bed




# only want to merge calls that are overlapping (redundant)
# assign cluster number to each variant 
cut -f 4 cluster_span.0.9.bed \
	| awk '{OFS="\t"; split($1,a,","); for(x=1;x<=length(a);x++){print NR,a[x]}}' \
	| awk '{OFS="\t"; split($2,a,"_"); print $1,a[4],a[5],a[6],$2}' \
	> cluster_numbers.txt


# get clustered variants that overlap
for x in {1..1019}
do
	cat cluster_numbers.txt | awk -v x=$x '$1==x' | cut -f 2- > tmp.txt
	bedtools intersect -a tmp.txt -b tmp.txt -wao | awk '$4!=$8' >> overlapping.txt
	rm -rf tmp.txt 
done
awk '{OFS="\t"; print $1,$2,$3,$4"\n"$5,$6,$7,$8}' overlapping.txt | sort -k4,4 | uniq > overlapping.bed

#re-generate groups, just with overlapping variants
zjoin -a cluster_numbers.txt -b overlapping.bed -1 5 -2 4 -w a > to_collapse.txt
awk '{OFS="\t"; print $1,$4-$3,$5}'  to_collapse.txt | groupBy -g 1 -c 2,3 -o collapse,collapse > groups_sizes.txt

# print out largest variant as one to keep
awk '{OFS="\t"; split($2,a,","); split($3,b,","); max=a[1]; i=1; for(x=2;x<=length(a);x++){if(a[x]>max){max=a[x]; i=x}}; print b[i]}' groups_sizes.txt > keep.txt

# remove smaller overlapping variants
zjoin -a keep.txt -b to_collapse.txt -1 1 -2 5 -V | cut -f 5 > remove.txt

# also keep non-overlapping variants
cat keep.txt remove.txt | zjoin -a stdin -b gs.merged.0.9.bed -1 1 -2 4 -V | cut -f 4 >> keep.txt


# remove from VCF file
zjoin -a /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gtex.genome_strip.filt_annot.merged.low_pass.vcf.gz -b remove.txt -1 3 -2 1 -p "#" -w a -v \
	| bgzip -c \
	> /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28/gtex.genome_strip.filt_annot.merged.no_redundancy.low_pass.vcf.gz










