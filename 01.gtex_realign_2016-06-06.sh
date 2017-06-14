# 2016-06-06

# ========================================
# GTEx batch realignment from BAM
# ========================================

pwd 
# /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06


# ========================================
# 1) Create the sample map (batch.1)
# ========================================

# batch 1: 62 samples

# make the batch
mkdir -p notes/batch1

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1*/*/*/*.bam`
do
	BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
	echo -e "$BAM\t$BAMPATH"
done | awk '{OFS="\t"; split($1,a,"."); print a[2]"_"a[3],$2}' >> notes/batch1/batch1.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch1/batch1.txt
head $BATCH
# GTEX-WHWD-0002_4        /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1035926/SRX1151833/SRR2165094_processed/G67881.GTEX-WHWD-0002.4.bam
# GTEX-14PHY-0002_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1396605/SRX1704483/SRR3382609_processed/G80551.GTEX-14PHY-0002.1.bam
# GTEX-1GMR3-0001_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1396609/SRX1704489/SRR3382626_processed/G87945.GTEX-1GMR3-0001.1.bam
# GTEX-1GMRU-0004_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1396772/SRX1704678/SRR3382878_processed/G87945.GTEX-1GMRU-0004.1.bam
# GTEX-15CHS-0001_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1396956/SRX1704851/SRR3383188_processed/G80551.GTEX-15CHS-0001.1.bam
# GTEX-14A6H-0003_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1396990/SRX1704880/SRR3383262_processed/G80551.GTEX-14A6H-0003.1.bam
# GTEX-XMD3-0003_1        /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1397529/SRX1705196/SRR3383940_processed/G80551.GTEX-XMD3-0003.1.bam
# GTEX-14PJ6-0001_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1397773/SRX1705469/SRR3384403_processed/G80551.GTEX-14PJ6-0001.1.bam
# GTEX-15ER7-0002_8       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1398174/SRX1706071/SRR3386142_processed/G80551.GTEX-15ER7-0002.8.bam
# GTEX-12C56-0003_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50512/reads-src/SRP012682/SRS1398176/SRX1706076/SRR3386165_processed/G80551.GTEX-12C56-0003.1.bam

###################################################################
# batch 2: 54 samples

mkdir -p notes/batch2

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch2/batch2.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch2/batch2.txt
head $BATCH
# GTEX-13W46-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423811/SRX1744729/SRR3478796_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597874/GTEX-13W46-0003.bam
# GTEX-13NZ8-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423812/SRX1744730/SRR3478797_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597473/GTEX-13NZ8-0004.bam
# GTEX-113JC-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423813/SRX1744732/SRR3478798_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596625/GTEX-113JC-0001.bam
# GTEX-WH7G-0004  /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423828/SRX1744750/SRR3478896_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596783/GTEX-WH7G-0004.bam
# GTEX-15SDE-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423829/SRX1744753/SRR3478900_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596132/GTEX-15SDE-0001.bam
# GTEX-13SLX-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423830/SRX1744754/SRR3478901_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595672/GTEX-13SLX-0003.bam
# GTEX-ZLFU-0004  /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423831/SRX1744755/SRR3478902_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597207/GTEX-ZLFU-0004.bam
# GTEX-1B996-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423832/SRX1744756/SRR3478903_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597218/GTEX-1B996-0003.bam
# GTEX-13CF3-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423833/SRX1744758/SRR3478906_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597965/GTEX-13CF3-0004.bam
# GTEX-1313W-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50517/reads-src/SRP012682/SRS1423834/SRX1744759/SRR3478907_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597635/GTEX-1313W-0004.bam

###################################################################
# batch 3: 51 samples

mkdir -p notes/batch3

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch3/batch3.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch3/batch3.txt
head $BATCH
# GTEX-15DYW-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425127/SRX1746422/SRR3481584_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596542/GTEX-15DYW-0002.bam
# GTEX-117XS-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425128/SRX1746423/SRR3481585_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598017/GTEX-117XS-0001.bam
# GTEX-13FHP-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425130/SRX1746425/SRR3481586_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598497/GTEX-13FHP-0004.bam
# GTEX-15ETS-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425131/SRX1746446/SRR3481640_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596253/GTEX-15ETS-0003.bam
# GTEX-12ZZX-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425132/SRX1746427/SRR3481588_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598707/GTEX-12ZZX-0001.bam
# GTEX-QMR6-1926  /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425133/SRX1746428/SRR3481590_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597213/GTEX-QMR6-1926.bam
# GTEX-13X6H-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425134/SRX1746429/SRR3481591_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597620/GTEX-13X6H-0001.bam
# GTEX-1399Q-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425135/SRX1746430/SRR3481592_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596297/GTEX-1399Q-0004.bam
# GTEX-12584-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425137/SRX1746437/SRR3481593_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598836/GTEX-12584-0004.bam
# GTEX-1GF9V-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50519/reads-src/SRP012682/SRS1425139/SRX1746441/SRR3481635_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598835/GTEX-1GF9V-0004.bam

###################################################################
# batch 4: 5 samples

mkdir -p notes/batch4

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50521/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch4/batch4.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch4/batch4.txt
head $BATCH
# GTEX-16NFA-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50521/reads-src/SRP012682/SRS1429226/SRX1751595/SRR3489517_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595983/GTEX-16NFA-0003.bam
# GTEX-145MH-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50521/reads-src/SRP012682/SRS1429228/SRX1751598/SRR3489806_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598815/GTEX-145MH-0001.bam
# GTEX-1GN1V-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50521/reads-src/SRP012682/SRS1429229/SRX1751600/SRR3489529_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596823/GTEX-1GN1V-0002.bam
# GTEX-14PJN-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50521/reads-src/SRP012682/SRS1430368/SRX1753490/SRR3493409_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596040/GTEX-14PJN-0002.bam
# GTEX-X4EP-0002  /gscmnt/gc2802/halllab/gtex_2016-06-01/50521/reads-src/SRP012682/SRS1430369/SRX1753491/SRR3493434_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598927/GTEX-X4EP-0002.bam

###################################################################
# batch 5: 67 samples

mkdir -p notes/batch5

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch5/batch5.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch5/batch5.txt
head $BATCH
# GTEX-ZYWO-0004  /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424220/SRX1745231/SRR3479746_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597420/GTEX-ZYWO-0004.bam
# GTEX-147JS-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424221/SRX1745230/SRR3479745_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598404/GTEX-147JS-0003.bam
# GTEX-13QJ3-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424222/SRX1745232/SRR3479747_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598892/GTEX-13QJ3-0002.bam
# GTEX-11NV4-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424223/SRX1745233/SRR3479748_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595970/GTEX-11NV4-0001.bam
# GTEX-13W3W-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424224/SRX1745234/SRR3479749_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597138/GTEX-13W3W-0004.bam
# GTEX-1F7RK-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424225/SRX1745235/SRR3479750_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597689/GTEX-1F7RK-0004.bam
# GTEX-14PKV-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424229/SRX1745241/SRR3479753_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596177/GTEX-14PKV-0002.bam
# GTEX-15G19-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424230/SRX1745242/SRR3479754_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598679/GTEX-15G19-0001.bam
# GTEX-ZQUD-0002  /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424243/SRX1745255/SRR3479774_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597357/GTEX-ZQUD-0002.bam
# GTEX-16Z82-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50518/reads-src/SRP012682/SRS1424244/SRX1745256/SRR3479775_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595777/GTEX-16Z82-0004.bam

###################################################################
# batch 6: 53 samples

mkdir -p notes/batch6

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch6/batch6.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch6/batch6.txt
head $BATCH
# GTEX-11ZUS-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423785/SRX1744692/SRR3478750_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596742/GTEX-11ZUS-0003.bam
# GTEX-131XF-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423787/SRX1766345/SRR3529270_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598744.1/GTEX-131XF-0004.bam
# GTEX-17HGU-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423788/SRX1744695/SRR3478752_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596105/GTEX-17HGU-0001.bam
# GTEX-13NYC-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423793/SRX1744703/SRR3478761_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596765/GTEX-13NYC-0004.bam
# GTEX-16BQI-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423795/SRX1744705/SRR3478763_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596309/GTEX-16BQI-0003.bam
# GTEX-15DDE-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423800/SRX1744711/SRR3478770_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598755/GTEX-15DDE-0001.bam
# GTEX-16NPV-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423804/SRX1744720/SRR3478781_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596965/GTEX-16NPV-0003.bam
# GTEX-13S86-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423809/SRX1744726/SRR3478793_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598989/GTEX-13S86-0003.bam
# GTEX-16YQH-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423810/SRX1744728/SRR3478795_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596518/GTEX-16YQH-0004.bam
# GTEX-13OVL-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50523/reads-src/SRP012682/SRS1423842/SRX1744769/SRR3478918_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598018/GTEX-13OVL-0003.bam

###################################################################
# batch 7: 61 samples

mkdir -p notes/batch7

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS14*/*/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done | awk '{OFS="\t"; split($1,a,"."); print a[2]"_"a[3],$2}' > notes/batch7/batch7.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch7/batch7.txt
head $BATCH
# GTEX-1A8G7-0003_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1404370/SRX1715101/SRR3404500_processed/G87945.GTEX-1A8G7-0003.1.bam
# GTEX-ZYFD-0001_1        /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1405703/SRX1716836/SRR3407178_processed/G77324.GTEX-ZYFD-0001.1.bam
# GTEX-15RJE-0003_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1411637/SRX1726035/SRR3439010_processed/G76371.GTEX-15RJE-0003.1.bam
# GTEX-ZXG5-0002_1        /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1411892/SRX1726325/SRR3439556_processed/G77324.GTEX-ZXG5-0002.1.bam
# GTEX-U8T8-2026_1        /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1411982/SRX1726458/SRR3439989_processed/G80551.GTEX-U8T8-2026.1.bam
# GTEX-1CAMQ-0003_2       /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1411995/SRX1726481/SRR3440041_processed/G87945.GTEX-1CAMQ-0003.2.bam
# GTEX-17KNJ-0003_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1412005/SRX1726491/SRR3440048_processed/G76371.GTEX-17KNJ-0003.1.bam
# GTEX-XMD2-0001_1        /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1412184/SRX1726780/SRR3440345_processed/G80551.GTEX-XMD2-0001.1.bam
# GTEX-1B97I-0004_2       /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1412763/SRX1727778/SRR3441539_processed/G87945.GTEX-1B97I-0004.2.bam
# GTEX-16A39-0004_1       /gscmnt/gc2802/halllab/gtex_2016-06-01/50516/reads-src/SRP012682/SRS1412764/SRX1727779/SRR3441544_processed/G76371.GTEX-16A39-0004.1.bam

###################################################################
# batch 8: 48 samples

mkdir -p notes/batch8

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch8/batch8.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch8/batch8.txt
head $BATCH
# GTEX-1C4CL-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424896/SRX1746063/SRR3480723_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598403/GTEX-1C4CL-0003.bam
# GTEX-XQ8I-0001  /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424955/SRX1746207/SRR3481286_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04599088/GTEX-XQ8I-0001.bam
# GTEX-16NPX-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424956/SRX1746208/SRR3481287_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596304/GTEX-16NPX-0002.bam
# GTEX-1GZHY-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424958/SRX1746210/SRR3481289_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597054/GTEX-1GZHY-0001.bam
# GTEX-1B97J-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424959/SRX1746211/SRR3481290_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595468/GTEX-1B97J-0003.bam
# GTEX-S33H-0004  /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424960/SRX1746215/SRR3481295_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597750/GTEX-S33H-0004.bam
# GTEX-1C2JI-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424961/SRX1746217/SRR3481297_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597723/GTEX-1C2JI-0003.bam
# GTEX-11XUK-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424962/SRX1746220/SRR3481300_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598124/GTEX-11XUK-0003.bam
# GTEX-14PII-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424963/SRX1746221/SRR3481301_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597266/GTEX-14PII-0004.bam
# GTEX-1C64N-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50520/reads-src/SRP012682/SRS1424966/SRX1746226/SRR3481305_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598190/GTEX-1C64N-0003.bam

###################################################################
# batch 9: 45 samples

mkdir -p notes/batch9

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch9/batch9.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch9/batch9.txt
head $BATCH
# GTEX-13112-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1423639/SRX1745990/SRR3480642_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596075/GTEX-13112-0004.bam
# GTEX-1A3MX-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1423748/SRX1744647/SRR3478691_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597568/GTEX-1A3MX-0002.bam
# GTEX-14E7W-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1423749/SRX1744649/SRR3478694_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596730/GTEX-14E7W-0001.bam
# GTEX-15SKB-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1423844/SRX1744771/SRR3478920_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04599136/GTEX-15SKB-0002.bam
# GTEX-1GF9U-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1423848/SRX1744778/SRR3478942_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597303/GTEX-1GF9U-0004.bam
# GTEX-1GMR8-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1424043/SRX1744998/SRR3479270_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595808/GTEX-1GMR8-0003.bam
# GTEX-YFC4-0001  /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1424044/SRX1744999/SRR3479271_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598775/GTEX-YFC4-0001.bam
# GTEX-YFCO-0001  /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1424151/SRX1745127/SRR3479630_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595607/GTEX-YFCO-0001.bam
# GTEX-12126-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1424152/SRX1745128/SRR3479631_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04599004/GTEX-12126-0001.bam
# GTEX-1F5PL-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50685/reads-src/SRP012682/SRS1424160/SRX1745138/SRR3479637_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596511/GTEX-1F5PL-0002.bam

###################################################################
# batch 10: 10 samples

mkdir -p notes/batch10

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch10/batch10.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch10/batch10.txt
head $BATCH
# GTEX-1C6VR-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424861/SRX1746023/SRR3480675_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598325/GTEX-1C6VR-0003.bam
# GTEX-178AV-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424862/SRX1746025/SRR3480678_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596989/GTEX-178AV-0004.bam
# GTEX-1B8SG-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424863/SRX1746027/SRR3480680_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596611/GTEX-1B8SG-0003.bam
# GTEX-11P81-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424864/SRX1746028/SRR3480681_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598309/GTEX-11P81-0002.bam
# GTEX-14PJ5-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424865/SRX1746031/SRR3480690_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596879/GTEX-14PJ5-0001.bam
# GTEX-1GL5R-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424866/SRX1746032/SRR3480691_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597446/GTEX-1GL5R-0002.bam
# GTEX-13CZV-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424867/SRX1746034/SRR3480694_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597836/GTEX-13CZV-0003.bam
# GTEX-15SB6-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424893/SRX1746059/SRR3480717_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597699/GTEX-15SB6-0002.bam
# GTEX-1AX9K-0001 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424894/SRX1746060/SRR3480718_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596912/GTEX-1AX9K-0001.bam
# GTEX-1A3MW-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50862/reads-src/SRP012682/SRS1424895/SRX1746061/SRR3480719_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598052/GTEX-1A3MW-0004.bam

###################################################################
# batch 11: 46 samples

mkdir -p notes/batch11

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch11/batch11.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch11/batch11.txt
head $BATCH
# GTEX-ZVP2-0002  /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1418168/SRX1766473/SRR3529803_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597939.1/GTEX-ZVP2-0002.bam
# GTEX-Y8E4-0001  /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1419033/SRX1766563/SRR3530190_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595641.1/GTEX-Y8E4-0001.bam
# GTEX-11WQK-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423469/SRX1766455/SRR3529722_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597290.1/GTEX-11WQK-0003.bam
# GTEX-146FR-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423574/SRX1744419/SRR3478276_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04597091/GTEX-146FR-0004.bam
# GTEX-RVPU-0001  /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423751/SRX1744651/SRR3478695_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598168/GTEX-RVPU-0001.bam
# GTEX-14ABY-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423753/SRX1744653/SRR3478696_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598420/GTEX-14ABY-0003.bam
# GTEX-1F88F-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423754/SRX1744656/SRR3478700_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598432/GTEX-1F88F-0003.bam
# GTEX-ZAB5-0004  /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423755/SRX1744657/SRR3478701_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598458/GTEX-ZAB5-0004.bam
# GTEX-14LLW-0002 /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423756/SRX1744658/SRR3478702_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04595625/GTEX-14LLW-0002.bam
# GTEX-S3XE-0003  /gscmnt/gc2802/halllab/gtex_2016-06-01/50686/reads-src/SRP012682/SRS1423757/SRX1744659/SRR3478703_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598401/GTEX-S3XE-0003.bam

###################################################################
# batch 12: 3 samples

mkdir -p notes/batch12

for BAMPATH in `ls /gscmnt/gc2802/halllab/gtex_2016-06-01/50970/reads-src/SRP012682/SRS1*/*/*/seq/submissions/epsilon9/prod/submissions/*/*.bam`
do
    BAM=`echo $BAMPATH | xargs -I{} basename {} .bam`
    echo -e "$BAM\t$BAMPATH"
done > notes/batch12/batch12.txt

BATCH=/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch12/batch12.txt
head $BATCH
# GTEX-1497J-0004 /gscmnt/gc2802/halllab/gtex_2016-06-01/50970/reads-src/SRP012682/SRS1423808/SRX1744727/SRR3478794_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598082/GTEX-1497J-0004.bam
# GTEX-13O21-0003 /gscmnt/gc2802/halllab/gtex_2016-06-01/50970/reads-src/SRP012682/SRS1425103/SRX1746385/SRR3481542_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04596449/GTEX-13O21-0003.bam
# GTEX-ZVT2-0003  /gscmnt/gc2802/halllab/gtex_2016-06-01/50970/reads-src/SRP012682/SRS1425552/SRX1747441/SRR3485450_processed/seq/submissions/epsilon9/prod/submissions/PRJNA75899_SAMN04598872/GTEX-ZVT2-0003.bam

# ----------------------------------------
# 2. Set up sample directories
# ----------------------------------------

# organize sample directory paths
for SAMPLE in `cat $BATCH | cut -f 1`
do
    mkdir -p /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE
    mkdir -p /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc
    mkdir -p /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log
done

# symlink the original BAM files
for SAMPLE in `cat $BATCH | cut -f 1`
do
    BAMPATH=`cat $BATCH | grep -m 1 $SAMPLE | cut -f 2`
    cd /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc
    ln -s $BAMPATH original.bam
    cd -
done

# ----------------------------------------
# 3. Realign the BAM files
# ----------------------------------------
# Use Haley's version of SpeedSeq for realignment, just to be safe

for SAMPLE in `cat $BATCH | cut -f 1`
do
    echo $SAMPLE
    bomb \
        -q long \
        -g /ajscott/gtex_realign \
        -m 48 \
        -t 8 \
        -J $SAMPLE.align \
        -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.align.%J.log \
        -e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.align.%J.log \
        "/gscmnt/sata849/info/speedseq_freeze/v3/speedseq/bin/speedseq realign \
            -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/$SAMPLE \
            -T /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/temp \
            -M 8 -n -t 8 \
            -v \
            /gscmnt/gc2719/halllab/genomes/human/GRCh37/hs37_ebv/hs37_ebv.fasta \
            /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/original.bam"
done

# align the 1 file that had to be downloaded as a fastq file
bomb \
    -q long \
    -g /ajscott/gtex_realign \
    -m 48 \
    -t 8 \
    -J GTEX-1GN2E-0003.align \
    -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/GTEX-1GN2E-0003/log/GTEX-1GN2E-0003.align.%J.log \
    -e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/GTEX-1GN2E-0003/log/GTEX-1GN2E-0003.align.%J.log \
    "/gscmnt/sata849/info/speedseq_freeze/v3/speedseq/bin/speedseq align \
        -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/GTEX-1GN2E-0003/GTEX-1GN2E-0003 \
        -R \"@RG\tID:GTEX-1GN2E-0003.S1\tPL:illumina\tLB:Solexa-315114\tSM:GTEX-1GN2E-0003\tCN:BI\" \
        -T /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/GTEX-1GN2E-0003/temp \
        -M 8 -t 8 \
        -v \
        /gscmnt/gc2719/halllab/genomes/human/GRCh37/hs37_ebv/hs37_ebv.fasta \
        /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/GTEX-1GN2E-0003/qc/SRR3538835_1.fq.gz \
        /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/GTEX-1GN2E-0003/qc/SRR3538835_2.fq.gz"

# ----------------------------------------
# 4. Flagstat the original BAM files
# ----------------------------------------

for SAMPLE in `cat $BATCH | cut -f 1`
do
    echo $SAMPLE
    bomb \
        -m 1 \
        -J $SAMPLE.flag \
        -g /ajscott/flagstat \
        -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.original.flagstat.%J.log \
        -e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.original.flagstat.%J.log \
        "samtools-1.1 flagstat /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/original.bam > /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/original.bam.flagstat"
done

# ----------------------------------------
# 5. Flagstat the realigned BAM files
# ----------------------------------------
for SAMPLE in `cat $BATCH | cut -f 1`
do
    echo $SAMPLE
    bomb \
        -m 1 \
        -J $SAMPLE.flag \
        -g /ajscott/flagstat \
        -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.realign.flagstat.%J.log \
        -e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.realign.flagstat.%J.log \
        "samtools-1.1 flagstat /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/$SAMPLE.bam > /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/$SAMPLE.bam.flagstat"
done

# ----------------------------------------
# 6. Compare the original and realigned flagstats
# ----------------------------------------

for SAMPLE in `cat $BATCH | cut -f 1`
do
    REALIGN_PAIRED=`cat /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/$SAMPLE.bam.flagstat | grep -m 1 "paired in sequencing" | awk '{ print $1+$3 }'`    
    ORIG_PAIRED=`cat /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/original.bam.flagstat | grep -m 1 "paired in sequencing" | awk '{ print $1+$3 }'`
    DIFF=$(( $REALIGN_PAIRED - $ORIG_PAIRED ))
    READ_LENGTH=`sambamba view /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/$SAMPLE.bam | head -n 10000 | awk 'BEGIN { MAX_LEN=0 } { LEN=length($10); if (LEN>MAX_LEN) MAX_LEN=LEN } END { print MAX_LEN }'`
    NONGAPGENOME=2867459933
    COV=`calc "$REALIGN_PAIRED*$READ_LENGTH/$NONGAPGENOME"`
    
    echo -e "$SAMPLE\t$REALIGN_PAIRED\t$ORIG_PAIRED\t$DIFF\t$READ_LENGTH\t$COV"
done > /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch13/flagstat_qc.txt

# ----------------------------------------
# 7. Calculate insert size for realigned BAM files
# ----------------------------------------

for SAMPLE in `cat $BATCH | cut -f 1`
do
	echo $SAMPLE
	bomb \
		-q long \
		-m 1 \
		-e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.ins.1.%J.log \
		-o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.ins.1.%J.log \
		"sambamba view -h -F 'paired and mate_is_reverse_strand and not (unmapped or mate_is_unmapped or reverse_strand or secondary_alignment or duplicate)' \
		/gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/$SAMPLE.bam | head -n 3000000 | awk '{OFS=\"\t\"; if (\$7==\"=\" && \$9>0 && \$9<1000) print \$9 }' \
		> /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/insert_dist_1.txt"
done


# ----------------------------------------
# 8. Generate histogram of insert size for each sample
# ----------------------------------------
for SAMPLE in `cat $BATCH | cut -f 1`
do
	echo $SAMPLE
	bomb \
        -q long \
		-m 1 \
		-e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.ins.1.hist.%J.log \
		-o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.ins.1.hist.%J.log \
		"Rscript /gscmnt/gc2802/halllab/ajscott/gtex_v7/gtex_insertDist_histogram.R \
        /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/insert_dist_1.txt $SAMPLE 1 /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/insert_dist.concordant.txt /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/insert_dist.discordant.txt"
done


# ----------------------------------------
# 9. Check that BAM header and EOF are present and okay
# ----------------------------------------

bomb -q hall-lab -e notes/samcheck.log -o notes/samcheck.log -m 1 "/gscmnt/gc2719/halllab/users/ajscott/bin/samtools quickcheck -v GTEX-*/*.bam > bad_bams.txt && echo 'all ok' || echo 'fail, see bad-bams.txt'"


# ----------------------------------------
# 10. Get coverage data for each sample
# ----------------------------------------
# note that earlier estimation method using (read_length*number_of_reads)/total_genome_length is a good approximation, so below is unnecessary

for SAMPLE in `cat $BATCH | cut -f 1`
do
    echo $SAMPLE
    bomb \
        -q long \
        -m 4 \
        -e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.coverage.%J.log \
        -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.coverage.%J.log \
        -g /ajscott/gtex_coverage \
        "bedtools genomecov -ibam /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/$SAMPLE.bam -bga \
        -g /gscmnt/gc2719/halllab/genomes/human/GRCh37/hs37_ebv/hs37_ebv.fasta.genome \
        | bgzip -c \
        > /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/$SAMPLE.coverage.bed.gz"
done

for SAMPLE in `cat $BATCH | cut -f 1`
do
    echo $SAMPLE
    bomb \
        -q long \
        -m 4 \
        -e /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.coverage.stats.%J.log \
        -o /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/log/$SAMPLE.coverage.stats.%J.log \
        "zcat /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/$SAMPLE.coverage.bed.gz \
        | cut -f 4 | zstats \
        > /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/$SAMPLE.coverage.stats.txt"
done

for SAMPLE in `cat $BATCH | cut -f 1`
do
    MEDIAN=`cat /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/$SAMPLE/qc/$SAMPLE.coverage.stats.txt | grep -m 1 "median:" | awk '{print $2}'`
    COV=`grep -m 1 $SAMPLE /gscmnt/gc2802/halllab/gtex_2016-06-01/realign_2016-06-06/notes/batch*/flagstat_qc.txt | awk '{print $6}'`
    echo -e "$SAMPLE\t$MEDIAN\t$COV"
done

######################################################################################################################################
# # ----------------------------------------
# # 7. Remove the original BAMs for completed samples
# # ----------------------------------------

# # Warning: use caution when deleting files
# FLAGSTAT_QC=/gscmnt/gc2802/halllab/gtex_realign_2015-03-16/notes/batch1/flagstat_qc.txt
# for SAMPLE in `cat $FLAGSTAT_QC | awk '{ if ($2==$3) print $1 }'`
# do
#     BAMPATH=`cat $BATCH | grep -m 1 $SAMPLE | cut -f 2`
#     echo $SAMPLE
    
#     rm $BAMPATH
# done





