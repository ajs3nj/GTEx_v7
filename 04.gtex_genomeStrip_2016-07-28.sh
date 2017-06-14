# 2016-07-28

# run genome strip on new GTEx dataset
# Genome STRiP version 2.00.1636

pwd
# /gscmnt/gc2802/halllab/gtex_2016-06-01/genomeStrip_2016-07-28

bomb -q long -m 8 -e gtex_gs.%J.log -o gtex_gs.%J.log "lsmake"