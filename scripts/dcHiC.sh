#!/bin/bash


#------------------------------------------------------------------------------------------------------------------#
# DCHIC PIPELINE FOR CALLING COMPARTMENT REPOSITIONING EVENTS                                                      #
#------------------------------------------------------------------------------------------------------------------#
# In this pipeline we use dcHiC to call regions of differential compartmentalization                               #
# in three tumors, comparing them with a control sample                                                            #
# The output of the dcHiC pipeline is converted to the CoRE format for direct comparison                           #
# with the diffComp software output.                                                                               #
#                                                                                                                  #
# The pipeline is composed by these steps:                                                                         #
# 1. We pre-process the .hic files of each sample (we have two replicates for each sample).                        #
#    This requires to convert the .hic files into .matrix files together with a bed file explaning                 #
#    the location of Hi-C bins. Bins intersecting a blacklist of bad bins are removed from the analysis            #
# 2. For each sample, we use dcHi-C to determine PCs and select the best ones                                      #
# 3. For each sample comparison, we use dcHi-C to get differential bins affected by compartment repositioning      #
# 4. For each sample comparison, we aggregate adjacent changing bins to form a "pseudo-CoRE". The output           #
#    has the same format as the CoREs returned by diffComp.                                                        #
#------------------------------------------------------------------------------------------------------------------#



BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"
DCHIC_PREPROCESS_CMD_PATH="/mnt/etemp/luca/software/dcHiC/utility/preprocess.py"
DCHIC_CMD_PATH="/mnt/etemp/luca/software/dcHiC/dchicf.r"
DCHIC_RESULTS_PATH="data/dchic"
DCHIC_RESOLUTION=50000
DCHIC_THREADS=10
DCHIC_INPUT_PATH="${DCHIC_RESULTS_PATH}/input.txt"
CORES_MAX_DISTANCE_MERGING=500000

HIC_PATH="/mnt/ndata/Juan/HiC-maps"
SAMPLES=(
    "RPE_C1_30:RPE_TP53_Ctrl"
    "RPE_C2_30:RPE_TP53_Ctrl"
    "TP53_20wTumo11_30:RPE_TP53_20w0T1"
    "TP53_20wTumo12_30:RPE_TP53_20w0T1"
    "TP53_20wTumo21_30:RPE_TP53_20w0T2"
    "TP53_20wTumo22_30:RPE_TP53_20w0T2"
    "TP53_20wTumo31_30:RPE_TP53_20w0T3"
    "TP53_20wTumo32_30:RPE_TP53_20w0T3"
)

COMPARISONS=(
    "RPE_TP53_Ctrl:RPE_TP53_20w0T1"
    "RPE_TP53_Ctrl:RPE_TP53_20w0T2"
    "RPE_TP53_Ctrl:RPE_TP53_20w0T3"
)

mkdir -p ${DCHIC_RESULTS_PATH}
> ${DCHIC_INPUT_PATH}


echo "---------------------------------------"
echo " PRE-PROCESSING HI-C SAMPLES FOR dcHiC"
echo "---------------------------------------"

echo "Downloading blacklisted regions"

if [[ ! -f ${DCHIC_RESULTS_PATH}/hg19-blacklist.v2.bed ]]; then
    wget ${BLACKLIST_URL} -O ${DCHIC_RESULTS_PATH}/hg19-blacklist.v2.bed.gz --quiet
    gzip -d ${DCHIC_RESULTS_PATH}/hg19-blacklist.v2.bed.gz
fi



for sample in "${SAMPLES[@]}"
do
    replicate=$(echo ${sample} | cut -d":" -f 1)
    experiment=$(echo ${sample} | cut -d":" -f 2)
    echo "- ${replicate}"
    echo "---- Experiment: ${experiment}"
    replicate_hic_path="${HIC_PATH}/${replicate}.hic"
    echo "---- Replicate path: ${replicate_hic_path}"

    replicate_mat_path=${DCHIC_RESULTS_PATH}/${replicate}_${DCHIC_RESOLUTION}.matrix
    replicate_bed_path=${DCHIC_RESULTS_PATH}/${replicate}_${DCHIC_RESOLUTION}_abs.bed


    echo "---- Preprocessing"
    if [[ ! -f ${replicate_mat_path} ]] || [[ ! -f ${replicate_bed_path} ]]; then
        python ${DCHIC_PREPROCESS_CMD_PATH} -input hic \
                                            -file ${replicate_hic_path} \
                                            -res ${DCHIC_RESOLUTION} \
                                            -prefix ${DCHIC_RESULTS_PATH}/${replicate} \
                                            -removeChr y,mt,m
    fi
    replicate_bed_with_black_path=${DCHIC_RESULTS_PATH}/${replicate}_${DCHIC_RESOLUTION}_abs_blacklisted.bed
    if [[ ! -f ${replicate_bed_with_black_path} ]]; then
        bedtools intersect -c \
                           -a ${replicate_bed_path} \
                           -b ${DCHIC_RESULTS_PATH}/hg19-blacklist.v2.bed \
                | awk -v FS='\t' -v OFS='\t' '{$5 = ($5 > 1 ? 1 : $5); print $1,$2,$3,$4,$5}' \
                > ${replicate_bed_with_black_path}
    fi

    echo -e "$(basename ${replicate_mat_path})\t$(basename ${replicate_bed_with_black_path})\t${replicate}\t${experiment}" >> ${DCHIC_INPUT_PATH}
done


echo "---------------"
echo " RUNNING DCHIC "
echo "---------------"

mkdir -p ${DCHIC_RESULTS_PATH}/cores

# cwd=$(pwd)
# cd ${DCHIC_RESULTS_PATH}
# Rscript ${DCHIC_CMD_PATH} --file ${cwd}/${DCHIC_INPUT_PATH} --pcatype cis --dirovwt T
# Rscript ${DCHIC_CMD_PATH} --file ${cwd}/${DCHIC_INPUT_PATH} --pcatype select --dirovwt T --genome hg19
# cd ${cwd}


for comparison in "${COMPARISONS[@]}"
do
    exp1=$(echo ${comparison} | cut -d":" -f 1)
    exp2=$(echo ${comparison} | cut -d":" -f 2)
    comparison_input="${DCHIC_RESULTS_PATH}/${exp1}_vs_${exp2}.txt"
    > ${comparison_input}
    cat ${DCHIC_INPUT_PATH} | grep ${exp1} >> ${comparison_input}
    cat ${DCHIC_INPUT_PATH} | grep ${exp2} >> ${comparison_input}

    # cd ${DCHIC_RESULTS_PATH}
    # Rscript ${DCHIC_CMD_PATH} --file $(basename ${comparison_input}) --pcatype analyze --dirovwt T --diffdir ${exp1}_vs_${exp2}
    # Rscript ${DCHIC_CMD_PATH} --file $(basename ${comparison_input}) --pcatype subcomp --dirovwt T --diffdir ${exp1}_vs_${exp2}
    # Rscript ${DCHIC_CMD_PATH} --file $(basename ${comparison_input}) --pcatype viz --diffdir ${exp1}_vs_${exp2} --genome hg19
    # cd ${cwd}


    comparison_fdr_path="${DCHIC_RESULTS_PATH}/DifferentialResult/${exp1}_vs_${exp2}/fdr_result/differential.intra_sample_combined.pcQnm.bedGraph"
    comparison_fdr_filtered_path="${DCHIC_RESULTS_PATH}/DifferentialResult/${exp1}_vs_${exp2}/fdr_result/differential.intra_sample_combined.Filtered.pcQnm.bedGraph"

    cat ${comparison_fdr_filtered_path} \
        | awk -v FS='\t' -v OFS='\t' 'NR>1{print $1,$2,$3,$(NF - 6),$(NF - 5),$(NF - 1)}' \
        | bedtools merge -d ${CORES_MAX_DISTANCE_MERGING} -c 4,5,6 -o mean \
        | awk -v FS='\t' -v OFS='\t' 'BEGIN{print "chr","start","end","value","pvalue"}{print $1,$2,$3,$5-$4,$6}' \
        > ${DCHIC_RESULTS_PATH}/cores/${exp1}_vs_${exp2}.tsv

    cat ${DCHIC_RESULTS_PATH}/cores/${exp1}_vs_${exp2}.tsv \
        | awk -v FS='\t' -v OFS='\t' 'NR>1{print $1,$2,$3,$1":"$2"-"$3,$4,".",$2,$3,($4 < 0 ? "0,0,255" : "255,0,0")}' \
        > ${DCHIC_RESULTS_PATH}/cores/${exp1}_vs_${exp2}.bed

done

# cd ${DCHIC_RESULTS_PATH}
# Rscript ${DCHIC_CMD_PATH} --file ${cwd}/${DCHIC_INPUT_PATH} --pcatype analyze --dirovwt T --diffdir ALL
# Rscript ${DCHIC_CMD_PATH} --file ${cwd}/${DCHIC_INPUT_PATH} --pcatype subcomp --dirovwt T --diffdir ALL
# Rscript ${DCHIC_CMD_PATH} --file ${cwd}/${DCHIC_INPUT_PATH} --pcatype viz --diffdir ALL --genome hg19
# cd ${cwd}
