#!/user/bin/env bash
### LAST VERSION UPDATE 27 MARCH 2023 (v2.6.6.1) - UPDATED AFTER IMPROVING PROJECT FOLDER STRUCTURE, ADDED FILEPATH ABBREVIATIONS TO MAKE CODE CLEANER, UPDATED ONCE MORE TO FIX FILEPATH AFTER LIFTING TO NEW LINUX
### THIS SCRIPT PERFORMS THE PRE-PROCESSING, MERGING, QC AND GWAS OF THE EIRA + SRQB WITH THE STR DATA TO CREATE THE PRIMARY DATA USED AS INPUT FOR THE SUBSEQUENT ANALYSIS ###

if ! [[ pwd == /home2/genetics/antobe/projects/OVERLAP_RA-CVD ]]; then cd /home2/genetics/antobe/projects/OVERLAP_RA-CVD; fi
if [[ -f DATA/OUTPUT/DATA/RA.bed ]] && [[ -f DATA/OUTPUT/DATA/RA.bim ]] && [[ -f DATA/OUTPUT/DATA/RA.fam ]]; then echo "FINISHED FILES ALREADY AVAILABILE AT 'DATA/OUTPUT/DATA/RA.p*'... EXITING BASH..."; exit 1; fi		#! Output message seems to suggest pfiles but looks for bfiles?

TMP=TEMPORARY/tmp-1    #SET FILEPATHS
EIRA_SRQB=/home2/genetics/antobe/data/EIRA-SRQB
STR=/home2/genetics/antobe/data/STR
SOFTWARE=../../software/
LOG=DATA/OUTPUT/LOG
RSQ=0.70    #SET THE IMPUTATION METRIC THRESHOLD

if [[ -f ${LOG}/PLog1.log ]]; then rm ${LOG}/PLog1.log; fi; touch ${LOG}/PLog1.log	#CREATE THE FILE THAT WILL HOST ALL THE PLINK OUTPUT
if [[ -f ${LOG}/EIRA-SRQB.log ]]; then rm ${LOG}/EIRA-SRQB.log; fi; touch ${LOG}/EIRA-SRQB.log
if [[ -f ${LOG}/SALTY.log ]]; then rm ${LOG}/SALTY.log; fi; touch ${LOG}/SALTY.log
if [[ -f ${LOG}/TwinGene.log ]]; then rm ${LOG}/TwinGene.log; fi; touch ${LOG}/TwinGene.log
if [[ -f ${LOG}/info.log ]]; then rm ${LOG}/info.log; fi; touch ${LOG}/info.log

### ### ### 1.1.1. EIRA-SPECIFIC FILTERING

if ! [[ -f ${TMP}/EIRA-SRQB_RSQ.bed ]] || ! [[ -f ${TMP}/EIRA-SRQB_RSQ.bim ]] || ! [[ -f ${TMP}/EIRA-SRQB_RSQ.fam ]]; then

    echo "PERFORMING INTRODUCTORY QUALITY CONTROL FOR THE EIRA-SRQB DATA... PLEASE HOLD..."
    tail -n +2 ${EIRA_SRQB}/All_chromosomes_info.txt | cut -f 1,7 | awk -v var="${RSQ}" '{ if ($2 >= var ) { print } }' > ${TMP}/RSQ_filt.txt    #ID SNPS WITH BELOW ACCEPTABLE IMPUTATION QUALITY
    ${SOFTWARE}/plink2 --bfile ${EIRA_SRQB}/eira-plus-others-imputed --extract ${TMP}/RSQ_filt.txt --make-bed --out ${TMP}/EIRA-SRQB_RSQ &>> ${LOG}/PLog1.log    #FILTER ON IMPUTATION QUALITY
    
    echo "RAW $(wc -l < ${EIRA_SRQB}/eira-plus-others-imputed.fam) $(wc -l < ${EIRA_SRQB}/eira-plus-others-imputed.bim)" >> ${LOG}/EIRA-SRQB.log    #COUNT THE REMAINING
    echo "IMP --||-- $(wc -l < ${TMP}/EIRA-SRQB_RSQ.bim)" >> ${LOG}/EIRA-SRQB.log
    echo "NAN --||-- --||--" >> ${LOG}/EIRA-SRQB.log    #ADD THIS TO MAKE IT COMPARABLE TO THE STR-PIPELINE...

fi

### ### ### 1.1.2. STR-SPECIFIC FILTERING

Rscript SCRIPTS/MISC/STRExtractRA.R    #RUN THIS TO ACQUIRE THE `STR_wRA.txt` FILE FOR FILTERING OF RA CASES
for FILENAME in SALTY TwinGene; do if ! [[ -f ${TMP}/${FILENAME}_RSQ.bed ]] || ! [[ -f ${TMP}/${FILENAME}_RSQ.bim ]] || ! [[ -f ${TMP}/${FILENAME}_RSQ.fam ]]; then    #SKIPPED THE INDIVIDUAL ONES IF FINISHED QC FILES ALREADY EXIST

    echo "PERFORMING INTRODUCTORY QUALITY CONTROL FOR THE ${FILENAME} DATA... PLEASE HOLD..."
    if [[ -f ${TMP}/${FILENAME}_pmerge_listfile.txt ]]; then rm ${TMP}/${FILENAME}_pmerge_listfile.txt; fi    #IN CASE IT ALREADY EXISTS IN SOME FORMAT DUE TO CRASHES - AVOID INFINITELY APPENDING
    N_SNP=0
    for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do    #FILTER ON IMPUTATION QUALITY STRAIGHT FROM THE VCF-FILES
	${SOFTWARE}/plink2 --vcf ${STR}/${FILENAME}_03_vcf/${CHR}.vcf.gz --exclude-if-info "INFO <= ${RSQ}" --maf 0.001 --make-pgen --out ${TMP}/${FILENAME}_fCHR${CHR} --threads 8 &>> ${LOG}/PLog1.log
	echo ${TMP}/${FILENAME}_fCHR${CHR} >> ${TMP}/${FILENAME}_pmerge_listfile.txt
	N_SNP_VCF=$(wc -l < ${TMP}/${FILENAME}_fCHR${CHR}.pvar); N_SNP=$((${N_SNP} + ${N_SNP_VCF}))    #COUNTING THE NUMBER OF SNPS IN THE RAW DATA - THERE IS AN ERROR OF ONE PER CHROMOSOME (BECAUSE THE HEADER IS COUNTED...)
	if [[ "${CHR}" =~ ^(4|8|12|16|20)$ ]]; then echo "COMPLETED WORK FOR CHR 1-${CHR}..."; fi    #GIVE OCCASSIONAL OUTPUT AS THIS IS INCREDIBLY SLOW
    done

    ${SOFTWARE}/plink2 --pmerge-list ${TMP}/${FILENAME}_pmerge_listfile.txt --merge-max-allele-ct 2 --make-bed --out ${TMP}/${FILENAME}_RSQ_wRA &>> ${LOG}/PLog1.log    #MERGE UP THE FILES
    ${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ_wRA --remove ${TMP}/STR_wRA.txt --make-bed --out ${TMP}/${FILENAME}_RSQ_NoSex &>> ${LOG}/PLog1.log    #FILTER OUT RA-CASES
    tail -n +2 ${STR}/${FILENAME}_04_pgen/${FILENAME}_HRC_imputation.psam | awk '{ print 0 "\t" $1 "_" $2 "\t" $3}' > ${TMP}/${FILENAME}_SEX.txt    #RE-OBTAIN SEX-INFORMATION LOST DURING CONVERSION
    ${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_RSQ_NoSex --update-sex ${TMP}/${FILENAME}_SEX.txt --make-bed --out ${TMP}/${FILENAME}_RSQ &>> ${LOG}/PLog1.log    #ADD SEX-INFORMATION BACK ONTO DATA

    echo "RAW $(wc -l < ${TMP}/${FILENAME}_RSQ_wRA.fam) ${N_SNP}" >> ${LOG}/${FILENAME}.log
    echo "IMP --||-- $(wc -l < ${TMP}/${FILENAME}_RSQ_wRA.bim)" >> ${LOG}/${FILENAME}.log
    echo "RA- $(wc -l < ${TMP}/${FILENAME}_RSQ_NoSex.fam) --||--" >> ${LOG}/${FILENAME}.log
    rm ${TMP}/${FILENAME}_fCHR*    #FREE UP SOME SPACE JUST IN CASE

fi; done

### ### ### 1.2. PRIMARY QUALITY CONTROL

for FILENAME in EIRA-SRQB SALTY TwinGene; do if ! [[ -f ${TMP}/${FILENAME}_QC.bed ]] || ! [[ -f ${TMP}/${FILENAME}_QC.bim ]] || ! [[ -f ${TMP}/${FILENAME}_QC.fam ]]; then    #SKIPPED THE INDIVIDUAL ONES IF FINISHED QC FILES ALREADY EXIST

    echo "PERFORMING PRIMARY QUALITY CONTROL FOR THE ${FILENAME} COHORT..."
    ${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_RSQ --check-sex 0.2 0.8 --out ${TMP}/${FILENAME}_RSQ_CHKSX &>> ${LOG}/PLog1.log
    grep PROBLEM ${TMP}/${FILENAME}_RSQ_CHKSX.sexcheck > ${TMP}/${FILENAME}_SXCHK_PROBLEM.txt
    ${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ --remove ${TMP}/${FILENAME}_SXCHK_PROBLEM.txt --make-bed --out ${TMP}/${FILENAME}_RSQ_SXCHK &>> ${LOG}/PLog1.log
    ${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ_SXCHK --mind 0.05 --chr 1-22 --make-bed --out ${TMP}/${FILENAME}_RSQ_SXCHK_MIND &>> ${LOG}/PLog1.log
    ${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ_SXCHK_MIND --geno 0.05 --make-bed --out ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO &>> ${LOG}/PLog1.log
    ${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO --maf 0.01 --make-bed --out ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO_MAF &>> ${LOG}/PLog1.log
    ${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO_MAF --hwe 1e-6 include-nonctrl --make-bed --out ${TMP}/${FILENAME}_QC &>> ${LOG}/PLog1.log
	
    echo "SEX $(wc -l < ${TMP}/${FILENAME}_RSQ_SXCHK.fam) --||--" >> ${LOG}/${FILENAME}.log
    echo "IND $(wc -l < ${TMP}/${FILENAME}_RSQ_SXCHK_MIND.fam) --||--" >> ${LOG}/${FILENAME}.log
    echo "GEN --||-- $(wc -l < ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO.bim)" >> ${LOG}/${FILENAME}.log
    echo "MAF --||-- $(wc -l < ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO_MAF.bim)" >> ${LOG}/${FILENAME}.log
    echo "HWE --||-- $(wc -l < ${TMP}/${FILENAME}_QC.bim)" >> ${LOG}/${FILENAME}.log
    rm ${TMP}/${FILENAME}_RSQ_*

fi; done

### ### ### 1.3.1. MERGING OF SALTY AND TWINGENE

echo "PERFORMING MERGING OF SALTY AND TWINGENE DATA..."

${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene_QC --set-missing-var-ids @:# --make-bed --out ${TMP}/TwinGene_idFix &>> ${LOG}/PLog1.log    #CONVERT SNP-IDS FOR SIMPLIFIED MERGING
${SOFTWARE}/plink2 --bfile ${TMP}/SALTY_QC --set-missing-var-ids @:# --make-bed --out ${TMP}/SALTY_idFix &>> ${LOG}/PLog1.log
${SOFTWARE}/plink --bfile ${TMP}/TwinGene_idFix --bmerge ${TMP}/SALTY_idFix --make-bed --out ${TMP}/TwinGene-SALTY_merged &>> ${LOG}/PLog1.log    #THIS WILL FAIL BUT OUTPUT A FILE FOR FILTERING 'POOR' SNPS

${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene_idFix --exclude ${TMP}/TwinGene-SALTY_merged-merge.missnp --make-bed --out ${TMP}/TwinGene_idFix_rmpoor &>> ${LOG}/PLog1.log    #CUT THOSE IDENTIFIED ABOVE
${SOFTWARE}/plink2 --bfile ${TMP}/SALTY_idFix --exclude ${TMP}/TwinGene-SALTY_merged-merge.missnp --make-bed --out ${TMP}/SALTY_idFix_rmpoor &>> ${LOG}/PLog1.log
${SOFTWARE}/plink --bfile ${TMP}/TwinGene_idFix_rmpoor --bmerge ${TMP}/SALTY_idFix_rmpoor --make-bed --out ${TMP}/TwinGene-SALTY_merged2 &>> ${LOG}/PLog1.log

cut -f 2 ${TMP}/TwinGene_idFix_rmpoor.bim > ${TMP}/TwinGene_SNPs.txt
cut -f 2 ${TMP}/SALTY_idFix_rmpoor.bim > ${TMP}/SALTY_SNPs.txt
grep -f ${TMP}/TwinGene_SNPs.txt ${TMP}/SALTY_SNPs.txt > ${TMP}/TwinGene-SALTY_SNPs.txt    #GRAB SNPS IN THE INTERSECTION OF TWINGENE AND SALTY
${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene-SALTY_merged2 --extract ${TMP}/TwinGene-SALTY_SNPs.txt --make-bed --out ${TMP}/TwinGene-SALTY_merged2_intersect &>> ${LOG}/PLog1.log    #EXCLUDE NON-OVERLAPPING SNPS

### ### ### 1.3.2. MERGING OF SALTY AND TWINGENE WITH EIRA AND SRQB

echo "PERFORMING MERGING OF SALTY AND TWINGENE DATA WITH EIRA AND SRQB DATA..."

${SOFTWARE}/plink2 --bfile ${TMP}/EIRA-SRQB_QC --set-all-var-ids @:# --make-bed --out ${TMP}/EIRA-SRQB_idFix &>> ${LOG}/PLog1.log    #CHANGE IDS TO REMOVE DEPENDENCE ON ALLELES
${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene-SALTY_merged2_intersect --set-all-var-ids @:# --make-bed --out ${TMP}/TwinGene-SALTY_merged2_intersect_idFix &>> ${LOG}/PLog1.log    #REMOVE RS-NUMBERS TO MAKE ALL COMPARABLE
${SOFTWARE}/plink --bfile ${TMP}/TwinGene-SALTY_merged2_intersect_idFix --bmerge ${TMP}/EIRA-SRQB_idFix --make-bed --out ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged &>> ${LOG}/PLog1.log

${SOFTWARE}/plink --bfile ${TMP}/TwinGene-SALTY_merged2_intersect_idFix --exclude ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged-merge.missnp --make-bed --out ${TMP}/TwinGene-SALTY_merged2_intersect_idFix_rmpoor &>> ${LOG}/PLog1.log
${SOFTWARE}/plink --bfile ${TMP}/EIRA-SRQB_idFix --exclude ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged-merge.missnp --make-bed --out ${TMP}/EIRA-SRQB_idFix_rmpoor &>> ${LOG}/PLog1.log
${SOFTWARE}/plink --bfile ${TMP}/TwinGene-SALTY_merged2_intersect_idFix_rmpoor --bmerge ${TMP}/EIRA-SRQB_idFix_rmpoor --make-bed --out ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged2 &>> ${LOG}/PLog1.log

cut -f 2 ${TMP}/TwinGene-SALTY_merged2_intersect_idFix_rmpoor.bim > ${TMP}/TwinGene-SALTY_SNPs2.txt
cut -f 2 ${TMP}/EIRA-SRQB_idFix_rmpoor.bim > ${TMP}/EIRA-SRQB_SNPs.txt
grep -f ${TMP}/TwinGene-SALTY_SNPs2.txt ${TMP}/EIRA-SRQB_SNPs.txt > ${TMP}/TwinGene-SALTY-EIRA-SRQB_SNPs.txt
${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged2 --extract ${TMP}/TwinGene-SALTY-EIRA-SRQB_SNPs.txt --make-bed --out ${TMP}/RA_merged &>> ${LOG}/PLog1.log

echo "MRG $(wc -l < ${TMP}/RA_merged.fam) $(wc -l < ${TMP}/RA_merged.bim)" >> ${LOG}/info.log
rm ${TMP}/EIRA-SRQB_*; rm ${TMP}/SALTY_*; rm ${TMP}/TwinGene*    #FOR STORAGE PURPOSES

### ### ### 1.4. POST-MERGING FINALIZING QC

echo "PERFORMING QUALITY CONTROL OF THE MERGED DATA AND CLEANING UP FOLDERS..."

${SOFTWARE}/plink2 --bfile ${TMP}/RA_merged --geno 0.05 --make-bed --out ${TMP}/RA_merged_GENO &>> ${LOG}/PLog1.log
${SOFTWARE}/plink2 --bfile ${TMP}/RA_merged_GENO --mind 0.05 --make-bed --out ${TMP}/RA_merged_GENO_MIND &>> ${LOG}/PLog1.log
${SOFTWARE}/plink2 --bfile ${TMP}/RA_merged_GENO_MIND --maf 0.01 --make-bed --out ${TMP}/RA_merged_GENO_MIND_MAF &>> ${LOG}/PLog1.log
${SOFTWARE}/plink --bfile ${TMP}/RA_merged_GENO_MIND_MAF --hwe 1e-6 include-nonctrl --make-bed --out DATA/OUTPUT/DATA/RA &>> ${LOG}/PLog1.log

echo "IND --||-- $(wc -l < ${TMP}/RA_merged_GENO.bim)" >> ${LOG}/info.log
echo "GEN $(wc -l < ${TMP}/RA_merged_GENO_MIND.fam) --||--" >> ${LOG}/info.log
echo "MAF --||-- $(wc -l < ${TMP}/RA_merged_GENO_MIND_MAF.bim)" >> ${LOG}/info.log
echo "HWE --||-- $(wc -l < DATA/OUTPUT/DATA/RA.bim)" >> ${LOG}/info.log

### ### ### 1.5. FILE CLEANUP AND EXITING BASH

cat ${LOG}/info.log >> ${LOG}/EIRA-SRQB.log; cat ${LOG}/info.log >> ${LOG}/SALTY.log; cat ${LOG}/info.log >> ${LOG}/TwinGene.log; rm ${LOG}/info.log
rm ${TMP}/*
echo "DONE! PLEASE FIND THE OUTPUT WITHIN THE 'DATA/OUPUT/DATA/' FOLDER..."

### TO DO:
# 1.1. Current merging approach is a bit rough. One could possibly do it all in one go. The fixing of SNP IDs could possibly be handled with less repetition. It works for now and is a bit complex
#      so I'll currently leave it as is.
### NOTES:
# 1.1. Version 1.0 created on 22 APRIL 2022.
#      Version 1.0 performs the pre-processing (currently, none) and QC of the EIRA + SRQb data. It has room for minor improvements but otherwise works as intended.
#      Subsequent versions will need improvements to allow for the same procedures being carried out in the STR data.
#      In general, it runs like the previous project script. Updates include better handling of PLINK logs (their output is now suppressed and re-directed to a single file containing all of them in the appropriate order).
#      I've also added continuous updates on the number of remaining individuals (N) and SNPs (M). These are continously updated into a file which should make it easier to follow.
# 1.2. Version 2.0. created on 16 MAY 2022.
#      Now contains pre-processing and QC of the SALTY and TwinGene data as well. Still, there's no combination of files into one distinct file-set, which will be added in the next version.
# 1.3. The current approach for merging comes with minor issues that are worth highlighting for transparency (see the associated .docx text for a more thorough description (currently resides in Main/Projekt/...)). 
#      Current approach simply cuts all problem SNPs without making any attempt at resolving conflicting SNP IDs. PLINK handles duplications poorly and so this is a simple way to make the code finish, losing only a handful of SNPs.
#      While this comes with the cost of (unnecessary) loss of data, a bigger issue is that PLINK 'combines' SNPs with duplicated IDs where the (unordered) set of alleles match.
#      For instance, a SNP rs123 A G will be combined with the SNP rs123 G A, while the SNPs rs456 C T and rs456 C G will be cut. What 'combined' means is somewhat elusive, making it somewhat troubling.
#      Nevertheless, such SNPs should be few and we'll leave it for later versions to fix it.
# 1.3. Version 2.6. created on 13 DEC 2022.
#      This version shortened the script somewhat, to make it more readable. It also added X-chromosome data to the STR-cohorts, meaning filtering on sex was now possible. Additionally, the script now actively excludes
#      participants of the STR-cohorts that were considered positive for RA (at the time of writing, this corresponds to having one ICD-10 code for an RA-related disease).
