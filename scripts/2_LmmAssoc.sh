#!/user/bin/env bash
### LAST VERSION UPDATE - 8 MAY 2023 (v.1.4.1) - UPDATED SOME SCRIPTS, FIXED SOME FILEPATHS, AFTER MOVING TO THE NEW CLUSTER
### THIS SCRIPT RUNS THE LINEAR MIXED MODEL AND SUBSEQUENT ASSOCIATION TESTING OF THE RA CASE-CONTROL SAMPLE PRODUCED FROM `1_preprocessing.sh`.

if ! [[ pwd == /home2/genetics/antobe/projects/OVERLAP_RA-CVD ]]; then cd /home2/genetics/antobe/projects/OVERLAP_RA-CVD; fi
if ! [[ -f DATA/OUTPUT/DATA/RA.bed ]] || ! [[ -f DATA/OUTPUT/DATA/RA.bim ]] || ! [[ -f DATA/OUTPUT/DATA/RA.fam ]]; then echo 'FAILED TO FIND INPUT FILES... PLEASE RUN `SCRIPTS/1_*.sh` PRIOR TO RUNNING THIS SCRIPT... EXITING BASH...'; exit 1; fi

TMP=TEMPORARY/tmp-2	#SET FILEPATHS
LOG=DATA/OUTPUT/LOG
DATA=DATA/OUTPUT/DATA
SOFTWARE=../../software/
CACHE=TEMPORARY/CACHE
RAW=/home2/genetics/antobe/data

if [[ -f ${LOG}/PLog2.log ]]; then rm ${LOG}/PLog2.log; fi; touch ${LOG}/PLog2.log	#CREATE FILES THAT WILL HOST OUTPUT

### ### ### 2.1. PRE-PROCESSING OF DATA

echo "CONSTRUCTING PHENOTYPE FILES..."
awk -v OFS='\t' '{ print $1, $2, $6}' ${DATA}/RA.fam > ${TMP}/PHENO_wNA_wUnAdj.txt
awk ' { $3 = ($3 == "-9" ? 1 : $3) } 1' ${TMP}/PHENO_wNA_wUnAdj.txt > ${TMP}/PHENO_wUnAdj.txt    #NOTE THAT THERE IS ONE EIRA/SRQB INDIVIDUAL WITH NA AS RA-STATUS - THIS BECOMES A CONTROL TOO - THIS CHANGES ALL WITH NA STATUS TO CONTROL
awk ' {$1 ; $2 ; $3 -= 1} 1 ' OFS='\t' ${TMP}/PHENO_wUnAdj.txt > ${TMP}/PHENO.txt    #SUBTRACT ONE FROM THE PHENOTYPE COLUMN TO APPEASE FORMAT WITHIN GCTA

echo 'COMPUTING PRINCIPAL COMPONENTS, PLEASE HOLD...'
${SOFTWARE}/plink2 --bfile ${DATA}/RA --indep-pairwise 1000 50 0.05 --out ${TMP}/RA_ldprune &>> ${LOG}/PLog2.log    #PARAMETERS IDENTICAL TO THOSE USED WTIHIN SIMULATIONS OF `Jiang, Zheng, Fang and Yang (2021)`
${SOFTWARE}/plink2 --bfile ${DATA}/RA --extract ${TMP}/RA_ldprune.prune.in --out ${TMP}/RA_PRUNED --make-pgen &>> ${LOG}/PLog2.log
${SOFTWARE}/plink2 --pfile ${TMP}/RA_PRUNED --pca 10 approx --out ${TMP}/PCA --threads 8 &>> ${LOG}/PLog2.log    #USING THE `Galinsky et al (2016)` (`fastPCA`) APPROXIMATION, SEE NOTES

Rscript SCRIPTS/MISC/PCScree.R ${TMP}/PCA.eigenval    
read ncomp    #DELETED THE CHECK FOR THE APPROPRIATE INPUT FORMAT - JUST LISTEN TO WHAT THE ECHO SAYS, DON'T BE AN IDIOT...
cov_fields=$((ncomp + 2))
tail -n +2 ${TMP}/PCA.eigenvec | cut -f 1-${cov_fields} > ${TMP}/PCA.txt

Rscript SCRIPTS/MISC/ExtractQCov.R 0    #EXTRACTS AGE AND LINKS IT ONTO THE ABOVE PCS
awk -v OFS='\t' '{ print $1, $2, $5}' ${DATA}/RA.fam > ${TMP}/SEX.txt    #EXTRACT SEX

### ### ### 2.2. fastGWA-GLMM WITH THE GCTA SOFTWARE ### ### ###

if ! [[ -f ${CACHE}/RA_GRM.grm.sp ]]; then
	echo "CONSTRUCTING GRM... PLEASE HOLD..."
	
	awk '{sub(/\:\w+\:\w+$/, "", $1); print $1}' ${RAW}/EIRA-SRQB/SNP-ID_linkage.txt > ${TMP}/OLD_ID.txt	#GET SNP ID WITHOUT ALLELES (ONLY CHR:BP)
	cut -f 2 ${RAW}/EIRA-SRQB/SNP-ID_linkage.txt > ${TMP}/NEW_ID.txt
	paste ${TMP}/OLD_ID.txt ${TMP}/NEW_ID.txt > ${TMP}/MOD_SNP-ID_linkage.txt
	${SOFTWARE}/plink2 --bfile ${DATA}/RA --update-name ${TMP}/MOD_SNP-ID_linkage.txt --out ${TMP}/RA_rsID --make-pgen &>> ${LOG}/PLog2.log    #UPDATE SNP IDS TO RS-IDS FOR LINKAGE TO HAPMAP3 DATA
	cut -f 2 ${RAW}/HAPMAP3/hapmap3_r3_b36_fwd.consensus.qc.poly.map > ${TMP}/hm3_snps.txt

	${SOFTWARE}/plink2 --pfile ${TMP}/RA_rsID --extract ${TMP}/hm3_snps.txt --out ${TMP}/RA_rsID_hm3 --make-pgen &>> ${LOG}/PLog2.log    #EXTRACT HM3 SNPS FROM OUR DATA
	${SOFTWARE}/plink2 --pfile ${TMP}/RA_rsID_hm3 --indep-pairwise 1000 100 0.90 --out ${TMP}/RA_rsID_hm3_ldprune &>> ${LOG}/PLog2.log
	${SOFTWARE}/plink2 --pfile ${TMP}/RA_rsID_hm3 --extract ${TMP}/RA_rsID_hm3_ldprune.prune.in --out ${TMP}/RA_rsID_hm3_PRUNED --make-bed &>> ${LOG}/PLog2.log    #PRE-FILTERING AD DONE FOR UKB IN `Jiang, Zheng, Fang and Yang (2021)`
	${SOFTWARE}/gcta-1.94.1 --bfile ${TMP}/RA_rsID_hm3_PRUNED --make-grm --sparse-cutoff 0.05 --thread-num 8 --out ${CACHE}/RA_GRM    #USE THE CUTOFF OF 0.05 AS RECOMMENDED IN `Jiang, Zheng and Qi et al (2019)`

	cat ${CACHE}/RA_GRM.grm.sp >> ${LOG}/PLog2.log    #SINCE WE DO NOT SILENCE THE LOG (IT SHOWS INCELY THE PROGRESS OF THE GRM) WE ALSO APPEND IT TO THE LOG-FILE FOR CLARITY
fi    

for PHENO in RA SPOS_RA SNEG_RA; do
	echo "RUNNING THE ASSOCIATION ANALYSIS FOR THE ${PHENO} PHENOTYPE, PLEASE HOLD... THIS MAY TAKE A WHILE..."
	
	Rscript SCRIPTS/MISC/ExtractSubPheno.R ${PHENO}
	${SOFTWARE}/plink2 --bfile ${DATA}/RA --keep ${TMP}/${PHENO}_TargetInds.txt --make-bed --out ${TMP}/${PHENO} &>> ${LOG}/PLog2.log
	FCOUNT="GWA $(wc -l < ${TMP}/${PHENO}.fam) $(wc -l < ${TMP}/${PHENO}.bim)"; echo ${FCOUNT} >> ${LOG}/EIRA-SRQB.log; echo ${FCOUNT} >> ${LOG}/TwinGene.log; echo ${FCOUNT} >> ${LOG}/SALTY.log
	${SOFTWARE}/gcta-1.94.1 --bfile ${TMP}/${PHENO} --grm-sparse ${CACHE}/RA_GRM --fastGWA-mlm-binary --geno 0.99 --pheno ${TMP}/PHENO.txt --qcovar ${TMP}/QCOV.txt --covar ${TMP}/SEX.txt --thread-num 8 --out ${DATA}/${PHENO}_GWAS &>> ${LOG}/PLog2.log

	Rscript SCRIPTS/MISC/Manhattan.R ${DATA}/${PHENO}_GWAS.fastGWA
done

rm ${TMP}/*
echo "DONE! PLEASE FIND ALL THE OUTPUT FILES AT 'DATA/OUTPUT/DATA/'"

### TO DO:
# 1.2. It might be possible to prepare the SNP ID linkage file in one go, instead of across three lines, but I have yet to check this more closely. It seems that all can be done either with awk or at least with pipe.
### NOTE:
# 2.1. Currently, I use a procedure similar to what was done within the fastGWA-GLMM paper simulations, when constructing the PCs. 
#      If the current procedure is deemed too slow, we can possibly expediate the computations by going into some other software (than PLINK) such as flashPCA2 which seems faster (also used in fastGWA-GLMM paper).
# 2.2. I also wish to raise a concern regarding the actual computations. I believe that the PCA should be performed within independent observations (hence the LD pruning) and I also believe that high relatedness
#      could bias the computations - something that may be an issue with the extensive twin-relationships in my STR data. As such I reckon we should make sure that including PCs as a covariate does not modify
#      our results in a major way, and we should be careful regarding their role within this analysis.
# 2.3. On a similar note, I use the `approx` option available within PLINK2 that uses the Galinsky et al (2016) 'fastPCA' algorithm, a MUCH faster option than the base algorithm. However, I'm not entirely
#      convinced on the validity of this one, and I might still want to check, just like in 2.2., how well it works for our purposes.
# 2.4. At the request of Helga, I tried running the code with an additional covariate, COHORT, which identifies whether an individual is an EIRA-SRQB, a TWINGENE or a SALTY participant. This code remains
#      but has been commented out to not be used when running the script. This model led to non-invertible matrices in all scenarios tried by me. Such results could be indicative of greater issues, but I
#      choose to believe that it is an issue with some type of perfect separation. I think adjusting for COHORT does more harm than good, even if it would be possible to obtain complete estiamtes with it as a covariate.
#      EDIT: I've actually removed this code now. It wasn't used and remained solely as a relic of a previously failed attempt. This way the script is tighter.
#      EDIT2: I believe that the pre-processing has changed since I last tested the above described chunk of code. It is possible that the same issue no longer occurs, but I have not seen any reason to test it yet.
