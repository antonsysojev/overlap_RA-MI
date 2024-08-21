#!/user/bin/env bash
### LAST VERSION UPDATE - 8 MAY 2023 (v1.1) - UPDATED AFTER MOVE TO NEW CLUSTER
### THIS SCRIPT EXTRACTS THE REFERENCE PANEL DATA FOR THE OVERLAP PROJECT FROM CLEANED GENOTYPE DATA

TMP=TEMPORARY/tmp-3
RAW=/home2/genetics/antobe/data
SOFTWARE=../../software/
touch ${TMP}/RPlog.log

### ### ### 3.1.1. COHORT SPECIFIC FILTERINGS

echo "BEGINNING COHORT-SPECIFIC QUALITY CONTROL..."

if ! [[ -f ${TMP}/EIRA-SRQB_RSQ.bed ]] || ! [[ -f ${TMP}/SALTY_RSQ.bed ]] || ! [[ -f ${TMP}/TwinGene_RSQ.fam ]]; then    #SKIPS THE CHUNK IF FILES ALREADY AVAILABLE

	awk -F ' ' '{ if ($6 != 2) print $1, $2}' ${RAW}/EIRA-SRQB/eira-plus-others-imputed.fam > ${TMP}/controls.txt    #FILTER EIRA-SRQB ON RA-CASES AND IMPUTATION QUALITY
	tail -n +2 ${RAW}/EIRA-SRQB/All_chromosomes_info.txt | cut -f 1,7 | awk '{ if ($2 >= 0.70) { print } }' > ${TMP}/RSQ_filt.txt
	${SOFTWARE}/plink2 --bfile ${RAW}/EIRA-SRQB/eira-plus-others-imputed --keep ${TMP}/controls.txt --extract ${TMP}/RSQ_filt.txt --make-bed --out ${TMP}/EIRA-SRQB_RSQ &>> ${TMP}/RPlog.log
	echo "COHORT-SPECIFIC QUALITY CONTROL COMPLETED FOR EIRA-SRQB..."

	Rscript SCRIPTS/MISC/STRExtractRA.R; mv TEMPORARY/tmp-1/STR_wRA.txt ${TMP}/STR_wRA.txt    #SIMPLE SOLUTION TO THE SCRIPT OUTPUTING THE FILE IN tmp-1
	for FILENAME in SALTY TwinGene; do    #FILTER SALTY AND TwinGene ON RA-CASES AND IMPUTATION QUALITY

		for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do
			${SOFTWARE}/plink2 --vcf ${RAW}/STR/${FILENAME}_03_vcf/${CHR}.vcf.gz --exclude-if-info "INFO <= 0.70" --maf 0.001 --make-pgen --out ${TMP}/${FILENAME}_fCHR${CHR} --threads 8 &>> ${TMP}/RPlog.log
			echo ${TMP}/${FILENAME}_fCHR${CHR} >> ${TMP}/${FILENAME}_pmerge_listfile.txt
		done

		${SOFTWARE}/plink2 --pmerge-list ${TMP}/${FILENAME}_pmerge_listfile.txt --merge-max-allele-ct 2 --make-bed --out ${TMP}/${FILENAME}_RSQ_wRA &>> ${TMP}/RPlog.log
		${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ_wRA --remove ${TMP}/STR_wRA.txt --make-bed --out ${TMP}/${FILENAME}_RSQ_NoSex &>> ${TMP}/RPlog.log
		tail -n +2 ${RAW}/STR/${FILENAME}_04_pgen/${FILENAME}_HRC_imputation.psam | awk '{ print 0 "\t" $1 "_" $2 "\t" $3}' > ${TMP}/${FILENAME}_SEX.txt
		${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_RSQ_NoSex --update-sex ${TMP}/${FILENAME}_SEX.txt --make-bed --out ${TMP}/${FILENAME}_RSQ &>> ${TMP}/RPlog.log
		rm ${TMP}/${FILENAME}_fCHR*
		echo "COHORT-SPECIFIC QUALITY CONTROL COMPLETED FOR ${FILENAME}..."

	done
fi

### ### ### 3.1.2. PRIMARY QUALITY CONTROL

echo "PERFORMING PRIMARY QUALITY CONTROL..."

for FILENAME in EIRA-SRQB SALTY TwinGene; do if ! [[ -f ${TMP}/${FILENAME}_QC.bed ]] || ! [[ -f ${TMP}/${FILENAME}_QC.bed ]] ||! [[ -f ${TMP}/${FILENAME}_QC.bed ]]; then

	${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_RSQ --check-sex 0.2 0.8 --out ${TMP}/${FILENAME}_RSQ_CHKSX &>> ${TMP}/RPlog.log
	grep PROBLEM ${TMP}/${FILENAME}_RSQ_CHKSX.sexcheck > ${TMP}/${FILENAME}_SXCHK_PROBLEM.txt
	${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ --remove ${TMP}/${FILENAME}_SXCHK_PROBLEM.txt --make-bed --out ${TMP}/${FILENAME}_RSQ_SXCHK &>> ${TMP}/RPlog.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_RSQ_SXCHK --mind 0.05 --geno 0.05 --maf 0.01 --chr 1-22 --make-bed --out ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO_MAF &>> ${TMP}/RPlog.log
	${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_RSQ_SXCHK_MIND_GENO_MAF --hwe 1e-6 include-nonctrl --make-bed --out ${TMP}/${FILENAME}_pQC &>> ${TMP}/RPlog.log

	${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_pQC --indep-pairwise 100 5 0.2 --out ${TMP}/${FILENAME}_pQC_LDPRUNE &>> ${TMP}/RPlog.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_pQC --extract ${TMP}/${FILENAME}_pQC_LDPRUNE.prune.in --make-bed --out ${TMP}/${FILENAME}_pQC_LDPRUNE &>> ${TMP}/RPlog.log
	${SOFTWARE}/plink --bfile ${TMP}/${FILENAME}_pQC_LDPRUNE --genome --min 0.10 --out ${TMP}/${FILENAME}_pQC_LDPRUNE_REL &>> ${TMP}/RPlog.log
	${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_pQC_LDPRUNE --remove ${TMP}/${FILENAME}_pQC_LDPRUNE_REL.genome --make-bed --out ${TMP}/${FILENAME}_pQC_LDPRUNE_REL &>> ${TMP}/RPlog.log

	PCA_INPUT=${FILENAME}_pQC_LDPRUNE_REL
	for ITER in $(seq 1 5); do
		if [[ ${DATA} =~ ^TwinGene$ ]]; then ${SOFTWARE}/plink2 --bfile ${TMP}/${PCA_INPUT} --pca 10 approx --out pca --threads 8 &>> ${TMP}/RPlog.log
		else ${SOFTWARE}/plink2 --bfile ${TMP}/${PCA_INPUT} --pca 10 --out pca --threads 8 &>> ${TMP}/RPlog.log; fi

		Rscript SCRIPTS/MISC/pcafiltration.R
		${SOFTWARE}/plink2 --bfile ${TMP}/${PCA_INPUT} --remove pca_outliers.txt --out ${TMP}/${PCA_INPUT}_PCA${ITER} --make-bed &>> ${TMP}/RPlog.log
		rm pca*; PCA_INPUT=${PCA_INPUT}_PCA${ITER}
	done

	${SOFTWARE}/plink2 --bfile ${TMP}/${FILENAME}_pQC --keep ${TMP}/${PCA_INPUT}.fam --make-bed --out ${TMP}/${FILENAME}_QC &>> ${TMP}/RPlog.log
	rm ${TMP}/${FILENAME}_RSQ*; rm ${TMP}/${FILENAME}_pQC*
	echo "PRE-MERGING COHORT-SPECIFIC QUALITY CONTROL COMPLETED FOR ${FILENAME}..."

fi; done

echo "PERFORMING MERGING OF TWINGENE WITH DATA"

${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene_QC --set-missing-var-ids @:# --make-bed --out ${TMP}/TwinGene_idFix &>> ${TMP}/RPlog.log    #CONVERT SNP-IDS FOR SIMPLIFIED MERGING
${SOFTWARE}/plink2 --bfile ${TMP}/SALTY_QC --set-missing-var-ids @:# --make-bed --out ${TMP}/SALTY_idFix &>> ${TMP}/RPlog.log
${SOFTWARE}/plink --bfile ${TMP}/TwinGene_idFix --bmerge ${TMP}/SALTY_idFix --make-bed --out ${TMP}/TwinGene-SALTY_merged &>> ${TMP}/RPlog.log    #THIS WILL FAIL BUT OUTPUT A FILE FOR FILTERING 'POOR' SNPS

${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene_idFix --exclude ${TMP}/TwinGene-SALTY_merged-merge.missnp --make-bed --out ${TMP}/TwinGene_idFix_rmpoor &>> ${TMP}/RPlog.log    #CUT THOSE IDENTIFIED ABOVE
${SOFTWARE}/plink2 --bfile ${TMP}/SALTY_idFix --exclude ${TMP}/TwinGene-SALTY_merged-merge.missnp --make-bed --out ${TMP}/SALTY_idFix_rmpoor &>> ${TMP}/RPlog.log
${SOFTWARE}/plink --bfile ${TMP}/TwinGene_idFix_rmpoor --bmerge ${TMP}/SALTY_idFix_rmpoor --make-bed --out ${TMP}/TwinGene-SALTY_merged2 &>> ${TMP}/RPlog.log

cut -f 2 ${TMP}/TwinGene_idFix_rmpoor.bim > ${TMP}/TwinGene_SNPs.txt
cut -f 2 ${TMP}/SALTY_idFix_rmpoor.bim > ${TMP}/SALTY_SNPs.txt
grep -f ${TMP}/TwinGene_SNPs.txt ${TMP}/SALTY_SNPs.txt > ${TMP}/TwinGene-SALTY_SNPs.txt    #GRAB SNPS IN THE INTERSECTION OF TWINGENE AND SALTY
${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene-SALTY_merged2 --extract ${TMP}/TwinGene-SALTY_SNPs.txt --make-bed --out ${TMP}/TwinGene-SALTY_merged2_intersect &>> ${TMP}/RPlog.log    #EXCLUDE NON-OVERLAPPING SNPS

echo "PERFORMING MERGING OF EIRA-SRQB WITH SALTY AND TWINGENE DATA"

${SOFTWARE}/plink2 --bfile ${TMP}/EIRA-SRQB_QC --set-all-var-ids @:# --make-bed --out ${TMP}/EIRA-SRQB_idFix &>> ${TMP}/RPlog.log    #CHANGE IDS TO REMOVE DEPENDENCE ON ALLELES
${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene-SALTY_merged2_intersect --set-all-var-ids @:# --make-bed --out ${TMP}/TwinGene-SALTY_merged2_intersect_idFix &>> ${TMP}/RPlog.log    #REMOVE RS-NUMBERS TO MAKE ALL COMPARABLE
${SOFTWARE}/plink --bfile ${TMP}/TwinGene-SALTY_merged2_intersect_idFix --bmerge ${TMP}/EIRA-SRQB_idFix --make-bed --out ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged &>> ${TMP}/RPlog.log

${SOFTWARE}/plink --bfile ${TMP}/TwinGene-SALTY_merged2_intersect_idFix --exclude ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged-merge.missnp --make-bed --out ${TMP}/TwinGene-SALTY_merged2_intersect_idFix_rmpoor &>> ${TMP}/RPlog.log
${SOFTWARE}/plink --bfile ${TMP}/EIRA-SRQB_idFix --exclude ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged-merge.missnp --make-bed --out ${TMP}/EIRA-SRQB_idFix_rmpoor &>> ${TMP}/RPlog.log
${SOFTWARE}/plink --bfile ${TMP}/TwinGene-SALTY_merged2_intersect_idFix_rmpoor --bmerge ${TMP}/EIRA-SRQB_idFix_rmpoor --make-bed --out ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged2 &>> ${TMP}/RPlog.log

cut -f 2 ${TMP}/TwinGene-SALTY_merged2_intersect_idFix_rmpoor.bim > ${TMP}/TwinGene-SALTY_SNPs2.txt
cut -f 2 ${TMP}/EIRA-SRQB_idFix_rmpoor.bim > ${TMP}/EIRA-SRQB_SNPs.txt
grep -f ${TMP}/TwinGene-SALTY_SNPs2.txt ${TMP}/EIRA-SRQB_SNPs.txt > ${TMP}/TwinGene-SALTY-EIRA-SRQB_SNPs.txt
${SOFTWARE}/plink2 --bfile ${TMP}/TwinGene-SALTY-EIRA-SRQB_merged2 --extract ${TMP}/TwinGene-SALTY-EIRA-SRQB_SNPs.txt --make-bed --out ${TMP}/REFPAN_MERGED &>> ${TMP}/RPlog.log
rm ${TMP}/EIRA-SRQB_idFix*; rm ${TMP}/TwinGene_idFix*; rm ${TMP}/SALTY_idFix*; rm ${TMP}/TwinGene-SALTY*;

### ### ### 3.1.4. POST-MERGING FINALIZATION QUALITY CONTROL

echo "PERFORMING POST-MERGING FINALIZATION QUALITY CONTROL"

${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_MERGED --geno 0.05 --mind 0.05 --maf 0.01 --make-bed --out ${TMP}/REFPAN_MERGED_GMM &>> ${TMP}/RPlog.log
${SOFTWARE}/plink --bfile ${TMP}/REFPAN_MERGED_GMM --hwe 1e-6 include-nonctrl --make-bed --out ${TMP}/REFPAN_MERGED_GMMH &>> ${TMP}/RPlog.log
${SOFTWARE}/plink --bfile ${TMP}/REFPAN_MERGED_GMMH --indep-pairwise 100 5 0.2 --out ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE &>> ${TMP}/RPlog.log
${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_MERGED_GMMH --extract ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE.prune.in --make-bed --out ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE &>> ${TMP}/RPlog.log
${SOFTWARE}/plink --bfile ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE --genome --min 0.10 --out ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE_REL &>> ${TMP}/RPlog.log
${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE --remove ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE_REL.genome --make-bed --out ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE_REL &>> ${TMP}/RPlog.log

PCA_INPUT=REFPAN_MERGED_GMMH_LDPRUNE_REL
for ITER in $(seq 1 5); do
	${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_MERGED_GMMH_LDPRUNE_REL --pca 10 approx --out pca --threads 8 &>> ${TMP}/RPlog.log
	Rscript SCRIPTS/MISC/pcafiltration.R
	${SOFTWARE}/plink2 --bfile ${TMP}/${PCA_INPUT} --remove pca_outliers.txt --out ${TMP}/${PCA_INPUT}_PCA${ITER} --make-bed &>> ${TMP}/RPlog.log
	rm pca*; PCA_INPUT=${PCA_INPUT}_PCA${ITER}
done
${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_MERGED_GMMH --keep ${TMP}/${PCA_INPUT}.fam --make-bed --out DATA/OUTPUT/DATA/REFPAN &>> ${TMP}/RPlog.log

#mv ${TMP}/RPLog.log DATA/OUTPUT/LOG/
#rm ${TMP}/*
echo "DONE! PLEASE FIND THE FINISHED REFERENCE PANEL AT 'DATA/OUTPUT/DATA/'"

### TO DO:
# 1.1. Storage and removal is not great right now, and we sometimes remove stuff that are kept for the checkpoints (if-checks that allow us to skip ahead) which make bug testing painful. Should be executed more precisely and more
#	carefully!
# 1.2. Add an IF-check that allows you to circumvent the merging if it already exists such files.
# 1.3. There are no safety checks built-in at the moment, meaning you need to go over the log-files (which are removed at the end of the script) if you wish to ascertain that everything worked out. You should implement some information
#	log that handles this better, so that you can check the removal of every individual step and make sure it looks OK... I have tested it carefully though, and since I've made no changes to data we should not have an issue with
#	this I believe.
### NOTES:
# 2.1.
