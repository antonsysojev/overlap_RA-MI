#!/user/bin/env bash
### LAST VERSION UPDATE 21 AUG 2023 (v1.1.3) - ADDED A SINGLE LINE IN 3.3.3. WHICH PERFORMS AN ADDITIONAL ANALYSIS
### THIS SCRIPT OBTAINS ALL LDSC GENOME-WIDE GENETIC CORRELATION ESTIMATES FOR THE OVERLAP PROJECT

if ! [[ pwd == /home2/genetics/antobe/projects/OVERLAP_RA-CVD ]]; then cd /home2/genetics/antobe/projects/OVERLAP_RA-CVD; fi

TMP=TEMPORARY/tmp-3	#SET FILEPATHS
SOFTWARE=../../software
CACHE=TEMPORARY/CACHE
DATA=DATA/OUTPUT/DATA
LOG=DATA/OUTPUT/LOG
RAW=DATA/RAW
RES=DATA/OUTPUT/RES

if [[ -f ${LOG}/PLog3.log ]]; then rm ${LOG}/PLog3.log; fi; touch ${LOG}/PLog3.log

### ### ### 3.1. CONSTRUCTING REFERENCE PANEL

if ! [[ -f ${DATA}/REFPAN.bed ]] || ! [[ -f ${DATA}/REFPAN.bim ]] || ! [[ -f ${DATA}/REFPAN.fam ]]; then echo "CONSTRUCTING REFERENCE PANEL DATA..."; bash SCRIPTS/MISC/ConstructRefPan.sh; fi

### ### ### 3.2. REFERENCE PANEL LD SCORES

source activate ldsc
for CHR in $(seq 1 22); do if ! [[ -f ${CACHE}/LD${CHR}.l2.ldscore.gz ]] || ! [[ -f ${CACHE}/LD${CHR}.l2.M ]] || ! [[ -f ${CACHE}/LD${CHR}.l2.M_5_50 ]]; then

	awk -v OFS='\t' '{print $1, $2, $3}' ${RAW}/genpos_cM/chr${CHR}.OMNI.interpolated_genetic_map | awk -v var="${CHR}" '{print var":"$2, $3}' > ${TMP}/cM${CHR}.txt	#CLEAN CM DATA
	sort -u ${TMP}/cM${CHR}.txt > ${TMP}/cM${CHR}_nodup.txt		#SIMPLE FILTERING ON DUPLICATES, WORKS FOR NOW
	${SOFTWARE}/plink2 --bfile ${DATA}/REFPAN --chr ${CHR} --make-bed --out ${TMP}/REFPAN_${CHR} &>> ${LOG}/PLog3.log
	${SOFTWARE}/plink --bfile ${TMP}/REFPAN_${CHR} --update-cm ${TMP}/cM${CHR}_nodup.txt --make-bed --out ${TMP}/REFPAN_${CHR}_cM &>> ${LOG}/PLog3.log
	awk '{ if ($3 == 0) print $2}' ${TMP}/REFPAN_${CHR}_cM.bim > ${TMP}/cM${CHR}_SNPsMissing.txt		#ID AND REMOVE THOSE THAT HAD NO cM DATA AVAILABLE
	${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_${CHR}_cM --exclude ${TMP}/cM${CHR}_SNPsMissing.txt --make-bed --out ${TMP}/LD_REFPAN${CHR} &>> ${LOG}/PLog3.log

	${SOFTWARE}/ldsc/ldsc.py --bfile ${TMP}/LD_REFPAN${CHR} --l2 --ld-wind-cm 1 --out ${CACHE}/LD${CHR} &>> ${LOG}/PLog3.log
	echo "LD SCORES DONE FOR CHR ${CHR}..."

fi; done

### ### ### 3.2. ACQUIRING GWAS DATA

if ! [[ -f ${DATA}/RA_GWAS.fastGWA ]] || ! [[ -f ${DATA}/SPOS_RA_GWAS.fastGWA ]] || ! [[ -f ${DATA}/SNEG_RA_GWAS.fastGWA ]]; then echo "FAILED TO FIND RA, SPOS_RA OR SNEG_RA GWAS DATA... EXITING BASH"; exit 1; fi
if ! [[ -f ${RAW}/RA_DECODE_GWAS.txt ]] || ! [[ -f ${RAW}/SPOS_RA_DECODE_GWAS.txt ]] || ! [[ -f ${RAW}/SNEG_RA_DECODE_GWAS.txt ]]; then echo "FAILED TO FIND DECODE RA, SPOS_RA OR SNEG_RA GWAS DATA... EXITING BASH"; exit 1; fi
if ! [[ -f ${RAW}/MI_GWAS.txt ]]; then wget -O ${RAW}/MI_GWAS.txt ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST011001-GCST012000/GCST011365/GCST011365_buildGRCh37.tsv; fi
if ! [[ -f ${RAW}/CRP_GWAS.txt ]]; then wget -O ${RAW}/CRP_GWAS.txt.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029070/GCST90029070_buildGRCh37.tsv.gz; gunzip ${RAW}/CRP_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/WHR_GWAS.txt ]]; then wget -O ${RAW}/WHR_GWAS.txt.gz https://zenodo.org/record/1251813/files/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1; gunzip ${RAW}/WHR_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/BMI_GWAS.txt ]]; then wget -O ${RAW}/BMI_GWAS.txt.gz https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1; gunzip ${RAW}/BMI_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/SBP_GWAS.txt ]]; then wget -O ${RAW}/SBP_GWAS.txt.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz; gunzip ${RAW}/SBP_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/DBP_GWAS.txt ]]; then wget -O ${RAW}/DBP_GWAS.txt.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/Evangelou_30224653_DBP.txt.gz; gunzip ${RAW}/DBP_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/PP_GWAS.txt ]]; then wget -O ${RAW}/PP_GWAS.txt.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006629/Evangelou_30224653_PP.txt.gz; gunzip ${RAW}/PP_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/LDL_GWAS.txt ]]; then wget -O ${RAW}/LDL_GWAS.txt.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002412/GCST90002412_buildGRCh37.tsv.gz; gunzip ${RAW}/LDL_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/SMKINIT_GWAS.txt ]]; then wget -O ${RAW}/SMKINIT_GWAS.txt.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007474/SmokingInitiation.txt.gz; gunzip ${RAW}/SMKINIT_GWAS.txt.gz; fi
if ! [[ -f ${RAW}/SMKAMT_GWAS.txt ]]; then wget -O ${RAW}/SMKAMT_GWAS.txt.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007459/CigarettesPerDay.txt.gz; gunzip ${RAW}/SMKAMT_GWAS.txt.gz; fi

### ### ### 3.3. OBTAINING GENETIC CORRELATION ESTIMATES
### ### ### ### 3.3.1. GENETIC CORRELATION BETWEEN RA AND EXPLORATORY PHENOTYPES

echo "ACQUIRING GENETIC CORRELATION ESTIMATES..."

echo -e "SNP\tA1\tA2" > ${TMP}/REFPAN_SNPLIST.txt; awk -v OFS='\t' '{ print $1":"$4, $5, $6 }' ${DATA}/REFPAN.bim >> ${TMP}/REFPAN_SNPLIST.txt    #GET REFERENCE PANEL SNPS FOR ALIGNMENT

Rscript SCRIPTS/MISC/MungeHelper.R ${DATA}/RA_GWAS.fastGWA
${SOFTWARE}/ldsc/munge_sumstats.py --sumstats ${TMP}/RA_pMunge.txt --merge-alleles ${TMP}/REFPAN_SNPLIST.txt --out ${TMP}/RA_GWAS_MUNGE --chunksize 500000 &>> ${LOG}/PLog3.log     #SEE NOTES

for PHENO in MI CRP WHR BMI SBP DBP PP LDL SMKINIT SMKAMT; do if ! [[ -f ${RES}/RA-${PHENO}_ldsc.log ]]; then

	Rscript SCRIPTS/MISC/MungeHelper.R ${RAW}/${PHENO}_GWAS.txt
	${SOFTWARE}/ldsc/munge_sumstats.py --sumstats ${TMP}/${PHENO}_pMunge.txt --merge-alleles ${TMP}/REFPAN_SNPLIST.txt --out ${TMP}/${PHENO}_GWAS_MUNGE --chunksize 500000 &>> ${LOG}/PLog3.log

	# SEE NOTES FOR COMMENT ON SAMPLE OVERLAP
	${SOFTWARE}/ldsc/ldsc.py --rg ${TMP}/RA_GWAS_MUNGE.sumstats.gz,${TMP}/${PHENO}_GWAS_MUNGE.sumstats.gz --ref-ld-chr ${CACHE}/LD --w-ld-chr ${CACHE}/LD --out ${RES}/RA-${PHENO}_ldsc &>> ${LOG}/PLog3.log
	echo "GENETIC CORRELATIONS BETWEEN RA AND ${PHENO} SUCCESFULLY COMPUTED!"
	rm ${TMP}/${PHENO}_GWAS*

fi; done

### ### ### ### 3.3.2. GENETIC CORRELATION BETWEEN ALTERNATIVE RA PHENOTYPES AND MI

Rscript SCRIPTS/MISC/MungeHelper.R ${RAW}/MI_GWAS.txt
${SOFTWARE}/ldsc/munge_sumstats.py --sumstats ${TMP}/MI_pMunge.txt --merge-alleles ${TMP}/REFPAN_SNPLIST.txt --out ${TMP}/MI_GWAS_MUNGE --chunksize 500000 &>> ${LOG}/PLog3.log

for PHENO in SPOS_RA SNEG_RA RA_DECODE SPOS_RA_DECODE SNEG_RA_DECODE; do if ! [[ -f ${RES}/${PHENO}-MI_ldsc.log ]]; then

	if [[ ${PHENO} == *"DECODE"* ]]; then Rscript SCRIPTS/MISC/MungeHelper.R ${RAW}/${PHENO}_GWAS.txt; fi
	if ! [[ ${PHENO} == *"DECODE"* ]]; then Rscript SCRIPTS/MISC/MungeHelper.R ${DATA}/${PHENO}_GWAS.fastGWA; fi
	${SOFTWARE}/ldsc/munge_sumstats.py --sumstats ${TMP}/${PHENO}_pMunge.txt --merge-alleles ${TMP}/REFPAN_SNPLIST.txt --out ${TMP}/${PHENO}_GWAS_MUNGE --chunksize 500000 &>> ${LOG}/PLog3.log
	
	${SOFTWARE}/ldsc/ldsc.py --rg ${TMP}/${PHENO}_GWAS_MUNGE.sumstats.gz,${TMP}/MI_GWAS_MUNGE.sumstats.gz --ref-ld-chr ${CACHE}/LD --w-ld-chr ${CACHE}/LD --out ${RES}/${PHENO}-MI_ldsc &>> ${LOG}/PLog3.log
	echo "GENETIC CORRELATIONS BETWEEN ${PHENO} AND MI SUCCESFULLY COMPUTED!"
	rm ${TMP}/${PHENO}_GWAS*

fi; done

### ### ### ### 3.3.3. GENETIC CORRELATION BETWEEN RA AND MI, USING DIFFERENT REFERENCE PANEL DATA

for CHR in $(seq 1 22); do if ! [[ -f ${CACHE}/LD_noMHC${CHR}.l2.ldscore.gz ]] || ! [[ -f ${CACHE}/LD_noMHC${CHR}.l2.M ]] || ! [[ -f ${CACHE}/LD_noMHC${CHR}.l2.M_5_50 ]]; then

	echo "COMPUTING LD SCORES FROM REFERENCE PANEL DATA (WITHOUT THE MHC)... PLEASE HOLD..."
	if [[ ${CHR} == 6 ]]; then
		${SOFTWARE}/plink2 --bfile ${DATA}/REFPAN --chr 6 --from-bp 25643792 --to-bp 33424951 --make-bed --out ${TMP}/MHC &>> ${LOG}/PLog3.log    #EXTRACT THE MHC
		${SOFTWARE}/plink2 --bfile ${DATA}/REFPAN --exclude ${TMP}/MHC.bim --make-bed --out ${TMP}/REFPAN_no_MHC &>> ${LOG}/PLog3.log

		awk -v OFS='\t' '{print $1, $2, $3}' DATA/RAW/genpos_cM/chr${CHR}.OMNI.interpolated_genetic_map | awk -v var="${CHR}" '{print var":"$2, $3}' > ${TMP}/cM${CHR}.txt	#CLEAN CM DATA
		sort -u ${TMP}/cM${CHR}.txt > ${TMP}/cM${CHR}_nodup.txt		#SIMPLE FILTERING ON DUPLICATES, WORKS FOR NOW
		${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_no_MHC --chr ${CHR} --make-bed --out ${TMP}/REFPAN_${CHR} &>> ${LOG}/PLog3.log
		${SOFTWARE}/plink --bfile ${TMP}/REFPAN_${CHR} --update-cm ${TMP}/cM${CHR}_nodup.txt --make-bed --out ${TMP}/REFPAN_${CHR}_cM &>> ${LOG}/PLog3.log
		awk '{ if ($3 == 0) print $2}' ${TMP}/REFPAN_${CHR}_cM.bim > ${TMP}/cM${CHR}_SNPsMissing.txt		#ID AND REMOVE THOSE THAT HAD NO cM DATA AVAILABLE
		${SOFTWARE}/plink2 --bfile ${TMP}/REFPAN_${CHR}_cM --exclude ${TMP}/cM${CHR}_SNPsMissing.txt --make-bed --out ${TMP}/LD_REFPAN${CHR} &>> ${LOG}/PLog3.log
	
		${SOFTWARE}/ldsc/ldsc.py --bfile ${TMP}/LD_REFPAN${CHR} --l2 --ld-wind-cm 1 --out ${CACHE}/LD_noMHC${CHR} &>> ${LOG}/PLog3.log
	else 
		cp ${CACHE}/LD${CHR}.l2.ldscore.gz ${CACHE}/LD_noMHC${CHR}.l2.ldscore.gz
		cp ${CACHE}/LD${CHR}.l2.M ${CACHE}/LD_noMHC${CHR}.l2.M
		cp ${CACHE}/LD${CHR}.l2.M_5_50 ${CACHE}/LD_noMHC${CHR}.l2.M_5_50
	fi
	echo "LD SCORES DONE FOR CHR ${CHR}..."

fi; done

if ! [[ -f ${CACHE}/LDscore.1.l2.ldscore.gz ]] || ! [[ -f ${CACHE}/LD_1KG_1.l2.M ]] || ! [[ -f ${CACHE}/LD_1KG_1.l2.M_5_50 ]]; then    #CHECKS ONLY FOR CHR 1 DATA WHICH IS NOT ROBUST
	wget -O ${TMP}/1000G_Phase3_ldscores.tgz https://zenodo.org/record/7768714/files/1000G_Phase3_ldscores.tgz?download=1; tar xzvf ${TMP}/1000G_Phase3_ldscores.tgz -C ${TMP}/	#SEE NOTES FOR COMMENT ON DOWNLOAD
	for CHR in $(seq 22); do
		gunzip ${TMP}/LDscore/LDscore.${CHR}.l2.ldscore.gz
		head -n +1 ${TMP}/LDscore/LDscore.${CHR}.l2.ldscore > ${CACHE}/LDscore.${CHR}.l2.ldscore
		awk -v FS='\t' -v OFS='\t' '{print $1, $1":"$3, $3, $4}' ${TMP}/LDscore/LDscore.${CHR}.l2.ldscore | tail -n +2 >> ${CACHE}/LDscore.${CHR}.l2.ldscore
		cp ${TMP}/LDscore/LDscore.${CHR}.l2.M ${CACHE}/; cp ${TMP}/LDscore/LDscore.${CHR}.l2.M_5_50 ${CACHE}/
	done
fi

${SOFTWARE}/ldsc/ldsc.py --rg ${TMP}/RA_GWAS_MUNGE.sumstats.gz,${TMP}/MI_GWAS_MUNGE.sumstats.gz --ref-ld-chr ${CACHE}/LD_noMHC --w-ld-chr ${CACHE}/LD_noMHC --out ${RES}/RA-MI-noMHC_ldsc &>> ${LOG}/PLog3.log
${SOFTWARE}/ldsc/ldsc.py --rg ${TMP}/RA_GWAS_MUNGE.sumstats.gz,${TMP}/MI_GWAS_MUNGE.sumstats.gz --ref-ld-chr ${CACHE}/LDscore. --w-ld-chr ${CACHE}/LDscore. --out ${RES}/RA-MI-1KG_ldsc &>> ${LOG}/PLog3.log
if ! [[ -f ${TMP}/RA_DECODE_GWAS_MUNGE.sumstats.gz ]]; then    #NEEDED SINCE THE PREVIOUS PART IS ONLY DONE IF RESULTS DO NOT YET EXIST
	Rscript SCRIPTS/MISC/MungeHelper.R ${RAW}/RA_DECODE_GWAS.txt
	${SOFTWARE}/ldsc/munge_sumstats.py --sumstats ${TMP}/RA_DECODE_pMunge.txt --merge-alleles ${TMP}/REFPAN_SNPLIST.txt --out ${TMP}/RA_DECODE_GWAS_MUNGE --chunksize 500000 &>> ${LOG}/PLog3.log
fi
${SOFTWARE}/ldsc/ldsc.py --rg ${TMP}/RA_DECODE_GWAS_MUNGE.sumstats.gz,${TMP}/MI_GWAS_MUNGE.sumstats.gz --ref-ld-chr ${CACHE}/LD_noMHC --w-ld-chr ${CACHE}/LD_noMHC --out ${RES}/RA_DECODE-MI-noMHC_ldsc &>> ${LOG}/PLog3.log

### ### ### 3.4. CLEANING OF RESULTS

echo "CLEANING RESULTS AND PREPARING IN A READABLE FORMAT..."
Rscript SCRIPTS/MISC/CleanResLDSC.R
echo "ALL ESTIMATES OF GENOME-WIDE GENETIC CORRELATION SUCCESFULLY COMPUTED! RESULTS ARE AVAILABLE AT 'DATA/OUTPUT/RES/'"
rm ${TMP}/*

### TO DO:
# 1.1. Use `conda init` code instead of `source activate` for setting up the LDSC environment. It works now but the alternative should be more stable.
# 1.2. The cM file could be constructed within one single awk-step instead of the current pipes. This would increase speed.
# 1.3. I've simply cut all the duplicate IDs here somewhere, but this could be improved to keep copies from the duplicates though it requries a bit more effort in coding.
# 1.4. Computing LD-scores is somewhat slow. Parallel processing could increase speed substantially but I didn't want to bother with it at this point.
# 1.5. Still missing a fully coded, and working helper script `CleanResLDSC.R`. Need to write this once I get access to the actualy results files.
# 1.6. The CACHE folder is a mess filled to the brim with 3 * 4 * 22 files containing LD scores with non-intuitive names. A better solution is to let this folder host subfolders, say
#	LD, LD_noMHC and LD_1KG, which all contain their respective files (named the same way). Then navigating will be easier while the LDSC input strings are a bit longer.
# 1.7. The final part of 3.3.3. should be conditional on the files existing: there is no point in doing them if the result files are already available.
### NOTES
# 2.1. Note the use of the optional flag --chunksize. If we do not use it, it defaults to a larger number (I believe ten times the given chunksize). This makes it so that a single munging takes
#	ages (130 minutes and counting, per my testing). With this change (which comes from suggestions per the GitHub page for LDSC) computational time reduces to minutes.
# 2.2. Downloading the LDSC premade LD-scores is no longer intuitive. The authors have recently taken down the link, and moved it to a page where you can download the data for a cost. However, going
#	from their github page to the place where the files used to be hosted, finds a README with the description of their new home. It leads us onwards to a zenodo.org page for S-LDSC (the version for data with annotations)
#	which hosts multiple filesets. Going through these, we find LD-scores within the `1000G_Phase3_ldscores.tgz` file. Unfortunately, there are not that much information on what these data are,
#	but it seems to be made from 1000 Genomes Phase 3 data, at least.
# 2.3. You might be wondering why we do noting to correct for sample overlap between the two input data sets. In HDL, we do this quite carefully, and in LAVA we account for it too. But in LDSC, there is actually
#	no option outside of constraining the intercept, a praxis that has, with time, been reconsidered from how it was initially described in the original paper. Reading up on it on GitHub, people seem to 
#	advise AGAINST constraining intercept, and it does not seem to be the basic approach. Constraining intercept is also not as intuitive to do as simply stating a number of overlapping individuals, but we need
#	to specify the whole thing, computed manually, and it depends on number of individuals (a number that differs over SNPs) and the phenotypic correlation (which is estimated within the algorithm). Thus,
#	I think it is all right to ignore it here, as we don't really have any computational issues and our numbers 'make sense'.
