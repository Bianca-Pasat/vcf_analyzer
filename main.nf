#!/usr/bin/env nextflow


// mount to container at last tool or transfer files to outputdir

ch_input = Channel.fromPath(params.input, checkIfExists: true)
(ch_vcf,ch_outvcf,ch_outvcfcopy,ch_vcfname1,ch_vcfname2) = ch_input.into(5)

//snpeff_annotate
if(params.snpeffdb=='snpEff_v4_3_hg38.zip'){
ch_snpeff= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpeff/snpEff_v4_3_hg38.zip"))
}
if(params.snpeffdb=='snpEff_v4_3_GRCh37.75.zip'){
ch_snpeff= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpeff/snpEff_v4_3_GRCh37.75.zip"))
}
if(params.snpeffdb=='snpEff_v4_3_GRCh38.86.zip'){
ch_snpeff= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpeff/snpEff_v4_3_GRCh38.86.zip"))
}
if(params.snpeffdb=='snpEff_v4_3_hg19.zip'){
ch_snpeff= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpeff/snpEff_v4_3_hg19.zip"))
}
if(params.snpeffdb=='snpEff_v4_3_mm10.zip'){
ch_snpeff= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpeff/snpEff_v4_3_mm10.zip"))
}

//snpsift_annotate
if(params.snpsiftdb1=='dbNSFP4.0b1a_hg38_complete.vcf.gz'){
ch_snpsiftdb1= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/dbNSFP4.0b1a_hg38_complete.vcf.gz"))
ch_snpsiftdb1_index= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/dbNSFP4.0b1a_hg38_complete.vcf.gz.tbi"))
}
if(params.snpsiftdb1=='dbNSFP4.0b1a_hg19_complete.vcf.gz'){
ch_snpsiftdb1= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/dbNSFP4.0b1a_hg19_complete.vcf.gz"))
ch_snpsiftdb1_index= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/dbNSFP4.0b1a_hg19_complete.vcf.gz.tbi"))
}
if(params.snpsiftdb2=='ExAC.r1.sites.vep.hg38.vcf.gz'){
ch_snpsiftdb2= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/ExAC.r1.sites.vep.hg38.vcf.gz"))
ch_snpsiftdb2_index= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/ExAC.r1.sites.vep.hg38.vcf.gz.tbi"))
}
if(params.snpsiftdb2=='ExAC.r1.sites.vep.hg19.vcf.gz'){
ch_snpsiftdb2= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/ExAC.r1.sites.vep.hg19.vcf.gz"))
ch_snpsiftdb2_index= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/ExAC.r1.sites.vep.hg19.vcf.gz.tbi"))
}

if(params.cosmicdb=='hg38_cosmic92_coding.vcf.gz'){
ch_cosmicdb= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/hg38_cosmic92_coding.vcf.gz"))
ch_cosmicdb_index= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/hg38_cosmic92_coding.vcf.gz.tbi"))
}
if(params.cosmicdb=='hg19_cosmic87_coding.vcf.gz'){
ch_cosmicdb= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/hg19_cosmic87_coding.vcf.gz"))
ch_cosmicdb_index= Channel.value(file("s3://e-nios-bucket-1/data/VCF/snpsift/hg19_cosmic87_coding.vcf.gz.tbi"))
}



process inputvalidator{

    
    containerOptions "-v $params.inputdir:/$params.outputdir"
    input:
    
    val(vcf) from ch_vcf
    val(outvcf) from ch_outvcf
    val(outvcfcopy) from ch_outvcfcopy
    val(exte) from params.exte
    output:
    file("*_validated.vcf") into ch_out_vcf
    file("*_log.txt") into ch_out_log
    

    script:
    """
    python /code/python_scripts/inputfile_validator.py /data/${vcf.baseName}.vcf $exte && cp /data/${outvcf.baseName}.vcf ${vcf.baseName}_validated.vcf && cp *.txt /data/results/${vcf.baseName}_log.txt 
    """
}


process vcftools1{
    
    tag "${base}"
    
    
    input:
    
    file(vcf_vcftools) from ch_out_vcf
    val(pass) from params.pass
    val(recode) from params.recode
    val(recode_INFO_all) from params.recode_INFO_all
    val(minall) from params.minall
    val(maxall) from params.maxall
    val(gq) from params.gq
    val(dp) from params.dp
     
    output:
    file("*_vcftools_output.vcf.recode.vcf") into ch_out_vcftools
    file("*_vcftools_copy2.vcf") into ch_out_vcftools_copy
    
    script:
    """
    
    vcftools --vcf $vcf_vcftools $pass $recode $recode_INFO_all $minall $maxall $gq $dp --out ${vcf_vcftools.baseName}_vcftools_output.vcf && ls && cp *.vcf.recode.vcf ${vcf_vcftools.baseName}_vcftools_copy2.vcf
    """
}

// doesn't process nexr 
process filtervalidator{
    tag "${base}"
    
    
    
    input:
    
    file(vcffilter) from ch_out_vcftools 
    file(outvcffilter) from ch_out_vcftools_copy 
    output:
    file("*.txt") into ch_out_filter_log
    file("*.vcf") into ch_out_filter_vcf
    

    script:
    """
    
    python /code/python_scripts/vcftools_filtering_validator.py '$vcffilter'
    cp $outvcffilter ${vcffilter.baseName}_fv.vcf 
    mv *.txt ${vcffilter.baseName}.filtered.log.txt
    """
}



/// put more parameters
process snpeffan{
    
    
    input:
    file(snpefin) from ch_out_filter_vcf
    val(assembly) from params.assembly
    file(snpeff_database) from ch_snpeff
    output:
    file("*_snpeff_annotated.vcf") into ch_snpeff_out
    
    script:
    """
    unzip -o $snpeff_database -d /opt/conda/share/snpeff-4.3.1k-0
    snpEff 'ANN' -Xmx4096M -o vcf -nodownload -noLog -noStats -noLof -hgvs -sequenceOntology -no-downstream -no-intergenic -no-intron -no-upstream -no CDS -no CHROMOSOME_LARGE_DELETION -no CODON_CHANGE -no CODON_INSERTION -no CODON_CHANGE_PLUS_CODON_INSERTION -no CODON_DELETION -no CODON_CHANGE_PLUS_CODON_DELETION -no DOWNSTREAM -no EXON -no EXON_DELETED -no GENE -no INTERGENIC -no INTERGENIC_CONSERVED -no INTRAGENIC -no INTRON -no INTRON_CONSERVED -no MICRO_RNA -no NON_SYNONYMOUS_START -no NON_SYNONYMOUS_STOP -no RARE_AMINO_ACID -no TRANSCRIPT -no REGULATION -no UPSTREAM -no UTR_3_PRIME -no UTR_3_DELETED -no UTR_5_PRIME -no UTR_5_DELETED -no NEXT_PROT -verbose '$assembly' $snpefin > ${snpefin.baseName}_snpeff_annotated.vcf
    """
}


process snpsiftfilter{
    
    
    input:
    file(snpeff_annotated) from ch_snpeff_out
    val(effect_expression) from params.effect_expression
    output:
    file("*_snpeff_eff_significant.vcf") into ch_snpsift_filter_out
    
    script:
    """
    SnpSift -Xmx4096M filter "$effect_expression" -f $snpeff_annotated > ${snpeff_annotated.baseName}_snpeff_eff_significant.vcf && ls 
    """
}

(ch_split_var1,ch_split_var2)=ch_snpsift_filter_out.into(2)

process split_variants1{
    
    tag "${base}"
    
    
    input:
    
    file(vcf_indels) from ch_split_var1
   
    val(pass) from params.pass
    val(recode) from params.recode
    val(recode_INFO_all) from params.recode_INFO_all
    val(minall) from params.minall
    val(maxall) from params.maxall
    val(gq) from params.gq
    val(dp) from params.dp
    val(keep_indels) from params.keep_indels
    
     
    output:
    file("${vcf_indels.baseName}_indels.vcf.recode.vcf") into ch_indels_out
    file("${vcf_indels.baseName}_indels_copy.vcf") into ch_indels_out_c
    file("${vcf_indels.baseName}_indels_copy2.vcf") into ch_indels_out_c2
    
    
    script:
    """
    vcftools $pass $recode $recode_INFO_all $minall $maxall $gq $dp $keep_indels --out ${vcf_indels.baseName}_indels.vcf --vcf $vcf_indels && cp ${vcf_indels.baseName}_indels.vcf.recode.vcf ${vcf_indels.baseName}_indels_copy.vcf && cp ${vcf_indels.baseName}_indels.vcf.recode.vcf ${vcf_indels.baseName}_indels_copy2.vcf
    """
}

process split_variants2{
    
    tag "${base}"
    
    
    input:
    
    
    file(vcf_snps) from ch_split_var2
    val(pass) from params.pass
    val(recode) from params.recode
    val(recode_INFO_all) from params.recode_INFO_all
    val(minall) from params.minall
    val(maxall) from params.maxall
    val(gq) from params.gq
    val(dp) from params.dp
    
    val(keep_snps) from params.keep_snps
     
    output:
    
    file("${vcf_snps.baseName}_snps.vcf.recode.vcf") into ch_snps_out
    file("${vcf_snps.baseName}_snps_copy.vcf") into ch_snps_out_c
    script:
    """

    vcftools $pass $recode $recode_INFO_all $minall $maxall $gq $dp $keep_snps --out ${vcf_snps.baseName}_snps.vcf --vcf $vcf_snps && cp ${vcf_snps.baseName}_snps.vcf.recode.vcf ${vcf_snps.baseName}_snps_copy.vcf
    """
}

process indels_validator{
    
    
    tag "${base}"
    
    input:
    
    file(in_indels) from ch_indels_out 
    
    output:
    file("*.functional_indels_log.txt") into ch_indelsvalidator_log
    
    

    script:
    """
    
    python /code/python_scripts/effect_indels_validator.py '$in_indels'
    
    mv *.txt ${in_indels.baseName}.functional_indels_log.txt
    """
}

process snps_validator{
    
    tag "${base}"
    
    
    
    input:
    
    file(in_snps) from ch_snps_out 
    
    output:
    file("*.functional_snps_log.txt") into ch_snpsvalidator_log
    
    

    script:
    """
    
    python /code/python_scripts/effect_snps_validator.py '$in_snps'
    
    mv *.txt ${in_snps.baseName}.functional_snps_log.txt
    """
}

process snpsift_annotate{

    input:
    file(snpsift_indels) from ch_indels_out_c
    file(snpsift_indels2) from ch_indels_out_c2
    file(snpsiftdb1_path) from ch_snpsiftdb1
    file(snpsiftdb2_path) from ch_snpsiftdb2
    file(snpsiftdb1_path_index) from ch_snpsiftdb1_index
    file(snpsiftdb2_path_index) from ch_snpsiftdb2_index
    output:
    file("*_indels_snpsift_annotated.vcf") into ch_snpsift_anno_indels_out
    file("*_snps_snpsift_annotated.vcf") into ch_snpsift_anno_snps_out
    
    script:
    """
    SnpSift -Xmx4096M annotate -noDownload -noLog -name dbNSFP_ -info rs_dbSNP151,1000Gp3_AF,ESP6500_AA_AF,ESP6500_EA_AF,ExAC_AF,gnomAD_exomes_AF,SIFT_pred,\
Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,CADD_phred,clinvar_clnsig,phyloP100way_vertebrate,phastCons100way_vertebrate \
-db $snpsiftdb1_path $snpsift_indels > ${snpsift_indels.baseName}_snps_snpsift_annotated.vcf
SnpSift -Xmx4096M annotate -noDownload -noLog -a -info KG_AF_GLOBAL,ESP_AF_GLOBAL,AF \
-db $snpsiftdb2_path $snpsift_indels2 > ${snpsift_indels.baseName}_indels_snpsift_annotated.vcf
 
    """
}

process snpsift_fully_annotate{

    input:
    file(snpsift2_indels) from ch_snpsift_anno_indels_out
    file(snpsift2_snps) from ch_snpsift_anno_snps_out
    file(cosmicdb_path) from ch_cosmicdb
    file(cosmicdb_path_index) from ch_cosmicdb_index
    output:
    file("*_indels_fully_annotated.vcf") into ch_snpsift_full_anno_indels_out
    file("*_snps_fully_annotated.vcf") into ch_snpsift_full_anno_snps_out
    
    script:
    """
    SnpSift -Xmx4096M annotate -noDownload -noLog -db $cosmicdb_path $snpsift2_snps > ${snpsift2_snps.baseName}_snps_fully_annotated.vcf
    SnpSift -Xmx4096M annotate -noDownload -noLog -db $cosmicdb_path $snpsift2_indels > ${snpsift2_indels.baseName}_indels_fully_annotated.vcf
    """
}

(ch_fulanno_snps,ch_fulanno_snps_bio)=ch_snpsift_full_anno_snps_out.into(2)
(ch_fulanno_indels,ch_fulanno_indels_bio)=ch_snpsift_full_anno_indels_out.into(2)

process snpsift_filter{

    input:
    file(indels) from ch_fulanno_indels
    file(snps) from ch_fulanno_snps
    val(snpsexpr) from params.snpsexpr
    val(indelsexpr) from params.indelsexpr
    output:
    file("*_sig_indels.vcf") into ch_snpsift_sig_indels_out
    file("*_sig_snps.vcf") into ch_snpsift_sig_snps_out
    
    script:
    """
    SnpSift -Xmx4096M filter '$snpsexpr' -f '$snps' > ${snps.baseName}_sig_snps.vcf
    SnpSift -Xmx4096M filter '$indelsexpr' -f '$indels' > ${indels.baseName}_sig_indels.vcf

    """
}


(ch_sig_snps_log_in,ch_sig_snps_bio,ch_sig_snps)=ch_snpsift_sig_snps_out.into(3)
(ch_sig_indels_log_in,ch_sig_indels_bio,ch_sig_indels)=ch_snpsift_sig_indels_out.into(3)


process indels_fun_validator{
    tag "${base}"
    
    input:
    
    file(in_indels) from ch_sig_indels_log_in
    
    output:
    file("*.txt") into ch_sig_indels_log_out
    

    script:
    """
    
    python /code/python_scripts/indels_functional_filtering_validator.py '$in_indels'
    
    mv *.txt ${in_indels.baseName}.sig_indels_log.txt 
    """
}

process snps_fun_validator{
    tag "${base}"
    
    
    
    input:
    
    file(in_snps) from ch_sig_snps_log_in
    
    output:
    file("*.txt") into ch_sig_snps_log_out
    
    

    script:
    """
    
    python /code/python_scripts/snps_functional_filtering_validator.py '$in_snps'
    
    mv *.txt ${in_snps.baseName}.sig_snps_log.txt
    """
}



process snpsift_extract_indels{
    
    
    input:
    file(snpsift_indels) from ch_sig_indels
    output:
    file("${snpsift_indels.baseName}_indels_extracted.vcf") into ch_extracted_indels
    
    script:
    """
    SnpSift -Xmx4096M extractFields $snpsift_indels "CHROM" "POS" "REF" "ID" "ALT" "ANN" > ${snpsift_indels.baseName}_indels_extracted.vcf 
    """
}

process snpsift_extract_snps{
    
    
    input:
    file(snpsift_snps) from ch_sig_snps
    output:
    file("${snpsift_snps.baseName}_snps_extracted.vcf") into ch_extracted_snps
    
    script:
    """
    SnpSift -Xmx4096M extractFields $snpsift_snps "CHROM" "POS" "REF" "ALT" "dbNSFP_rs_dbSNP151" "dbNSFP_clinvar_clnsig" "ANN" > ${snpsift_snps.baseName}_snps_extracted.vcf 
    """
}



process final_tsvs{
    
    
    input:
    val(vcf) from ch_vcfname2
    file(indels_extracted) from ch_extracted_indels
    file(snps_extracted) from ch_extracted_snps
    
    output:
    file("*_indels_extracted.tsv") into ch_indels_tsv
    file("*_snps_extracted.tsv") into ch_snps_tsv
    file("*_Genes_BIMinput.txt") into ch_genes_bio
    
    script:
    """
    
    mv $snps_extracted ${vcf.baseName}_snps.vcf && 
    mv $indels_extracted ${vcf.baseName}_indels.vcf &&
    python /code/code2/final_outputs.py --flags /code/python_scripts/flags.txt --snps ${vcf.baseName}_snps.vcf --indels ${vcf.baseName}_indels.vcf --outsnps ${vcf.baseName}_snps_extracted.tsv --outindels ${vcf.baseName}_indels_extracted.tsv && cp Genes_BIMinput.txt ${vcf.baseName}_Genes_BIMinput.txt
    
    """
}

process prepare_bio{
    containerOptions "-v $params.inputdir:/$params.outputdir"
    publishDir "$params.outputdir/results", mode: 'copy'
    input:
    file(indels_final) from ch_indels_tsv
    file(snps_final) from ch_snps_tsv
    file(snps_fully_annotated) from ch_fulanno_snps_bio
    file(indels_fully_annotated) from ch_fulanno_indels_bio
    file(snps_significant) from ch_sig_snps_bio
    file(indels_significant) from ch_sig_indels_bio
    file(genes_bio) from ch_genes_bio
    val(vcf) from ch_vcfname1
    file(log1) from ch_out_log
    file(log2) from ch_out_filter_log
    file(log3) from ch_snpsvalidator_log //same as effect
    file(log4) from ch_indelsvalidator_log
    file(log5) from ch_sig_indels_log_out //same as functional 
    file(log6) from ch_sig_snps_log_out
    output:
    file("${vcf.baseName}_log.txt") 
    file("*_annotated_snps.vcf") 
    file("*_final_snps.vcf") 
    file("*_final_indels.vcf") 
    file("*_snps.tsv")
    file("*_indels.tsv") 
    file("*_genes_BIMinput.txt") 
    script:
    """
    
    cat $log1 $log2  $log3 $log4  $log5 $log6> ${vcf.baseName}.txt
    python /code/code2/name_handling.py --tag ${vcf.baseName} --genes $genes_bio --log ${vcf.baseName}.txt --indelstsv $indels_final --snpstsv $snps_final --snpsvcf $snps_fully_annotated --snps2vcf $snps_significant --indelsvcf $indels_fully_annotated --indels2vcf $indels_significant --outdir "" 
    
    """
}













