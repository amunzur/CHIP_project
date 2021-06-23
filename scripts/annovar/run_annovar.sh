#!/bin/bash

vcf_file=$1 # mutect result
type=$2 # are we dealing with ctDNA or WBC samples? given as a string, in the 2nd argument

if [[ "$type" == "ctDNA" ]]
then 
	output_file="/groups/wyattgrp/users/amunzur/chip_project/annovar_results/ctDNA_prognosis/$vcf_file.ANNOTATED"
	vcf_file_avinput="/groups/wyattgrp/users/amunzur/chip_project/annovar_inputs/ctDNA_prognosis/$vcf_file.avinput"
else
	output_file="/groups/wyattgrp/users/amunzur/chip_project/annovar_results/GU_finland_download/$vcf_file.ANNOTATED"
	vcf_file_avinput="/groups/wyattgrp/users/amunzur/chip_project/annovar_inputs/GU_finland_download/$vcf_file.avinput" 
fi

# convert to a compatible formate for annovar
perl /groups/wyattgrp/software/annovar/annovar/convert2annovar.pl -format vcf4 $vcf_file -outfile $vcf_file_avinput --withfreq
# run the annotation 
perl /groups/wyattgrp/software/annovar/annovar/table_annovar.pl $vcf_file_avinput /groups/wyattgrp/software/annovar/annovar/humandb/ -buildver hg38 -out $output_file -remove -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel -operation g,g,f,f,f,f,f,f,f,f,f,f,f -nastring .
