type="GU_finland_download/"
key="-WBC"

input_bed="/groups/wyattgrp/users/amunzur/chip_project/bed_variants/" # where the bed files originated from variants are kept
input_bam="/groups/wyattgrp/users/amunzur/chip_project/finland_bams/" # the bam files are kept here and are used to calculate the pileup
main_dir="/groups/wyattgrp/users/amunzur/chip_project/variant_bed_files/" # bed files based on called variants will be saved here
tmp_dir="/groups/wyattgrp/users/amunzur/chip_project/tmp/" # to save random tmp stuff

output=${main_dir}${type} # this depends on the type given by the user
mkdir -p ${output} # make it if it doesnt already exist
##########################################
# MAKE THE BED FILES 
##########################################
cd $input_bed
for file in *.tsv
do
    echo "${file}"
    cut -f -2 $file > "${output}$file" # get columns up to 2, which has chrom and start pos
    cut -f 2 $file > "${tmp_dir}col"

    paste "${output}$file" "${tmp_dir}col" > "${output}temp_${file}" # paste an extra col here
    # awk '{$3+=1; print;}' "${output}temp_${file}" > "${output}temp2_${file}" # add 1 to stop locations
    awk '{$2-=1; print;}' "${output}temp_${file}" > "${output}temp2_${file}" # subtract 1 from start locations
    tail -n +2 "${output}temp2_${file}" > "${output}RG_${file}"

done

cd $output
rm GU*
rm temp*

##########################################
# RUN MPILEUP - can run this for ctDNA or WBC bams by changing the keyword
##########################################
path_to_bams=${input_bam}${type}
cd $path_to_bams

for bam in RG_GU*.bam
do 
    echo "${bam}"
    sample=`echo "${bam}" | awk -F ${key} '{print $1}'` # get the sample name with some string cleaning
    echo $bed_name

    echo ${path_to_bams}${bam} > ${output}${sample}.bamList
    samtools mpileup -A -b ${output}${sample}.bamList -f /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa -l ${output}${sample}.tsv > ${output}${sample}.pileup
    
done

# a sample command:
# samtools mpileup -b /groups/wyattgrp/users/amunzur/chip_project/finland_bams/ctDNA_prognosis_ORIGINAL/GU-16-016-cfDNA-UMI-2016Dec12.bam \
# -f /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa \
# -l /groups/wyattgrp/users/amunzur/chip_project/variant_bed_files/RG_GU-16-016.tsv > tmp