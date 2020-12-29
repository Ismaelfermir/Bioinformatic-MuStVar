#!/bin/bash
#	Author: Ismael Fernandez-Miranda
#	Modified year: 2020

echo -e '\nMake sure that your normal and tumor alignmet files are named as:\n'
echo -e '-  Sample1.normal.bam'
echo -e '-  Sample1.tumor.bam'

tool_directory=$( pwd )
echo -e '\nTool directory:'
sleep 2
cd ~
echo -e '\nCreating temporal files...'
mkdir ~/.temporal
mkdir ~/.temporal/VarScan2
mkdir ~/.temporal/VarScan2/filter
mkdir ~/.temporal/Strelka2
mkdir ~/.temporal/Strelka2/filter
mkdir ~/.temporal/Mutect2
mkdir ~/.temporal/Mutect2/filter
temporal=$( echo '~/.temporal' )
sleep 2

echo -e '\nCreating database...'
mkdir ~/.temporal/DB
mkdir ~/.temporal/output
sleep 2

echo -e '\nWork directory:'
sleep 2
echo $2
cd $2
sleep 2

echo -e '\nThreads:'
sleep 2
echo $1
sleep 2

echo -e '\nReference fasta:'
sleep 2
echo $3
sleep 2

echo -e '\nPanel bed file:'
sleep 2
echo $4
sleep 2

cp $4 ~/.temporal/

echo -e '\nSorting alingment files...\n'

for file in *.bam
do
	file_name=$( echo $file | sed 's:.bam::g' )
	echo -e 'File present. Sorting alignment for:\n'$file_name
	/home/tools/samtools-1.10/samtools sort -@ $1 $file -o $file_name'.sort.bam'
	echo -e 'Done!\n'
done

for file in *sort.bam
do
	file_name=$( echo $file | sed 's:.sort.bam::g' )
	echo -e 'File present. Indexing bam for:\n'$file_name
	/home/tools/samtools-1.10/samtools index $file $file.bai
	echo -e 'Done!\n'
done

echo -e '\nCalling SNVs and Indels with Mutect2...\n'

for file in *tumor.sort.bam
do
	file_name=$( echo $file | sed 's:tumor.sort.bam::g' )
	file_normal=$( echo $file | sed 's:tumor:normal:g' )
	output_name=$( echo $file | sed 's:.tumor.sort.bam::g' )
	echo -e 'Extracting the BAM header name for the normal sample:'
	normal_name=$( /home/tools/samtools-1.10/samtools view -H $file_normal | head -30 | grep 'SM:' | awk -F\: '{print $5}' | awk '{print $1}' )
	echo -e $normal_name
	echo -e '\nFile present. Mutect2:\n'$file_name
	/home/tools/gatk-4.1.9.0/gatk Mutect2 \
	-R $3 \
	-I $file \
	-I $file_normal \
	-normal $normal_name \
	--intervals $4 \
	-O $output_name'.somatic.vcf'
	echo -e 'Done!\n'
done

mkdir Mutect2
mv *.somatic.vcf* Mutect2

echo -e '\nLaunching annotation for Mutect2 SNVs and Indels...\n'
cd $2Mutect2

for file in *.somatic.vcf
do
	file_name=$( echo $file | sed 's:.somatic.vcf::g' )
	echo -e 'File present. Making annotation for:\n'$file_name
	perl /home/tools/annovar/table_annovar.pl --buildver hg38 --remove --outfile $file_name'_somatic' --protocol refGene,cytoBand,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_eur,avsnp150,cosmic70,clinvar_20190305 --operation g,r,f,f,f,f,f,f  --nastring '.' --vcfinput $file /home/tools/annovar/humandb38/
	echo -e 'Done!\n'
done

cp *multianno.txt ~/.temporal/Mutect2/

cd ~/.temporal/Mutect2

for file in *multianno.txt
do
	file_name=$( echo $file | sed 's:_snp.hg38_multianno.txt::g' )
	file_name=$( echo $file | sed 's:_indel.hg38_multianno.txt::g' )
	echo -e 'File present. Preprocessing for filter...\n'$file_name
	cut -d$'\t' -f 1-22 $file | sed 's:\t\.:\tNA:g' > $file'_completa.txt'
	cut -d$'\t' -f 24-35 $file | sed 's:\t\.:\tNA:g' > $file'_filtrada.txt'
	echo -e 'Done!\n'
done



echo -e '------------------------------------------------------'



echo -e '\nCalling SNVs and Indels with Strelka2...\n'
cd $2
mkdir Strelka2

for file in *tumor.sort.bam
do
	file_name=$( echo $file | sed 's:.tumor.sort.bam::g' )
	file_normal=$( echo $file | sed 's:tumor:normal:g' )
	directory=$( echo $file_name )
	cd $2Strelka2 && mkdir $directory
	echo -e 'File present. Preprocessing...'
	echo -e '\n'$file_name
	echo -e $2Strelka2/$directory
	sleep 30
	cd $2
	/home/tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam $file_normal --tumorBam $file --referenceFasta $3 --runDir $2Strelka2/$directory
	echo 'Starting with variant calling process...\n'$file_name
	cd $2Strelka2/$directory && ./runWorkflow.py -m local -j $1 --quiet
	cd $2Strelka2/$directory/results/variants
	gunzip somatic.snvs.vcf.gz
	gunzip somatic.indels.vcf.gz
	echo -e '\nPostprocessing file...\n'
	sleep 2
	sed '50i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' somatic.indels.vcf > somatic.indels.GT.vcf
	sed -ri 's|DP:DP2:|GT:DP:DP2:|g' somatic.indels.GT.vcf
	sed -ri 's|(:BCN50\t)|\10/0:|g' somatic.indels.GT.vcf
	sed -ri 's|(:BCN50\t[^\t]*\t)|\10/1:|g' somatic.indels.GT.vcf
	
	sed '50i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' somatic.snvs.vcf > somatic.snvs.GT.vcf
	sed -ri 's|DP:FDP:|GT:DP:FDP:|g' somatic.snvs.GT.vcf
	sed -ri 's|(:TU\t)|\10/0:|g' somatic.snvs.GT.vcf
	sed -ri 's|(:TU\t[^\t]*\t)|\10/1:|g' somatic.snvs.GT.vcf
	
	echo -e 'Done!\n'
	sleep 2
	
	echo -e '\nLauncing annotation for Strelka2 SNVs and Indels...\n'
	perl /home/tools/annovar/table_annovar.pl --buildver hg38 --remove --outfile $file_name"_snp" --protocol refGene,cytoBand,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_eur,avsnp150,cosmic70,clinvar_20190305 --operation g,r,f,f,f,f,f,f  --nastring "." --vcfinput somatic.snvs.GT.vcf /home/tools/annovar/humandb38/
	
	perl /home/tools/annovar/table_annovar.pl --buildver hg38 --remove --outfile $file_name"_indel" --protocol refGene,cytoBand,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_eur,avsnp150,cosmic70,clinvar_20190305 --operation g,r,f,f,f,f,f,f  --nastring "." --vcfinput somatic.indels.GT.vcf /home/tools/annovar/humandb38/
	
	cp *_snp.hg38_multianno.txt ~/.temporal/Strelka2
	echo -e 'Done!\n'
	sleep 2
done

cd ~/.temporal/Strelka2

for file in *multianno.txt
do
	file_name=$( echo $file | sed 's:_snp.hg38_multianno.txt::g' )
	file_name=$( echo $file | sed 's:_indel.hg38_multianno.txt::g' )
	echo -e 'File present. Preprocessing for filter...\n'$file_name
	cut -d$'\t' -f 1-22 $file | sed 's:\t\.:\tNA:g' > $file'_completa.txt'
	cut -d$'\t' -f 24-35 $file | sed 's:\t\.:\tNA:g' > $file'_filtrada.txt'
	echo -e 'Done!\n'
done



echo -e '\n------------------------------------------------------'



echo -e '\nCreating mpileup file from alignment...\n'

cd $2
mkdir VarScan2

for file in *tumor.sort.bam
do
	file_name=$( echo $file | sed 's:.tumor.sort.bam::g' )
	normal_file=$( echo $file | sed 's:tumor:normal:g' )
	echo -e 'File present...\n'$file_name
	/home/tools/samtools-1.10/samtools mpileup --fasta-ref $3 --positions $4 $normal_file $file > $file_name.mpileup
	echo -e 'Done!\n'
done


echo -e '\nCalling SNVs and Indels with VarScan2...\n'
cd $2
sleep 2

for file in *.mpileup
do
	file_name=$( echo $file | sed 's:.mpileup::g' )
	echo -e 'File present. Calling somatic variants...\n'$file_name
	java -jar /home/tools/VarScan-2.3.9/VarScan.v2.3.9.jar somatic $file --output-snp $file_name'.snp' --output-indel $file_name'.indel' --mpileup 1 --min-coverage 1 --min-coverage-normal 1 --min-coverage--turmor 1 --min-var-freq 0 --min-var-freq 0 --somatic-p-value 0.05 --min-avg-qual 20 --strand-filter 1 --validation 1
	echo -e 'Done!\n'
done

cp *indel ~/.temporal/VarScan2 && cp *snp ~/.temporal/VarScan2
mv *indel VarScan2 && mv *snp VarScan2
rm *mpileup




echo -e '\n\n\n------------------------------------------------------\n\n\n'



echo -e 'Filtering variant files with R...\n'

cd $tool_directory

/home/tools/R-3.6.3/bin/Rscript --vanilla $tool_directory/MuStVar.filter.R

/home/tools/R-3.6.3/bin/Rscript --vanilla $tool_directory/Compare.callers.R

echo -e '\nFinal mutational files are in:\n'
echo -e $2

cd ~/.temporal/output && cp * $2
cd ~/.temporal/ && cp --recursive DB/ $2





