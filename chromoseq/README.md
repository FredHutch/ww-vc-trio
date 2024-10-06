## Reference fasta file
/fh/fast/paguirigan_a/pub/ReferenceDataSets/genome_data/human/hg38/Gencode_GRCh38.primary_assembly.genome.fa

## Building Annovar Database
ml annovar/20200607-GCCcore-11.2.0-Perl-5.34.0

protocols=("refGene" "knownGene" "cosmic70" "esp6500siv2_all" "clinvar_20180603" "gnomad211_exome")
for prot in ${protocols[@]};do
	annotate_variation.pl --buildver hg38 --downdb --webfrom annovar $prot humandb/
done

table_annovar.pl RO20053_JR-WGS_230705_A00613_0568_BHW3JTDSX5.combined_variants_tagged.vcf.gz humandb/ \
	--buildver hg38 --out anno --remove \
	--protocol refGene,knownGene,cosmic70,esp6500siv2_all,clinvar_20180603,gnomad211_exome \
	--operation g,f,f,f,f,f --nastring . --vcfinput --nopolish
