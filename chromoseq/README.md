- `prepareBAM`: inital alignment, recalibration and merging of per sample CRAM files
- `chromoseq`: variant calling workflow
- `workflow_scripts`: helper scripts in chromoseq WDL
- `input_files`: static data files used for chromoseq
- `organize_output.sh`: group output files from chromoseq into per sample directory
- `analysis.ipynb`: generate the final output TSV file with annotation
- `04.14.csv`: gene panel annotation per sample

## prepareBAM workflow
```mermaid
graph LR
    A[Start] --> B[Scatter: cram_file]
    B --> C[BwaMem]
    C --> D[MarkDuplicates]
    D --> E[ApplyBaseRecalibrator]
    E --> F[mergeBAMsPerSample_renameReadGroup]
    F --> G["Output: sample.bam (with .bai)"]

    subgraph "For each cram_file"
        B
        C
        D
        E
    end
```

## Chromoseq_custom workflow
```mermaid
graph TD
    A[Start] --> B[prepare_bed]
    B --> C[count_reads]
    A --> D[gene_qc]
    A --> E[sv_qc]
    A --> F[run_manta]
    A --> G[run_ichor]
    A --> H[run_varscan]
    A --> I[run_pindel_flt3itd]
    
    H --> J[combine_variants]
    I --> J
    
    J --> K[annotate_variants]
    K --> L[annovar]
    J --> L
    
    C --> G
    K & L & F & G --> N[End]
```

### Helping code snippets
```
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
```