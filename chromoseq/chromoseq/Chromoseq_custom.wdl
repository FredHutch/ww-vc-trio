version 1.0

struct Samples {
  File Cram
  File CramIndex
  String Name
}


workflow ChromoSeq {
  input{
  # String Cram
  # String CramIndex 
  # String Name

  Array[Samples] MultipleSamples
  String Gender
  String MappingSummary
  String? CoverageSummary
  String OutputDir
  
  String Translocations
  String GenesBed
  
  String Cytobands
  String SVDB

  String CustomAnnotationVcf 
  String CustomAnnotationIndex 
  String CustomAnnotationParameters 

  String HotspotVCF
  String MantaConfig
  
  File Reference
  File ReferenceIndex
  File ReferenceBED
  File ReferenceBEDIndex
  String VEP

  File gcWig
  File mapWig
  String ponRds
  String centromeres
  String DupholdStatic

  String genomeStyle
  String genome

  String tmp
  
  Float minVarFreq = 0.02
  Int MinReads = 3
  Float varscanPvalindel = 0.1
  Float varscanPvalsnv = 0.01

  Int MinCNASize = 2000000
  Float MinCNAabund = 5.0
  Int LowCNASize = 50000000
  Float LowCNAabund = 10.0
  
  String chromoseq_docker

  File ichorToVCF
  File addReadCountsToVcfCRAM3

  }


scatter (sample in MultipleSamples){
      
  call prepare_bed {
    input:
      Bedpe=Translocations,
      Bed=GenesBed,
      Reference=ReferenceBED,
      tmp=tmp,
      docker=chromoseq_docker
  }

  call count_reads {
    input:
      Bam=sample.Cram,
      BamIndex=sample.CramIndex,
      list_chr=prepare_bed.chroms,
      ReferenceBED=ReferenceBED,
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      tmp=tmp,
      docker=chromoseq_docker,
  }
  
  call cov_qc as gene_qc {
    input:
      Cram=sample.Cram,
      CramIndex=sample.CramIndex,
      Name=sample.Name,
      Bed=GenesBed,
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      tmp=tmp,
      docker=chromoseq_docker
  }

  call cov_qc as sv_qc {
    input:
      Cram=sample.Cram,
      CramIndex=sample.CramIndex,
      Name=sample.Name,
      Bed=prepare_bed.svbed,
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      tmp=tmp,
      docker=chromoseq_docker
  }
  
  call run_manta {
    input:
      Bam=sample.Cram,
      BamIndex=sample.CramIndex,
      Config=MantaConfig,
      Reference=Reference,
      ReferenceIndex=ReferenceIndex,
      ReferenceBED=ReferenceBED,
      ReferenceBEDIndex=ReferenceBEDIndex,
      Name=sample.Name,
      tmp=tmp,
      docker=chromoseq_docker
    
  }

  call run_ichor {
    input:
      Bam=sample.Cram,
      BamIndex=sample.CramIndex,
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      ReferenceBED=ReferenceBED,
      CountFiles=count_reads.counts_bed,
      gender=Gender,
      gcWig=gcWig,
      mapWig=mapWig,
      ponRds=ponRds,
      centromeres=centromeres,
      Name=sample.Name,
      genomeStyle = genomeStyle,
      genome = genome,
      tmp=tmp,
      docker=chromoseq_docker 
  }
  
  call run_varscan {
    input:
      Bam=sample.Cram,
      BamIndex=sample.CramIndex,
      CoverageBed=GenesBed,
      MinFreq=minVarFreq,
      pvalsnv=varscanPvalsnv,
      pvalindel=varscanPvalindel,
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      Name=sample.Name,
      tmp=tmp,
      docker=chromoseq_docker,
      ReferenceIndex=ReferenceIndex
    
  }
  
  call run_pindel_region as run_pindel_flt3itd {
    input:
      Bam=sample.Cram,
      BamIndex=sample.CramIndex,
      Reg='chr13:28033987-28034316',
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      Name=sample.Name,
      genome=genome,
      tmp=tmp,
      docker=chromoseq_docker
    
  }

  call combine_variants {
    input:
      VCFs=[run_varscan.varscan_snv_file,
      run_varscan.varscan_indel_file, 
      run_pindel_flt3itd.pindel_vcf_file,
      HotspotVCF],
      Bam=sample.Cram,
      BamIndex=sample.CramIndex,
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      Name=sample.Name,
      MinReads=MinReads,
      MinVAF=minVarFreq,
      tmp=tmp,
      docker=chromoseq_docker,
      addReadCountsToVcfCRAM3=addReadCountsToVcfCRAM3
    
  }
  
  call annotate_variants {
    input:
      Vcf=combine_variants.combined_vcf_file,
      refFasta=Reference,
      ReferenceIndex=ReferenceIndex,
      Vepcache=VEP,
      Cytobands=Cytobands,
      CustomAnnotationVcf=CustomAnnotationVcf,
      CustomAnnotationIndex=CustomAnnotationIndex,
      CustomAnnotationParameters=CustomAnnotationParameters,
      Name=sample.Name,
      tmp=tmp,
      docker=chromoseq_docker
    
  }

  call annovar {
    input:
      combine_variants_vcf=combine_variants.combined_vcf_file,
      annotate_variants_vcf=annotate_variants.annotated_filtered_vcf, 
      genome=genome
  }
  
  # call annotate_svs {
  #   input:
  #     Vcf=run_manta.vcf,
  #     CNV=run_ichor.seg,
  #     refFasta=Reference,
  #     ReferenceIndex=ReferenceIndex,
  #     Vepcache=VEP,
  #     SVAnnot=SVDB,
  #     Translocations=Translocations,
  #     Cytobands=Cytobands,
  #     minCNAsize=MinCNASize,
  #     minCNAabund=MinCNAabund,
  #     lowCNAsize=LowCNASize,
  #     lowCNAabund=LowCNAabund,    
  #     Name=sample.Name,
  #     gender=Gender,
  #     tmp=tmp,
  #     docker=chromoseq_docker,
  #     ichorToVCF=ichorToVCF
    
  # }
  
  # call make_report {
  #   input:
  #     SVVCF=annotate_svs.vcf,
  #     GeneVCF=annotate_variants.annotated_filtered_vcf,
  #     KnownGenes=prepare_bed.genes,
  #     GeneQC=gene_qc.qc_out,
  #     SVQC=sv_qc.qc_out,
  #     MappingSummary=MappingSummary,
  #     CoverageSummary=CoverageSummary,
  #     Name=sample.Name,
  #     docker=chromoseq_docker,
  #     tmp=tmp
  # }

  }
  

  output {
    # Array[File] varscan_snv = run_varscan.varscan_snv_file
    # Array[File] varscan_indel = run_varscan.varscan_indel_file 
    # Array[File] pindel_flt3 = run_pindel_flt3itd.pindel_vcf_file
    Array[File] manta = run_manta.vcf
    Array[File] ichor_params = run_ichor.params
    Array[File] ichor_seg = run_ichor.seg
    Array[File] ichor_genomewide_pdf = run_ichor.genomewide_pdf
    Array[File] ichor_allgenomewide_pdf = run_ichor.allgenomewide_pdf
    Array[File] ichor_correct_pdf = run_ichor.correct_pdf
    # Array[File] ichor_rdata = run_ichor.rdata
    Array[File] ichor_wig = run_ichor.wig
    Array[File] combined_vcf = combine_variants.combined_vcf_file
    Array[File] annotate_variants_annotated_vcf = annotate_variants.annotated_vcf
    Array[File] annotate_variants_annotated_filtered_vcf = annotate_variants.annotated_filtered_vcf 
    Array[File] annovar_CombVar_vcf = select_all(annovar.annovar_CombVar_vcf)
    Array[File] annovar_CombVar_table = select_all(annovar.annovar_CombVar_table)
    Array[File] annovar_AnnotVar_vcf = select_all(annovar.annovar_AnnotVar_vcf)
    Array[File] annovar_AnnotVar_table = select_all(annovar.annovar_AnnotVar_table)

  }
  
}

task count_reads {
  input {
    File Bam
    File BamIndex
    File ReferenceBED
    File refFasta
    File ReferenceIndex
    String tmp
    String docker
    Array[String] list_chr
  }
  
  command <<<
    set -eo pipefail && \
    for Chrom in ~{sep=' ' list_chr}
    do
      (/usr/local/bin/bedtools makewindows -b ~{ReferenceBED} -w 500000 | awk -v OFS="\t" -v C="${Chrom}" '$1==C && NF==3' > ~{tmp}/${Chrom}.windows.bed) && \
      /usr/local/bin/samtools view -@ 5 -b -f 0x2 -F 0x400 -q 20 -T ~{refFasta} ~{Bam} ${Chrom} | /usr/local/bin/intersectBed -sorted -nobuf -c -bed -b stdin -a ~{tmp}/${Chrom}.windows.bed > ${Chrom}.counts.bed &
    done;
    wait
  >>>

  runtime {
    docker: docker
    cpu: "10"
    memory: "24 G"
  }
  output {
    Array[File] counts_bed =  ["chr1.counts.bed", "chr2.counts.bed", "chr3.counts.bed", "chr4.counts.bed", "chr5.counts.bed", "chr6.counts.bed", "chr7.counts.bed", "chr8.counts.bed", "chr9.counts.bed", "chr10.counts.bed",
                               "chr11.counts.bed", "chr12.counts.bed", "chr13.counts.bed", "chr14.counts.bed", "chr15.counts.bed", "chr16.counts.bed", "chr17.counts.bed", "chr18.counts.bed", "chr19.counts.bed", "chr20.counts.bed",
                               "chr21.counts.bed", "chr22.counts.bed", "chrX.counts.bed", "chrY.counts.bed"]
  }
}

task prepare_bed {
  input {  
    File Bedpe
    File Bed
    File Reference
    String tmp
    String docker
  }
  command <<<
    awk -v OFS="\t" '{ split($7,a,"_"); print $1,$2,$3,a[1],".",$9; print $4,$5,$6,a[2],".",$10; }' ~{Bedpe} | sort -u -k 1,1V -k 2,2n > sv.bed && \
    ((cat sv.bed | cut -f 4) && (cat ~{Bed} | cut -f 6)) > genes.txt && \
    gunzip -c ~{Reference} | cut -f 1 > chroms.txt
  >>>

  runtime {
    docker: docker
    cpu: "1"
    memory: "4 G"
  }

  output {
    File svbed = "sv.bed"
    File genes = "genes.txt"
    Array[String] chroms = read_lines("chroms.txt")
  }
}

task cov_qc {
  input {
    File Cram
    File CramIndex
    File Bed
    File refFasta
    File ReferenceIndex
    String Name
    String tmp
    String docker
  }

  String bed_basename = basename(Bed, ".bed")
  String covqc_out = Name + "." + bed_basename + ".covqc.txt"
  String region_dist_out = Name + ".mosdepth." + bed_basename + ".region.dist.txt"
  # String bed_basename = basename(Bed, ".bed")
  

  command <<<
    set -eo pipefail && \
    mosdepth -n -f "~{refFasta}" -t 4 -i 2 -x -Q 20 -b "~{Bed}" --thresholds 10,20,30,40 "~{Name}" "~{Cram}" && \
    bedtools intersect -header -b "~{Name}.regions.bed.gz" -a "~{Name}.thresholds.bed.gz" -wo | \
    awk -v OFS="\t" '{ if (NR==1){ print $0,"%"$5,"%"$6,"%"$7,"%"$8,"MeanCov"; } else { print $1,$2,$3,$4,$5,$6,$7,$8,sprintf("%.2f\t%.2f\t%.2f\t%.2f",$5/$NF*100,$6/$NF*100,$7/$NF*100,$8/$NF*100),$(NF-1); } }' > "~{covqc_out}" && \
    mv "~{Name}.mosdepth.region.dist.txt" "~{region_dist_out}"
  >>>
  
  runtime {
    docker: docker
    cpu: "6"
    memory: "32 G"
  }
  
  output {
    File global_dist = "${Name}.mosdepth.global.dist.txt"
    File qc_out = "${covqc_out}"
    File region_dist = "${region_dist_out}"
  }

}

task run_manta {
  input {
    File Bam
    File BamIndex 
    File Config
    File Reference
    File ReferenceIndex
    File ReferenceBED
    File ReferenceBEDIndex
    String Name
    String tmp
    String docker
  }

  command <<<
    set -eo pipefail && \
    /usr/local/src/manta/bin/configManta.py --config=~{Config} --tumorBam=~{Bam} --referenceFasta=~{Reference} --runDir=manta --callRegions=~{ReferenceBED} --outputContig && \
    ./manta/runWorkflow.py -m local -q research-hpc -j 4 -g 32 && \
    zcat ./manta/results/variants/tumorSV.vcf.gz | /bin/sed 's/DUP:TANDEM/DUP/g' > fixed.vcf && \
    /usr/local/bin/duphold_static -v fixed.vcf -b ~{Bam} -f ~{Reference} -t 4 -o ~{Name}.manta.tumorSV.vcf && \
    bgzip ~{Name}.manta.tumorSV.vcf && /usr/bin/tabix ~{Name}.manta.tumorSV.vcf.gz
  >>>
  runtime {
    docker: docker
    cpu: "6"
    memory: "32 G"
  }
  output {
    File vcf = "${Name}.manta.tumorSV.vcf.gz"
    File index = "${Name}.manta.tumorSV.vcf.gz.tbi"
  }
}


task run_ichor {
  input {
    File Bam
    File BamIndex
    File ReferenceBED
    Array[File] CountFiles
    File refFasta
    File ReferenceIndex
    String Name
    String gender
    String genome
    String genomeStyle
    File gcWig
    File mapWig
    File ponRds
    File centromeres
    
    Int? minCNAsize
    Float? lowAbundVal
    Int? lowAbundCnaSize
    
    String? tmp
    String docker
  }

  command <<<
    set -eo pipefail
    cat ~{sep=" " CountFiles} | sort -k 1,1V -k 2,2n | \
    awk -v window=500000 'BEGIN { chr=""; } { if ($1!=chr){ printf("fixedStep chrom=%s start=1 step=%d span=%d\n",$1,window,window); chr=$1; } print $4; }' > "~{Name}.ichor.tumor.wig" && \
    /usr/local/bin/Rscript /usr/local/bin/ichorCNA/scripts/runIchorCNA.R --id ~{Name} \
    --WIG "~{Name}.ichor.tumor.wig" --ploidy "c(2)" --normal "c(0.1,0.5,.85)" --maxCN 3 \
    --gcWig "~{gcWig}" \
    --mapWig "~{mapWig}" \
    --centromere "~{centromeres}" \
    --normalPanel "~{ponRds}" \
    --genomeBuild "~{genome}" \
    --sex "~{gender}" \
    --includeHOMD False --chrs "c(1:22, \"X\", \"Y\")" --chrTrain "c(1:22)" --fracReadsInChrYForMale 0.0005 \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --txnE 0.999999 --txnStrength 1000000 --genomeStyle ~{genomeStyle} --outDir ./ --libdir /usr/local/bin/ichorCNA/ && \
    awk -v G=~{gender} '$2!~/Y/ || G=="male"' "~{Name}.seg.txt" > "~{Name}.ichor.segs.txt" && \
    mv ~{Name}/*.pdf .
  >>>
  
  runtime {
    docker: docker
    cpu: "6"
    memory: "16 G"
  }
  
  output {
    File params = "${Name}.params.txt"
    File seg = "${Name}.ichor.segs.txt"
    File genomewide_pdf = "${Name}_genomeWide.pdf"
    File allgenomewide_pdf = "${Name}_genomeWide_all_sols.pdf"
    File correct_pdf = "${Name}_genomeWideCorrection.pdf"
    File rdata = "${Name}.RData"
    File wig = "${Name}.ichor.tumor.wig"
  }
}


task run_varscan {
  input {
    File Bam
    File BamIndex
    Int? MinCov
    Float? MinFreq
    Int? MinReads
    Float? pvalindel
    Float? pvalsnv
    File CoverageBed
    File refFasta
    File ReferenceIndex
    String Name
    String? tmp
    String docker
  }

  command <<<
    /usr/local/bin/samtools mpileup -f ~{refFasta} -l ~{CoverageBed} ~{Bam} > mpileup.out && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2snp mpileup.out --min-coverage ~{default=6 MinCov} --min-reads2 ~{default=3 MinReads} \
    --min-var-freq ~{default="0.02" MinFreq} --p-value ~{default="0.01" pvalsnv} --output-vcf | bgzip -c > ~{Name}.varscan.snv.vcf.gz && tabix ~{Name}.varscan.snv.vcf.gz && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2indel mpileup.out --min-coverage ~{default=6 MinCov} --min-reads2 ~{default=3 MinReads} \
    --min-var-freq ~{default="0.02" MinFreq} --p-value ~{default="0.1" pvalindel} --output-vcf | bgzip -c > ~{Name}.varscan.indel.vcf.gz && tabix ~{Name}.varscan.indel.vcf.gz
  >>>
  
  runtime {
    docker: docker
    cpu: "4"
    memory: "16 G"
  }
  output {
    File varscan_snv_file = "${Name}.varscan.snv.vcf.gz"
    File varscan_indel_file = "${Name}.varscan.indel.vcf.gz"
    File varscan_snv_file_idx = "${Name}.varscan.snv.vcf.gz.tbi"
    File varscan_indel_file_idx = "${Name}.varscan.indel.vcf.gz.tbi"

  }
}

task run_pindel_region {
  input {
    File Bam
    File BamIndex
    String Reg
    Int? Isize
    Int? MinReads
    File refFasta
    File ReferenceIndex
    String Name
    String? tmp
    String genome
    String docker
  }

  command <<<
    (set -eo pipefail && /usr/local/bin/samtools view -T ~{refFasta} ~{Bam} ~{Reg} | /opt/pindel-0.2.5b8/sam2pindel - ~{tmp}/in.pindel ~{default=250 Isize} tumor 0 Illumina-PairEnd) && \
    /usr/local/bin/pindel -f ~{refFasta} -p ~{tmp}/in.pindel -c ~{Reg} -o ~{tmp}/out.pindel && \
    /usr/local/bin/pindel2vcf -P ~{tmp}/out.pindel -G -r ~{refFasta} -e ~{default=3 MinReads} -R ~{default="hg38" genome} -d ~{default="hg38" genome} -v ~{tmp}/pindel.vcf && \
    /bin/sed 's/END=[0-9]*\;//' ~{tmp}/pindel.vcf | bgzip -c > ~{Name}.pindel.vcf.gz && tabix ~{Name}.pindel.vcf.gz
  >>>
  
  runtime {
    docker: docker
    cpu: "4"
    memory: "16 G"
  }
  output {
    File pindel_vcf_file = "${Name}.pindel.vcf.gz"
    File pindel_vcf_file_idx = "${Name}.pindel.vcf.gz.tbi"
  }
}

task combine_variants {
  input {
    Array[File] VCFs
    File Bam
    File BamIndex
    File refFasta
    File ReferenceIndex
    String Name
    Int MinReads
    Float MinVAF
    String? tmp
    String docker
    File addReadCountsToVcfCRAM3
  }

  command {
    /usr/bin/tabix ~{VCFs[0]}
    /usr/bin/tabix ~{VCFs[1]}
    /usr/bin/tabix ~{VCFs[2]}
    /usr/bin/tabix ~{VCFs[3]}
    /opt/conda/envs/python2/bin/bcftools merge --force-samples -O z ~{sep=" " VCFs} | \
    /opt/conda/envs/python2/bin/bcftools norm -m- -f ~{refFasta} -O z > ~{tmp}/combined.vcf.gz && /usr/bin/tabix -p vcf ~{tmp}/combined.vcf.gz
    /opt/conda/bin/python ~{addReadCountsToVcfCRAM3} -n ~{MinReads} -v ~{MinVAF} -r ~{refFasta} ~{tmp}/combined.vcf.gz ~{Bam} ~{Name} | \
    bgzip -c > ~{Name}.combined_variants_tagged.vcf.gz && \
    /usr/bin/tabix -p vcf ~{Name}.combined_variants_tagged.vcf.gz
  }
  runtime {
    docker: docker
    cpu: "4"
    memory: "12 G"
  }
  output {
    File combined_vcf_file = "${Name}.combined_variants_tagged.vcf.gz"
  }

}

task annotate_variants {
  input {
    File Vcf
    File refFasta
    File ReferenceIndex
    File Vepcache
    File Cytobands
    File CustomAnnotationVcf
    File CustomAnnotationIndex
    String CustomAnnotationParameters
    Float? maxAF
    String Name
    String? tmp
    String docker
  }

  command {
    set -eo pipefail && \
    /usr/bin/tabix ~{Cytobands} && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf --vcf --fasta ~{refFasta} --hgvs --symbol --term SO --per_gene -o ~{Name}.annotated_variants.vcf \
    -i ~{Vcf} --custom ~{Cytobands},cytobands,bed --custom ~{CustomAnnotationVcf},~{CustomAnnotationParameters} --offline --cache --max_af --dir ~{Vepcache} --cache_version 90 && \
    /opt/htslib/bin/bgzip -c ~{Name}.annotated_variants.vcf > ~{Name}.annotated_variants.vcf.gz && \
    /usr/bin/tabix -p vcf ~{Name}.annotated_variants.vcf.gz && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/ensembl-vep/filter_vep -i ~{Name}.annotated_variants.vcf.gz --format vcf -o ~{Name}.annotated_variants_filtered.vcf \
    --filter "(MAX_AF < ${default='0.001' maxAF} or not MAX_AF) or MYELOSEQ_TCGA_AC or MYELOSEQ_MDS_AC" && \
    /opt/htslib/bin/bgzip -c ~{Name}.annotated_variants_filtered.vcf > ~{Name}.annotated_variants_filtered.vcf.gz && \
    /usr/bin/tabix -p vcf ~{Name}.annotated_variants_filtered.vcf.gz
  }
  runtime {
    docker: docker
    cpu: "1"
    memory: "32 G"
  }
  output {
    File annotated_vcf = "${Name}.annotated_variants.vcf.gz"
    File annotated_filtered_vcf = "${Name}.annotated_variants_filtered.vcf.gz"
  }
}

task annovar {
  input {
    File combine_variants_vcf
    File annotate_variants_vcf
    String genome
  }

  String sampleName = basename(combine_variants_vcf, ".combined_variants_tagged.vcf.gz")
  String annovar_CombVar_out = "~{sampleName}.combined_variants"
  String annovar_AnnotVar_out = "~{sampleName}.annotated_variants"

  command <<<
    # set -eo pipefail && \
    if zgrep -v "^#" ~{combine_variants_vcf} | grep -q "."; then
      /annovar/table_annovar.pl ~{combine_variants_vcf} /annovar/humandb/ --buildver ~{genome} --outfile ~{annovar_CombVar_out} --remove --protocol refGene,knownGene,cosmic70,esp6500siv2_all,clinvar_20180603,gnomad211_exome --operation g,f,f,f,f,f --nastring . --vcfinput
    else
      touch "~{annovar_CombVar_out}.~{genome}_multianno.vcf" "~{annovar_CombVar_out}.~{genome}_multianno.txt" && echo "Combine variants file is empty"
    fi && \
    if zgrep -v "^#" ~{annotate_variants_vcf} | grep -q "."; then
      /annovar/table_annovar.pl ~{annotate_variants_vcf} /annovar/humandb/ --buildver ~{genome} --outfile ~{annovar_AnnotVar_out} --remove --protocol refGene,knownGene,cosmic70,esp6500siv2_all,clinvar_20180603,gnomad211_exome --operation g,f,f,f,f,f --nastring . --vcfinput
    else
      touch "~{annovar_AnnotVar_out}.~{genome}_multianno.vcf" "~{annovar_AnnotVar_out}.~{genome}_multianno.txt" && echo "Annotate variants file is empty"
    fi
  >>>

  runtime {
    docker: "getwilds/annovar:${genome}"
    memory: "4 GB"
    cpu: "2"
  }

  output{
    # File? annovar_CombVar_vcf = glob("${annovar_CombVar_out}.${genome}_multianno.vcf")
    # File? annovar_CombVar_table = glob("${annovar_CombVar_out}.${genome}_multianno.txt")
    # File? annovar_AnnotVar_vcf = glob("${annovar_AnnotVar_out}.${genome}_multianno.vcf")  
    # File? annovar_AnnotVar_table = glob("${annovar_AnnotVar_out}.${genome}_multianno.txt")
    File? annovar_AnnotVar_vcf = "${annovar_AnnotVar_out}.${genome}_multianno.vcf"
    File? annovar_AnnotVar_table = "${annovar_AnnotVar_out}.${genome}_multianno.txt"
    File? annovar_CombVar_vcf = "${annovar_CombVar_out}.${genome}_multianno.vcf"
    File? annovar_CombVar_table = "${annovar_CombVar_out}.${genome}_multianno.txt"
    
    # File annovar_manta_vcf = "${annovar_manta_out}.${genome}_multianno.vcf"
    # File annovar_manta_table = "${annovar_manta_out}.${genome}_multianno.txt"

  }

}


# task annotate_svs {
#   input {
#     File Vcf
#     File CNV
#     File refFasta
#     File ReferenceIndex
#     File Vepcache
#     String Name
#     String gender
#     File SVAnnot
#     File Translocations
#     File Cytobands
#     Int? minCNAsize
#     Float? minCNAabund
#     Int? lowCNAsize
#     Float? lowCNAabund
#     File ichorToVCF
    
#     String? tmp
#     String docker
#   }

#   command {
#     set -eo pipefail && \
#     # perl /usr/local/bin/ichorToVCF.pl -g ${gender} -minsize ${minCNAsize} -minabund ${minCNAabund} -lowsize ${lowCNAsize} -lowabund ${lowCNAabund} -r ${refFasta} ${CNV} | bgzip -c > cnv.vcf.gz && \
#     perl ~{ichorToVCF} -g ~{gender} -minsize ~{minCNAsize} -minabund ~{minCNAabund} -lowsize ~{lowCNAsize} -lowabund ~{lowCNAabund} -r ~{refFasta} ~{CNV} | bgzip -c > cnv.vcf.gz && \
#     /opt/htslib/bin/tabix -p vcf cnv.vcf.gz && \
#     /opt/conda/envs/python2/bin/bcftools query -l cnv.vcf.gz > name.txt && \
#     perl /usr/local/bin/FilterManta.pl -a ~{minCNAabund} -r ~{refFasta} -k ~{Translocations} ~{Vcf} filtered.vcf && \
#     /opt/conda/envs/python2/bin/svtools afreq filtered.vcf | \
#     /opt/conda/envs/python2/bin/svtools vcftobedpe -i stdin | \
#     /opt/conda/envs/python2/bin/svtools varlookup -d 200 -c BLACKLIST -a stdin -b ~{SVAnnot} | \
#     /opt/conda/envs/python2/bin/svtools bedpetovcf | \
#     /usr/local/bin/bedtools sort -header -g ~{ReferenceIndex} -i stdin | bgzip -c > filtered.tagged.vcf.gz && \
#     /opt/conda/envs/python2/bin/bcftools reheader -s name.txt filtered.tagged.vcf.gz > filtered.tagged.reheader.vcf.gz && \
#     /opt/htslib/bin/tabix -p vcf filtered.tagged.reheader.vcf.gz && \
#     /opt/conda/envs/python2/bin/bcftools concat -a cnv.vcf.gz filtered.tagged.reheader.vcf.gz | \
#     /usr/local/bin/bedtools sort -header -g ~{ReferenceIndex} -i stdin > svs.vcf && \
#     /opt/conda/envs/python2/bin/python /usr/local/src/manta/libexec/convertInversion.py /usr/local/bin/samtools ~{refFasta} svs.vcf | bgzip -c > svs.vcf.gz && \
#     /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl --format vcf --vcf --fasta ~{refFasta} --per_gene --symbol --term SO -o ~{Name}.svs_annotated.vcf -i svs.vcf.gz --custom ~{Cytobands},cytobands,bed --cache --dir ~{Vepcache} && \
#     /opt/htslib/bin/bgzip -c ~{Name}.svs_annotated.vcf > ~{Name}.svs_annotated.vcf.gz && \
#     /opt/htslib/bin/tabix -p vcf ~{Name}.svs_annotated.vcf.gz
#   }
  
#   runtime {
#     docker: docker
#     cpu: "1"
#     memory: "24 G"
#   }
  
#   output {
#     File vcf = "${Name}.svs_annotated.vcf.gz"
#     File vcf_index = "${Name}.svs_annotated.vcf.gz.tbi"
#   }
# }

# task make_report {
#   input {
#     File SVVCF
#     File GeneVCF
#     File KnownGenes
#     File MappingSummary
#     String? CoverageSummary
#     String SVQC
#     String GeneQC
#     String Name
#     String tmp
#     String docker
#     Int? MinGeneCov
#     Int? MinFracGene20
#     Int? MinRegionCov
#     Int? MinFracRegion10
#   }

#   command <<<
#     cat ${MappingSummary} ${CoverageSummary} | cut -d ',' -f 3,4 | sort -u > qc.txt && \
#     /opt/conda/bin/python /usr/local/bin/make_report3.py ${Name} ${GeneVCF} ${SVVCF} ${KnownGenes} "qc.txt" ${GeneQC} ${SVQC} > "${Name}.chromoseq.txt"
#   >>>
  
#   runtime {
#     docker: docker
#   }
  
#   output {
#     File report = "${Name}.chromoseq.txt"
#   }
  
# }

# task make_igv {
#   String Name
#   String docker
  
#   command {
#     cat <<EOF > ${Name}.igv.xml
#     <?xml version="1.0" encoding="UTF-8"?>
#     <Session genome="hg38" locus="All" version="3">
#     <Resources>
#     <Resource name="Structural variants" path="${Name}.svs_annotated.vcf.gz"/>
#     <Resource name="Gene variants" path="${Name}.annotated_filtered.vcf.gz"/>
#     <Resource name="Log2Ratio CN" path="${Name}.l2r.bw"/>
#     <Resource name="Copy Number Est." path="${Name}.cn.bw"/>
#     <Resource name="Copy Number Call" path="${Name}.cnv.bed"/>
#     <Resource name="Ensemble Genes" path="http://www.broadinstitute.org/igvdata/annotations/hg38/EnsemblGenes.ensGene"/>
#     </Resources>
#     </Session>
#     EOF    
#   }
  
#   runtime {
#     docker: docker
#   }
  
#   output {
#     File igv_xml = "${Name}.igv.xml"
#   }  
# }

