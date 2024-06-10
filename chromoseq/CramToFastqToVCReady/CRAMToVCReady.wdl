version 1.0

struct referenceGenome {
	File ref_fasta
	File ref_fasta_index
	File ref_dict
	File ref_amb
	File ref_ann
	File ref_bwt
	File ref_pac
	File ref_sa
	String ref_name
}


workflow mutation_calling {
  input {
	Array[File] tumorCram
	referenceGenome refGenome
	File dbSNP_vcf
	File dbSNP_vcf_index
	File known_indels_sites_VCFs
	File known_indels_sites_indices

  }
 
  # Scatter for "tumor" samples   
  scatter (tumorSamples in tumorCram){

	call BwaMem as tumorBwaMem {
      input:
        input_cram = tumorSamples,
        refGenome = refGenome
    }
    
	call MarkDuplicates as tumorMarkDuplicates{
      input:
        input_bam = tumorBwaMem.analysisReadySorted
    }

	call ApplyBaseRecalibrator as tumorApplyBaseRecalibrator{
      input:
        input_bam = tumorMarkDuplicates.markDuplicates_bam,
        input_bam_index = tumorMarkDuplicates.markDuplicates_bai,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        refGenome = refGenome
      }
  }

output {
    # Array[File] tumoralignedBamSorted = tumorBwaMem.analysisReadySorted
    # Array[File] tumorMarkDuplicates_bam = tumorMarkDuplicates.markDuplicates_bam
    # Array[File] tumorMarkDuplicates_bai = tumorMarkDuplicates.markDuplicates_bai
    Array[File] tumoranalysisReadyBam = tumorApplyBaseRecalibrator.recalibrated_bam
    Array[File] tumoranalysisReadyIndex = tumorApplyBaseRecalibrator.recalibrated_bai
  }
}


# TASK DEFINITIONS

task BwaMem{
  
  ## Convert unmapped CRAM -> fastq -> aligned SAM -> aligned BAM -> sorted aligned BAM

	input{
		File input_cram
		referenceGenome refGenome
    Int threads = 24
	}

	String base_file_name = basename(input_cram, ".cram")
	String ref_fasta_local = basename(refGenome.ref_fasta)
	String read_group_id = "ID:" + base_file_name
	String sample_name = "SM:" + base_file_name
	String platform = "illumina"
	String platform_info = "PL:" + platform   # Create the platform information

  command <<<
  set -eo pipefail && \
  cp ~{refGenome.ref_fasta} . && \
  cp ~{refGenome.ref_fasta_index} . && \
  cp ~{refGenome.ref_dict} . && \
  cp ~{refGenome.ref_amb} . && \
  cp ~{refGenome.ref_ann} . && \
  cp ~{refGenome.ref_bwt} . && \
  cp ~{refGenome.ref_pac} . && \
  cp ~{refGenome.ref_sa} . && \
  samtools fastq -@ ~{threads-1} ~{input_cram} | \
  bwa mem -p -v 3 -t ~{threads} -M -R '@RG\t~{read_group_id}\t~{sample_name}\t~{platform_info}' ~{ref_fasta_local} - | \
  samtools view -bS -@ ~{threads-1} - | \
  samtools sort -@ ~{threads-1} -o ~{base_file_name}.sorted_query_aligned.bam -
  >>>

  output {
    File analysisReadySorted = "~{base_file_name}.sorted_query_aligned.bam"
  }
  
  runtime {
    memory: "48 GB"
    cpu: "~{threads}"
    docker: "ghcr.io/getwilds/bwa:0.7.17"
  }
}



# Mark duplicates (not SPARK, for some reason that does something weird)
task MarkDuplicates{
	input{
		File input_bam
	}

	String base_file_name = basename(input_bam, ".sorted_query_aligned.bam")
	String output_bam = "~{base_file_name}.duplicates_marked.bam"
	String output_bai = "~{base_file_name}.duplicates_marked.bai"
	String metrics_file = "~{base_file_name}.duplicate_metrics"

  command <<<
    gatk MarkDuplicates \
      --INPUT ~{input_bam} \
      --OUTPUT ~{output_bam} \
      --METRICS_FILE ~{metrics_file} \
      --CREATE_INDEX true \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
      --VALIDATION_STRINGENCY SILENT
  >>>

  runtime {
    docker: "ghcr.io/getwilds/gatk:4.3.0.0"
    memory: "48 GB"
    cpu: 16
  }

  output {
    File markDuplicates_bam = "~{output_bam}"
    File markDuplicates_bai = "~{output_bai}"
    File duplicate_metrics = "~{metrics_file}"
  }
}




task ApplyBaseRecalibrator{
  input{
    File input_bam
    File input_bam_index
    File dbSNP_vcf
    File dbSNP_vcf_index
    File known_indels_sites_VCFs
    File known_indels_sites_indices
    referenceGenome refGenome
  }
  
  String base_file_name = basename(input_bam, ".duplicates_marked.bam")
  String ref_fasta_local = basename(refGenome.ref_fasta)
  String dbSNP_vcf_local = basename(dbSNP_vcf)
  String known_indels_sites_VCFs_local = basename(known_indels_sites_VCFs)

  command <<<
    set -eo pipefail && \
    cp ~{refGenome.ref_fasta} . && \
    cp ~{refGenome.ref_fasta_index} . && \
    cp ~{refGenome.ref_dict} . && \
    cp ~{dbSNP_vcf} . && \
    cp ~{dbSNP_vcf_index} . && \
    cp ~{known_indels_sites_VCFs} . && \
    cp ~{known_indels_sites_indices} . && \
    samtools index ~{input_bam} && \
    gatk --java-options "-Xms8g" BaseRecalibrator -R ~{ref_fasta_local} -I ~{input_bam} -O ~{base_file_name}.recal_data.csv \
        --known-sites ~{dbSNP_vcf_local} --known-sites ~{known_indels_sites_VCFs_local} && \
    gatk --java-options "-Xms8g" ApplyBQSR -bqsr ~{base_file_name}.recal_data.csv -I ~{input_bam} \
        -O ~{base_file_name}.recal.bam -R ~{ref_fasta_local} && \
    # finds the current sort order of this bam file
    samtools view -H ~{base_file_name}.recal.bam | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > ~{base_file_name}.sortOrder.txt
  >>>

  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bai"
    File sortOrder = "~{base_file_name}.sortOrder.txt"
  }

  runtime {
    memory: "36 GB"
    cpu: 12
    docker: "ghcr.io/getwilds/gatk:4.3.0.0"
  }
}
