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
  File dbSNP_vcf_index
  File dbSNP_vcf
  File known_indels_sites_indices
  File known_indels_sites_VCFs
}


workflow mutation_calling {
  input {
	Array[File] tumorCram
	referenceGenome refGenome
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
        # input_bam_index = tumorMarkDuplicates.markDuplicates_bai,
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
  
  ## Convert unmapped CRAM -> fastq -> aligned SAM -> sorted aligned BAM

	input{
		File input_cram
		referenceGenome refGenome
    Int threads = 16
	}

	String base_file_name = basename(input_cram, ".cram")
	String read_group_id = "ID:" + base_file_name
	String sample_name = "SM:" + base_file_name
	String platform = "illumina"
	String platform_info = "PL:" + platform   # Create the platform information

  command <<<
  set -eo pipefail && \
  samtools fastq -@ ~{threads-1} ~{input_cram} | \
  bwa mem -p -v 3 -t ~{threads} -M -R '@RG\t~{read_group_id}\t~{sample_name}\t~{platform_info}' ~{refGenome.ref_fasta} - > ~{base_file_name}.sam && \
  # samtools view -bS -@ ~{threads-1} -o ~{base_file_name}.aligned.bam ~{base_file_name}.sam && \
  # samtools sort -@ ~{threads-1} -o ~{base_file_name}.sorted_query_aligned.bam ~{base_file_name}.aligned.bam
  samtools sort -@ ~{threads-1} -o ~{base_file_name}.sorted_query_aligned.bam ~{base_file_name}.sam
  >>>

  output {
    File analysisReadySorted = "~{base_file_name}.sorted_query_aligned.bam"
  }
  
  runtime {
    memory: "64 GB"
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
    cpu: 8
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
    # File input_bam_index
    referenceGenome refGenome
  }
  
  String base_file_name = basename(input_bam, ".duplicates_marked.bam")

  command <<<
    set -eo pipefail && \
    samtools index ~{input_bam} && \
    gatk --java-options "-Xms8g" BaseRecalibrator -R ~{refGenome.ref_fasta} -I ~{input_bam} -O ~{base_file_name}.recal_data.csv --known-sites ~{refGenome.dbSNP_vcf} --known-sites ~{refGenome.known_indels_sites_VCFs} && \
    gatk --java-options "-Xms8g" ApplyBQSR -bqsr ~{base_file_name}.recal_data.csv -I ~{input_bam} -O ~{base_file_name}.recal.bam -R ~{refGenome.ref_fasta} && \
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
    cpu: 6
    docker: "ghcr.io/getwilds/gatk:4.3.0.0"
  }
}
