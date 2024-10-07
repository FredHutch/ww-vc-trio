version 1.0

# struct SampleIDs {
#     String Sample
# }

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
 
workflow prepareBAM {
    input{
        Map[String, Array[String]] tumorSamples
        Array[String] sample_ids
        referenceGenome refGenome
    }

    scatter (sample_id in sample_ids) {

        scatter (cram_file in tumorSamples[sample_id]) {
            call BwaMem {
                input:
                    input_cram = cram_file,
                    refGenome = refGenome,
                    sample = sample_id
            }

            call MarkDuplicates {
                input:
                    input_bam = BwaMem.analysisReadySorted
            }

            call ApplyBaseRecalibrator {
                input:
                    input_bam = MarkDuplicates.markDuplicates_bam,
                    input_bai = MarkDuplicates.markDuplicates_bai,
                    refGenome = refGenome
                }            
        }

        Array[File] bam_files = ApplyBaseRecalibrator.recalibrated_bam

        call mergeBAMsPerSample_renameReadGroup {
            input:
                bam_files = bam_files,
                sample = sample_id
        }

    }

    output {
        Array[File] final_bam = mergeBAMsPerSample_renameReadGroup.output_bam
        Array[File] final_bai = mergeBAMsPerSample_renameReadGroup.output_bai

    }

}
    




task BwaMem{
  
  ## Convert unmapped CRAM -> fastq -> aligned SAM -> sorted aligned BAM

	input{
		File input_cram
		referenceGenome refGenome
    Int threads = 20
    String sample
	}

	String base_file_name = basename(input_cram, ".cram")
    # String base_file_name = ~{sample}
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
  samtools sort -@ ~{threads-1} -o ~{base_file_name}.sorted_aligned.bam ~{base_file_name}.sam
  >>>

  output {
    File analysisReadySorted = "~{base_file_name}.sorted_aligned.bam"
  }
  
  runtime {
    memory: "64 GB"
    cpu: "~{threads}"
    docker: "getwilds/bwa:0.7.17"
  }
}


task MarkDuplicates{
	input{
		File input_bam
	}

	String base_file_name = basename(input_bam, ".sorted_aligned.bam")
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
    docker: "getwilds/gatk:4.3.0.0"
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
    File input_bai
    referenceGenome refGenome
  }
  
  String base_file_name = basename(input_bam, ".duplicates_marked.bam")

  command <<<
    set -eo pipefail && \
    gatk --java-options "-Xms8g" BaseRecalibrator -R ~{refGenome.ref_fasta} -I ~{input_bam} -O ~{base_file_name}.recal_data.csv --known-sites ~{refGenome.dbSNP_vcf} --known-sites ~{refGenome.known_indels_sites_VCFs} && \
    gatk --java-options "-Xms8g" ApplyBQSR -bqsr ~{base_file_name}.recal_data.csv -I ~{input_bam} -O ~{base_file_name}.recal.bam -R ~{refGenome.ref_fasta} --create-output-bam-index true && \
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
    docker: "getwilds/gatk:4.3.0.0"
  }
}


task mergeBAMsPerSample_renameReadGroup{
    input{
        Array[File] bam_files
        String sample
        Int num_threads = 23
    }

    command <<<
        samtools merge -f -@ ~{num_threads} ~{sample}.merged.bam ~{sep=' ' bam_files} && \
        gatk --java-options "-Xms14g" AddOrReplaceReadGroups -I ~{sample}.merged.bam -O ~{sample}.bam --CREATE_INDEX true -PL ILLUMINA -SM ~{sample} -LB lib1 -PU unit1
        #samtools addreplacerg -@ ~{num_threads} -r "@RG\tID:~{sample}\tSM:~{sample}\tPL:ILLUMINA" -o ~{sample}.aligned.merged.sameRG.bam -
    >>>

    output {
        File output_bam = "~{sample}.bam"
        File output_bai = "~{sample}.bai"

    }

    runtime {
        memory: "36 GB"
        cpu: "~{num_threads}"
        docker: "getwilds/gatk:4.3.0.0"
    }

}


task MarkDuplicates_ApplyBaseRecalibrator {
    input{
		File input_bam
    referenceGenome refGenome
	}

  String base_file_name = basename(input_bam, ".sorted_aligned.bam")
	String output_bam = "~{base_file_name}.duplicates_marked.bam"
	String output_bai = "~{base_file_name}.duplicates_marked.bai"
	String metrics_file = "~{base_file_name}.duplicate_metrics"

  command <<<
        set -eo pipefail && \
        gatk MarkDuplicates --INPUT ~{input_bam} --OUTPUT "~{base_file_name}.duplicates_marked.bam" --METRICS_FILE ~{metrics_file} --CREATE_INDEX true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --VALIDATION_STRINGENCY SILENT && \
        gatk --java-options "-Xms8g" BaseRecalibrator -R ~{refGenome.ref_fasta} -I ~{base_file_name}.duplicates_marked.bam -O ~{base_file_name}.recal_data.csv --known-sites ~{refGenome.dbSNP_vcf} --known-sites ~{refGenome.known_indels_sites_VCFs} && \
        gatk --java-options "-Xms12g" ApplyBQSR -bqsr ~{base_file_name}.recal_data.csv -I ~{base_file_name}.duplicates_marked.bam -O ~{base_file_name}.recal.bam -R ~{refGenome.ref_fasta} --create-output-bam-index true && \
        samtools view -H ~{base_file_name}.recal.bam | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > ~{base_file_name}.sortOrder.txt
  >>>

  output {

    File finalBAM = "~{base_file_name}.recal.bam"
    File finalBAI = "~{base_file_name}.recal.bai"

  }

  runtime{
    memory: "48 GB"
    cpu: 12
    docker: "getwilds/gatk:4.3.0.0"
  }

}
