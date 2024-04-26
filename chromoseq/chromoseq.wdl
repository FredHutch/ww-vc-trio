
workflow ChromoSeq {

    String Translocations
    String GenesBed
    String ReferenceBED
    String tmp
    String OutputDir

    call prepare_bed {
        input: 
            Bedpe=Translocations,
            Bed=GenesBed,
            Reference=ReferenceBED,
            tmp=tmp,
            OutputDir=OutputDir,
    }

    scatter (chr in prepare_bed.chroms){
    call count_reads {
      input: Bam=Cram,
      BamIndex=CramIndex,
      ReferenceBED=ReferenceBED,
      refFasta=Reference,
      refIndex=ReferenceIndex,
      Chrom=chr,
      tmp=tmp,
      docker=chromoseq_docker
    }   
  }

}


task prepare_bed {
  String Bedpe
  String Bed
  String Reference
  String? tmp
  String OutputDir

  command <<<
    awk -v OFS="\t" '{ split($7,a,"_"); print $1,$2,$3,a[1],".",$9; print $4,$5,$6,a[2],".",$10; }' ${Bedpe} | sort -u -k 1,1V -k 2,2n > sv.bed
    ((cat sv.bed | cut -f 4) && (cat ${Bed} | cut -f 6)) > genes.txt
    gunzip -c ${Reference} | cut -f 1 > chroms.txt
    # cp sv.bed genes.txt chroms.txt ${OutputDir}/
  >>>

  runtime {
    docker_image: "ubuntu:noble-20240423"
    cpu: "1"
    memory: "4 G"
  }

  output {
    File svbed = "sv.bed"
    File genes = "genes.txt"
    File chroms = "chroms.txt"
    # Array[String] chroms = read_lines("chroms.txt")
  }
}


task count_reads {
  String Bam
  String BamIndex
  String ReferenceBED
  String Chrom
  String jobGroup
  String refFasta
  String refIndex
  String tmp
  String docker
  
  command {
    set -eo pipefail
    bedtools makewindows -b ${ReferenceBED} -w 500000 | awk -v OFS="\t" -v C="${Chrom}" '$1==C && NF==3' > ${tmp}/${Chrom}.windows.bed
    samtools view -b -f 0x2 -F 0x400 -q 20 -T ${refFasta} ${Bam} ${Chrom} | intersectBed -sorted -nobuf -c -bed -b stdin -a ${tmp}/${Chrom}.windows.bed > ${Chrom}.counts.bed
  }

  runtime {
    docker_image: docker
    cpu: "1"
    memory: "8 G"
    job_group: jobGroup
  }
  output {
    File counts_bed = "${Chrom}.counts.bed"
  }
}


