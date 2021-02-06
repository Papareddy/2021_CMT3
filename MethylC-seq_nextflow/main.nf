#!/usr/bin/env nextflow

/**********************
* parameters
**********************/
params.output = "./results"
params.files = "./*.fastq"
params.annot = "/groups/nodine/lab/members/Ranjith/annotations_ranj"


/********************
* set up channels
********************/
files = Channel.fromPath(params.files).map { file -> [ id:file.baseName,file:file] }
.ifEmpty { error "Cannot find any reads matching: ${params.files}" }



files.into{for_trim;for_fastQC}

/**********************
* fastqc
********************/
process fastqc {
    tag "$name"
    publishDir "${params.output}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(fq) from for_fastQC

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $fq
    """
}

/********************
* trim_galore
*********************/
process trim {
publishDir "$params.output/trimmed", mode: 'copy'
input:

  set val(id), file(fq) from for_trim

  output:
  set id, file("${id}_trimmed.fq") into bismark_align
  file("${id}.txt") into trimmed_log

script:
 """
  trim_galore --clip_R1 7 --three_prime_clip_R1 7 $fq 2> ${id}.txt
 """
}

/********************
* bismark align
*********************/
process align {
publishDir "$params.output/alignments", mode: 'copy'
  input:

  set id, file(fq) from bismark_align

  output:
  set id, file("*.bam") into dedup_bam
  set id, file("*report.txt") into alignment_log


  script:
   """
     bismark $params.annot/bismark_genome/ --bam  --score_min L,0,-0.2  --multicore 45 --se $fq
   """
  }

/********************
* bismark deduplicate
*********************/
process dedup {
publishDir "$params.output/deduplicate", mode: 'copy'
  input:

  set id, file(bam) from dedup_bam

  output:
  set id, file("*.bam") into tosort
  set id, file("*report.txt") into dup_log


  script:
   """
    deduplicate_bismark $bam
   """
  }


/********************
* samsort
*********************/

process samsort {
publishDir "$params.output/sortedbams", mode: 'copy'
  input:

  set id, file(bam) from tosort

  output:
  set id, file("*Chsort.bam") into allC

  script:
   """
    samtools sort -@ 8  $bam -o ${id}_Chsort.bam
   """
  }

/********************
* call_allC
*********************/

process call_allC {
  tag "$name"
      publishDir "$params.output/mpy", mode: 'copy',
          saveAs: {filename ->
              if (filename =~ '^allc' ) "methylpy/$filename"
              else if (filename =~ '^conversion' ) "info/$filename"
              else if (filename =~ '^log' ) "info/log.${name}.txt"
            }

 input:
  set val(name), file(bam) from allC

  output:
      file "*tsv" into for_CG,for_CHG,for_CHH
      file("conversion_rate_${name}.txt") into methylpy_conv_rate
      file("log.txt") into methylpy_log_file


  script:
      """
      export TMPDIR=\$(pwd)
      samtools index -b $bam
          methylpy call-methylation-state \
          --input-file $bam \
          --paired-end FALSE \
          --sample ${name} \
          --num-procs 16 \
	  --compress-output False \
          --unmethylated-control "chloroplast" --min-cov 1 \
          --binom-test False \
          --ref-fasta $params.annot/TAIR10.fa > log.txt 2>&1
    cat log.txt | grep "non-conversion rate" > conversion_rate_${name}.txt
      """
      }
