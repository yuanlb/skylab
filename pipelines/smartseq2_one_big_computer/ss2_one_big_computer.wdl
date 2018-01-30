task SmartSeq2 {

  # load annotation
  File gtf_file
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat

  #load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index

  # ref index name
  String hisat2_ref_name
  String hisat2_ref_trans_name

  # samples
  String stranded
  String sample_name
  String output_name
  File fastq1
  File fastq2

  command <<<
    # This combines all of the tasks of the current smartseq2 single cell
    # workflow into a single task. The purpose is to evaluate how per-task
    # overhead affects the performance of the workflow. 
    #
    # The SS2 workflow has six tasks
    # 1. hisat2.HISAT2PE 
    # 2.   picard.CollectMultipleMetrics
    # 3.   picard.CollectRnaMetrics
    # 4.   picard.CollectDuplicationMetrics
    # 5. hisat2.HISAT2rsem
    # 6.   rsem.RsemExpression
    #
    # There are three docker images that have to be combined for this to work:
    # 1. secondary-analysis-picard
    # 2. secondary-analysis-hisat2
    # 3. secondary-analysis-rsem
    # 
    # I don't think the source of these images is documented , but I'm
    # assuming they're in skylab/docker.
    
    set -u -x -e -o pipefail
    
    # We'll make use of fifos to minimize unnecessary disk IO 
    mkdir .fifos 
    
    # This is a common step for both hisat2 tasks
    if [[ "${fastq1}" != *.fastq.gz ]]; then
        FQ1="${fastq1}".fastq.gz
        mv "${fastq1}" "${fastq1}".fastq.gz
    else
        FQ1="${fastq1}"
    fi
    if [[ "${fastq2}" != *.fastq.gz ]]; then
        FQ2="${fastq2}".fastq.gz
        mv "${fastq2}" "${fastq2}".fastq.gz
    else
        FQ2="${fastq2}"
    fi
    
    # Extract the references for HISAT2 and RSEM
    tar xf "${hisat2_ref_trans_index}" &
    extract_trans_index_pid=$!
    tar xf "${rsem_ref_index}" &
    extract_rsem_index_pid=$!
    tar xf "${hisat2_ref_index}" &
    extract_index_pid=$!

    ##### 
    # Run hisat2.HISAT2rsem
    #####
    
    # We'll start rsem first because  
    # We need to wait for the extraction of the reference to be done before we
    # can start.
    wait "$extract_trans_index_pid" || true
    
    mkfifo .fifos/hisat2_rsem_output.sam
    hisat2 -t \
      -x "${hisat2_ref_trans_name}"/"${hisat2_ref_trans_name}" \
      -1 "$FQ1" \
      -2 "$FQ2" \
      --rg-id="${sample_name}" --rg SM:"${sample_name}" --rg LB:"${sample_name}" \
      --rg PL:ILLUMINA --rg PU:"${sample_name}" \
      --new-summary --summary-file "${output_name}"_rsem.log \
      --met-file "${output_name}"_rsem.hisat2.met.txt --met 5 \
      -k 10 \
      --mp 1,1 \
      --np 1 \
      --score-min L,0,-0.1 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      --seed 12345 \
      -p $(nproc) -S .fifos/hisat2_rsem_output.sam &
    
    # Sadly, RSEM doesn't seem to actually work reading from a named pipe or
    # stdin. For the former it just hangs, for the latter it complains about
    # the header (even if I add -h to samtools)
    samtools view -bS - < .fifos/hisat2_rsem_output.sam > "${output_name}_rsem.bam"
    
    #####
    ## Run rsem.RsemExpression
    #####
    
    wait "$extract_rsem_index_pid" || true
    rsem-calculate-expression \
      --bam \
      --paired-end \
       -p $(($(nproc) - 2)) \
      --time --seed 555 \
      --calc-pme \
      --single-cell-prior \
      "${output_name}_rsem.bam" \
      rsem/rsem_trans_index  \
      "${output_name}_rsem" &
   
    ##### 
    # Run hisat2.HISAT2PE
    #####
    wait "$extract_index_pid" || true
    # We don't use the SAM file produced by hisat2, so we'll pipe it directly
    # to samtools sort using this fifo
    mkfifo .fifos/hisat2_output.sam

    hisat2 -t \
      -x "${hisat2_ref_name}"/"${hisat2_ref_name}" \
      -1 "$FQ1" \
      -2 "$FQ2" \
      --rg-id="${sample_name}" --rg SM:"${sample_name}" --rg LB:"${sample_name}" \
      --rg PL:ILLUMINA --rg PU:"${sample_name}" \
      --new-summary --summary-file "${output_name}_qc.log" \
      --met-file "${output_name}_qc.hisat2.met.txt" --met 5 \
      --seed 12345 \
      -k 10 \
      --secondary \
      -p 1 -S .fifos/hisat2_output.sam &
    
    # Sort the output of hisat2. Write the results to another named pipe that
    # we'll use to distribute results to picard and to a SAM->BAM process
    mkfifo .fifos/sorted_hisat2_output.sam
    samtools sort -O SAM .fifos/hisat2_output.sam > .fifos/sorted_hisat2_output.sam &
    
    # Read from the sorted output and send it three places:
    # 1. A named pipe to picard.CollectMultipleMetrics
    # 2. A named pipe to picard.CollectRnaMetrics
    # 3. samtools view to make the BAM file
    #
    # Note that this working on the theory that the two picard processed can
    # work faster on piped, uncompressed input.
    # Also, MarkDuplicates cannot take streaming input, so we don't run it yet.
    # Finally, this is a useless cat but I can't remember where the & goes
    # otherwise
    mkfifo .fifos/collect_multiple_metrics_input.sam
    mkfifo .fifos/collect_rna_metrics_input.sam
    cat .fifos/sorted_hisat2_output.sam | \
        tee .fifos/collect_multiple_metrics_input.sam .fifos/collect_rna_metrics_input.sam |
        samtools view -bh - > "${output_name}_qc.bam" &
    bam_pid=$!
    
    ##### 
    # Run picard.CollectMultipleMetrics
    #####

    # This task and the subsequent one receive input streamed from samtools
    # sort.
    java -Xmx6g -jar /usr/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT=/dev/stdin \
      OUTPUT="${output_name}_qc" \
      FILE_EXTENSION=".txt" \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectGcBiasMetrics \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=QualityScoreDistribution \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=CollectQualityYieldMetrics \
      REFERENCE_SEQUENCE="${genome_ref_fasta}" \
      ASSUME_SORTED=true \
      < .fifos/collect_multiple_metrics_input.sam &
    
    ##### 
    # Run picard.CollectRnaMetrics
    #####
    java -Xmx3g -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT=/dev/stdin \
      OUTPUT="${output_name}_qc.rna_metrics.txt" \
      REF_FLAT="${gene_ref_flat}" \
      RIBOSOMAL_INTERVALS="${rrna_intervals}" \
      STRAND_SPECIFICITY=${stranded} \
      CHART_OUTPUT="${output_name}_qc.rna.coverage.pdf" \
      < .fifos/collect_rna_metrics_input.sam &
    

    ##### 
    # Run picard.MarkDuplicates
    #####
    
    # This can't start until the BAM is complete, so wait for that process.
    # Wait will fail if it's already done, but that's okay.
    wait "$bam_pid" || true

    java -Xmx6g -XX:ParallelGCThreads=2  -jar /usr/picard/picard.jar  MarkDuplicates \
       VALIDATION_STRINGENCY=SILENT  \
       INPUT="${output_name}_qc.bam" \
       OUTPUT="${output_name}_qc.MarkDuplicated.bam" \
       ASSUME_SORTED=true \
       METRICS_FILE="${output_name}_qc.duplicate_metrics.txt" \
       REMOVE_DUPLICATES=false &

    # Not totally sure what this is about, but it seems important.
    touch "${output_name}_qc.rna.coverage.pdf"

    wait
  >>>

  output {
    File aligned_bam = "${output_name}_qc.bam"
    File alignment_summary_metrics = "${output_name}_qc.alignment_summary_metrics.txt"
    File bait_bias_detail_metrics = "${output_name}_qc.bait_bias_detail_metrics.txt"
    File bait_bias_summary_metrics = "${output_name}_qc.bait_bias_summary_metrics.txt"
    File base_call_dist_metrics = "${output_name}_qc.base_distribution_by_cycle_metrics.txt"
    File base_call_pdf = "${output_name}_qc.base_distribution_by_cycle.pdf"
    File dedup_metrics = "${output_name}_qc.duplicate_metrics.txt"
    File error_summary_metrics = "${output_name}_qc.error_summary_metrics.txt"
    File gc_bias_detail_metrics = "${output_name}_qc.gc_bias.detail_metrics.txt"
    File gc_bias_dist_pdf = "${output_name}_qc.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_name}_qc.gc_bias.summary_metrics.txt"
    File insert_size_hist = "${output_name}_qc.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_name}_qc.insert_size_metrics.txt"
    File hisat2_logfile = "${output_name}_qc.log"
    File hisat2_metfile = "${output_name}_qc.hisat2.met.txt"
    File pre_adapter_details_metrics = "${output_name}_qc.pre_adapter_detail_metrics.txt"
    File quality_by_cycle_metrics = "${output_name}_qc.quality_by_cycle_metrics.txt"
    File quality_by_cycle_pdf = "${output_name}_qc.quality_by_cycle.pdf"
    File quality_distribution_dist_pdf = "${output_name}_qc.quality_distribution.pdf"
    File quality_distribution_metrics = "${output_name}_qc.quality_distribution_metrics.txt"
    File rna_coverage = "${output_name}_qc.rna.coverage.pdf"
    File rna_metrics = "${output_name}_qc.rna_metrics.txt"

    File aligned_trans_bam = "${output_name}_rsem.bam"
    File hisat2tran_logfile = "${output_name}_rsem.log"
    File hisat2tran_metfile = "${output_name}_rsem.hisat2.met.txt"
    File rsem_cnt_log = "${output_name}_rsem.stat/${output_name}_rsem.cnt"
    File rsem_model_log = "${output_name}_rsem.stat/${output_name}_rsem.model"
    File rsem_theta_log = "${output_name}_rsem.stat/${output_name}_rsem.theta"
    File rsem_gene_results = "${output_name}_rsem.genes.results"
    File rsem_isoform_results = "${output_name}_rsem.isoforms.results"
    File rsem_time_log = "${output_name}_rsem.time"
  }

  runtime {
    docker:"marcusczi/ss2_combo"
    memory:"25 GB"
    disks: "local-disk 40 HDD"
    cpu: "4"
  }

}

workflow ss2wf {
  File gtf_file
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat

  #load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index

  # ref index name
  File hisat2_ref_name
  File hisat2_ref_trans_name

  # samples
  String stranded
  String sample_name
  String output_name
  File fastq1
  File fastq2

  call SmartSeq2 {
    input: 
      gtf_file = gtf_file,
      genome_ref_fasta = genome_ref_fasta,
      rrna_intervals = rrna_intervals,
      gene_ref_flat = gene_ref_flat,
      hisat2_ref_index = hisat2_ref_index,
      hisat2_ref_trans_index = hisat2_ref_trans_index,
      rsem_ref_index = rsem_ref_index,
      hisat2_ref_name = hisat2_ref_name,
      hisat2_ref_trans_name = hisat2_ref_trans_name,
      stranded = stranded,
      sample_name = sample_name,
      output_name = output_name,
      fastq1 = fastq1,
      fastq2 = fastq2
  }

  output {
    File aligned_bam = SmartSeq2.aligned_bam
    File alignment_summary_metrics = SmartSeq2.alignment_summary_metrics
    File bait_bias_detail_metrics = SmartSeq2.bait_bias_detail_metrics
    File bait_bias_summary_metrics = SmartSeq2.bait_bias_summary_metrics
    File base_call_dist_metrics = SmartSeq2.base_call_dist_metrics
    File base_call_pdf = SmartSeq2.base_call_pdf
    File dedup_metrics = SmartSeq2.dedup_metrics
    File error_summary_metrics = SmartSeq2.error_summary_metrics
    File gc_bias_detail_metrics = SmartSeq2.gc_bias_detail_metrics
    File gc_bias_dist_pdf = SmartSeq2.gc_bias_dist_pdf
    File gc_bias_summary_metrics = SmartSeq2.gc_bias_summary_metrics
    File insert_size_hist = SmartSeq2.insert_size_hist
    File insert_size_metrics = SmartSeq2.insert_size_metrics
    File hisat2_logfile = SmartSeq2.hisat2_logfile
    File hisat2_metfile = SmartSeq2.hisat2_metfile
    File pre_adapter_details_metrics = SmartSeq2.pre_adapter_details_metrics
    File quality_by_cycle_metrics = SmartSeq2.quality_by_cycle_metrics
    File quality_by_cycle_pdf = SmartSeq2.quality_by_cycle_pdf
    File quality_distribution_dist_pdf = SmartSeq2.quality_distribution_dist_pdf
    File quality_distribution_metrics = SmartSeq2.quality_distribution_metrics
    File rna_coverage = SmartSeq2.rna_coverage
    File rna_metrics = SmartSeq2.rna_metrics

    File aligned_trans_bam = SmartSeq2.aligned_trans_bam
    File hisat2tran_logfile = SmartSeq2.hisat2tran_logfile
    File hisat2tran_metfile = SmartSeq2.hisat2tran_metfile
    File rsem_cnt_log = SmartSeq2.rsem_cnt_log
    File rsem_model_log = SmartSeq2.rsem_model_log
    File rsem_theta_log = SmartSeq2.rsem_theta_log
    File rsem_gene_results = SmartSeq2.rsem_gene_results
    File rsem_isoform_results = SmartSeq2.rsem_isoform_results
    File rsem_time_log = SmartSeq2.rsem_time_log
  }
}
