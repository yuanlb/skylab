task DropSeqToolsDigitalExpression {
  File bam_input
  File whitelist

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.12"
  Int machine_mem_mb = 7500
  Int cpu = 1
  Int disk = ceil((size(bam_input, "G") + size(whitelist, "G")) * 4.0) + 1
  Int preemptible = 0

  meta {
    description: "Constructs a tab-delimited gzipped count matrix from a bam file containing reads marked with cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
  }

  parameter_meta {
    bam_input: "input bam file marked with cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
    whitelist: "10x genomics cell barcode whitelist for 10x V2. Only CB found in this list are included in the count matrix"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    DigitalExpression \
      SUMMARY=digital_expression_summary.txt \
      OUTPUT=digital_expression.txt.gz \
      CELL_BARCODE_TAG=CB \
      MOLECULAR_BARCODE_TAG=UB \
      GENE_EXON_TAG=GE \
      CELL_BC_FILE=${whitelist} \
      USE_STRAND_INFO=false \
      INPUT=${bam_input}
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File matrix_output = "digital_expression.txt.gz"
    File matrix_summary = "digital_expression_summary.txt"
  }
}


task CreateSparseCountMatrix {
  File bam_input
  File gtf_file

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.2"
  Int machine_mem_mb = 7500
  Int cpu = 1
  Int disk = ceil((size(bam_input, "G") + size(gtf_file, "G")) * 4.0) + 1
  Int preemptible = 0

  meta {
    description: "Constructs a compressed sparse row matrix from a bam file containing reads marked with cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
  }

  parameter_meta {
    bam_input: "input bam file marked with cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
    gtf_file: "the annotation file that was used to align the bam file passed to this function"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    CreateCountMatrix \
      --bam-file ${bam_input} \
      --output-prefix sparse_counts \
      --gtf-annotation-file ${gtf_file} \
      --cell-barcode-tag CB \
      --molecule-barcode-tag UB \
      --gene-id-tag GE

  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File sparse_count_matrix = "sparse_counts.npz"
    File row_index = "sparse_counts_row_index.npy"
    File col_index = "sparse_counts_col_index.npy"
  }
}

task MergeCountFiles {
  Array[File] sparse_count_matrices
  Array[File] row_indices
  Array[File] col_indices

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.2"
  Int machine_mem_mb = 7500
  Int cpu = 1
  Int disk = 20  # todo find out how to make this adaptive with Array[file] input
  Int preemptible = 0

  meta {
    description: "Constructs a compressed sparse row matrix by concatenating multiple input matrices"
  }

  parameter_meta {
    sparse_count_matrices: "array of count matrices in csr format, saved as numpy archives (.npy)"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    prefixes=$(python <<CODE
    matrices = ["${sep='", "' sparse_count_matrices}"]
    print(' '.join(m.replace('.npz', '') for m in matrices))
    CODE)

    MergeCountMatrices -o sparse_counts -i $prefixes
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File sparse_count_matrix = "sparse_counts.npz"
    File col_index = "sparse_counts_col_index.npy"
    File row_index = "sparse_counts_row_index.npy"
  }
}
