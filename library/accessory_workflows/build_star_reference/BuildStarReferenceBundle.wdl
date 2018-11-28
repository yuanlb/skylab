task BuildStarReference{
  File ref_fasta
  File gtf_file
  String ref_name
  Int machine_cpu_cores=16
  Int machine_mem_gb=50
  command {
    mkdir star
    STAR --runMode genomeGenerate \
      --genomeDir star \
      --genomeFastaFiles ${ref_fasta} \
      --sjdbGTFfile ${gtf_file} \
      --sjdbOverhang 100 \
      --runThreadN ${machine_cpu_cores} \
    tar -cvf "${ref_name}.tar" star
  }
  output {
    File starRef = "${ref_name}.tar"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0"
    memory: "${machine_mem_gb} GB"
    disks :"local-disk 100 HDD"
    cpu: ${machine_mem_gb}
  }
}


workflow StarRef {
  File fasta
  File gtf
  String ref_name

  call BuildStarReference {
    input:
      ref_fasta = fasta,
      gtf_file = gtf,
      ref_name = ref_name
  }
  output {
    File star_ref = BuildStarReference.starRef
  }
}
