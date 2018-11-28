version 1.0

task BuildStarReference{
  File ref_fasta
  File gtf_file
  String ref_name
  Int machine_cpu_cores=16
  Int machine_mem_gb=50
  command {

    ## Name of the input fasta file
    ref_fasta_inflated_filename=ref_fasta_inflated.fasta
    ## If the input reference is bzip2 compressed then decompress
    if [[ ~{ref_fasta} =~ \.bz2 ]]
    then
	bunzip2 -c ~{ref_fasta} > ${ref_fasta_inflated_filename}
    else 
    	mv ~{ref_fasta} ${ref_fasta_inflated_filename}
    fi

    ## Name of the input gtf file
    gtf_file_inflated_filename=gtf_ref_inflated.gtf
    ## If the input gtf annotaiton is bzip2 compressed then decompress
    if [[ ~{gtf_file} =~ \.bz2 ]]
    then
	bunzip2 -c ~{gtf_file} > ${gtf_file_inflated_filename}
    else 
    	 mv ~{gtf_file} ${gtf_file_inflated_filename}
    fi 

    ## Make output directory
    mkdir star

    ## Run STAR in genomeGenerate mode
    STAR --runMode genomeGenerate \
      --genomeDir star \
      --genomeFastaFiles ${ref_fasta_inflated_filename} \
      --sjdbGTFfile ${gtf_file_inflated_filename} \
      --sjdbOverhang 100 \
      --runThreadN ~{machine_cpu_cores} \

    ## Tar the index up
    tar -cvf "~{ref_name}.tar" star
  }
  output {
    File starRef = "~{ref_name}.tar"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0"
    memory: "~{machine_mem_gb} GB"
    disks :"local-disk 100 HDD"
    cpu: ~{machine_cpu_cores}
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
