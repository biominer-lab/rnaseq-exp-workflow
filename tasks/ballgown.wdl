task ballgown {
    File gene_abundance
    Array[File] ballgown
    String sample_id
    String docker
    String cluster
    String disk_size

    command <<<
      mkdir -p /cromwell_root/tmp/${sample_id}
      cp -r ${sep=" " ballgown} /cromwell_root/tmp/${sample_id}
      ballgown /cromwell_root/tmp/${sample_id} ${sample_id}.txt
    >>>
    
    runtime {
      docker: docker
      cluster: cluster
      systemDisk: "cloud_ssd 40"
      dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
    }
    
    output {
      File mat_expression = "${sample_id}.txt"
    }
}