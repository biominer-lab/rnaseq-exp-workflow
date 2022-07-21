task samtools {
    File sam
    String sample_id
    String bam = sample_id + ".bam"
    String sorted_bam = sample_id + ".sorted.bam"
    String percent_bam = sample_id + ".percent.bam"
    String sorted_bam_index = sample_id + ".sorted.bam.bai"
    String ins_size = sample_id + ".ins_size"
    String docker
    String cluster
    String disk_size
    Int insert_size

    command <<<
       set -o pipefail
       set -e
       /opt/conda/bin/samtools view -bS ${sam} > ${bam}
       /opt/conda/bin/samtools sort -m 1000000000 ${bam} -o ${sorted_bam}
       /opt/conda/bin/samtools index ${sorted_bam}
       /opt/conda/bin/samtools view -bs 42.1 ${sorted_bam} > ${percent_bam}
       /opt/conda/bin/samtools stats -i ${insert_size} ${sorted_bam} |grep ^IS|cut -f 2- > ${sample_id}.ins_size
    >>>

    runtime {
       docker: docker
       cluster: cluster
       systemDisk: "cloud_ssd 40"
       dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
    }

    output {
      File out_bam = sorted_bam
      File out_percent = percent_bam
      File out_bam_index = sorted_bam_index
      File out_ins_size = ins_size
    }

}

