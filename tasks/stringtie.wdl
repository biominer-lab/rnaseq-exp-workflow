task stringtie {
    File bam
    File gtf
    String docker
    String sample_id
    String cluster
    String disk_size
    Int minimum_length_allowed_for_the_predicted_transcripts
    Int Junctions_no_spliced_reads
    Float minimum_isoform_abundance
    Float maximum_fraction_of_muliplelocationmapped_reads

    command <<<
	nt=$(nproc)
	mkdir ballgown
	/opt/conda/bin/stringtie -e \
                 -B \
                 -p $nt \
                 -f ${minimum_isoform_abundance} \
                 -m ${minimum_length_allowed_for_the_predicted_transcripts} \
                 -a ${Junctions_no_spliced_reads} \
                 -M ${maximum_fraction_of_muliplelocationmapped_reads} \
                 -G ${gtf} \
                 --rf \
                 -o ballgown/${sample_id}/${sample_id}.gtf \
                 -C ${sample_id}.cov.ref.gtf \
                 -A ${sample_id}.gene.abundance.txt \
                 ${bam}
	
    >>>
    
    runtime {
      docker: docker
      cluster: cluster
      systemDisk: "cloud_ssd 40"
      dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
    }
    
    output {
      File covered_transcripts = "${sample_id}.cov.ref.gtf"
      File gene_abundance = "${sample_id}.gene.abundance.txt"
      Array[File] ballgown = ["ballgown/${sample_id}/${sample_id}.gtf", "ballgown/${sample_id}/e2t.ctab", "ballgown/${sample_id}/e_data.ctab", "ballgown/${sample_id}/i2t.ctab", "ballgown/${sample_id}/i_data.ctab", "ballgown/${sample_id}/t_data.ctab"]
      File genecount = "${sample_id}_genecount.csv"
    }
}
