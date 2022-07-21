task qualimapRNAseq {
	File bam
	File gtf
	String bamname = basename(bam,".bam")
	String docker
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e
		nt=$(nproc)
		/opt/qualimap/qualimap rnaseq -bam ${bam} -outformat HTML -outdir ${bamname}_RNAseq -gtf ${gtf} -pe --java-mem-size=32G
		tar -zcvf ${bamname}_RNAseq_qualimap.tar.gz ${bamname}_RNAseq
		
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File rnaseq_zip = "${bamname}_RNAseq_qualimap.tar.gz"
	}
}