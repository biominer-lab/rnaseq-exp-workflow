task multiqc {

	Array[File] read1_zip
	Array[File] read2_zip

	Array[File] txt1
	Array[File] txt2

	Array[File] bamqc_zip
	Array[File] rnaseq_zip

	String docker
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e
		mkdir -p /cromwell_root/tmp/fastqc
		mkdir -p /cromwell_root/tmp/fastqscreen
		mkdir -p /cromwell_root/tmp/bamqc
		mkdir -p /cromwell_root/tmp/rnaseq

		cp ${sep=" " read1_zip} ${sep=" " read2_zip} /cromwell_root/tmp/fastqc
		cp ${sep=" " txt1} ${sep=" " txt2} /cromwell_root/tmp/fastqscreen
		for i in ${sep=" " bamqc_zip}
		do
		  tar -zxvf $i -C /cromwell_root/tmp/bamqc
		done
		
		for i in ${sep=" " rnaseq_zip}
		do
		  tar -zxvf $i -C /cromwell_root/tmp/rnaseq
		done
		

		multiqc /cromwell_root/tmp/
		cat multiqc_data/multiqc_fastq_screen.txt > multiqc_fastq_screen.txt
		cat multiqc_data/multiqc_fastqc.txt > multiqc_fastqc.txt
		cat multiqc_data/multiqc_general_stats.txt > multiqc_general_stats.txt
		cat multiqc_data/multiqc_qualimap_bamqc_genome_results.txt > multiqc_qualimap_bamqc_genome_results.txt

	
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File multiqc_html = "multiqc_report.html"
		Array[File] multiqc_txt = glob("multiqc_data/*")
		File multiqc_fastq_screen = "multiqc_fastq_screen.txt"
		File multiqc_fastqc = "multiqc_fastqc.txt"
		File multiqc_general_stats = "multiqc_general_stats.txt"
		File bamqc_genome_results = "multiqc_qualimap_bamqc_genome_results.txt"
	}
}