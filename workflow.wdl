import "./tasks/fastp.wdl" as fastp
import "./tasks/hisat2.wdl" as hisat2
import "./tasks/samtools.wdl" as samtools
import "./tasks/stringtie.wdl" as stringtie
import "./tasks/fastqc.wdl" as fastqc
import "./tasks/fastqscreen.wdl" as fastqscreen
import "./tasks/qualimapBAMqc.wdl" as qualimapBAMqc
import "./tasks/qualimapRNAseq.wdl" as qualimapRNAseq
import "./tasks/ballgown.wdl" as ballgown

workflow {{ project_name }} {
	File read1
	File read2
	File idx
	File screen_ref_dir
	File fastq_screen_conf
	File gtf
	String sample_id
	String fastp_docker
	String adapter_sequence
	String adapter_sequence_r2
	String fastp_cluster
	String umi_loc
	String idx_prefix
	String pen_intronlen
	String fastqc_cluster_config
	String fastqc_docker
	String fastqscreen_docker
	String fastqscreen_cluster_config
	String hisat2_docker
	String hisat2_cluster
	String qualimapBAMqc_docker
	String qualimapBAMqc_cluster_config
	String qualimapRNAseq_docker
	String qualimapRNAseq_cluster_config
	String samtools_docker
	String samtools_cluster
	String stringtie_docker
	String stringtie_cluster
	String multiqc_cluster_config
	String multiqc_docker
	Int multiqc_disk_size
	Int trim_front1 
	Int trim_tail1 
	Int max_len1 
	Int trim_front2 
	Int trim_tail2  
	Int max_len2 
	Int disable_adapter_trimming
	Int length_required
	Int umi_len
	Int UMI
	Int qualified_quality_phred
	Int length_required1
	Int disable_quality_filtering
	Int pen_cansplice
	Int pen_noncansplice
	Int min_intronlen
	Int max_intronlen
	Int maxins
	Int minins
	Int fastqc_disk_size
	Int fastqscreen_disk_size
	Int qualimapBAMqc_disk_size
	Int qualimapRNAseq_disk_size
	Int insert_size
	Int minimum_length_allowed_for_the_predicted_transcripts
	Int Junctions_no_spliced_reads
	Float minimum_isoform_abundance
	Float maximum_fraction_of_muliplelocationmapped_reads
	String ballgown_docker
	String ballgown_cluster
	String disk_size

	call fastp.fastp as fastp {
		input: 
		sample_id=sample_id, 
		read1 = read1, 
		read2 = read2,
		docker = fastp_docker,
		cluster = fastp_cluster,
		disk_size = disk_size,
		adapter_sequence = adapter_sequence,
		adapter_sequence_r2 = adapter_sequence_r2,
		umi_loc = umi_loc,
		trim_front1 = trim_front1,
		trim_tail1 = trim_tail1, 
		max_len1  = max_len1,
		trim_front2  = trim_front2,
		trim_tail2   = trim_tail2,
		max_len2  = max_len2,
		disable_adapter_trimming = disable_adapter_trimming,
		length_required = length_required,
		umi_len = umi_len,
		UMI = UMI,
		qualified_quality_phred = qualified_quality_phred,
		length_required1 = length_required1,
		disable_quality_filtering = disable_quality_filtering
		}

	call fastqc.fastqc as fastqc {
		input:
		read1 = fastp.Trim_R1, 
		read2 = fastp.Trim_R2,
		docker = fastqc_docker,
		cluster_config = fastqc_cluster_config,
		disk_size = fastqc_disk_size
	}

	call fastqscreen.fastq_screen as fastqscreen {
		input:
		read1 = fastp.Trim_R1, 
		read2 = fastp.Trim_R2,
		screen_ref_dir = screen_ref_dir,
		fastq_screen_conf = fastq_screen_conf,
		docker = fastqscreen_docker,
		cluster_config = fastqscreen_cluster_config,
		disk_size = fastqscreen_disk_size
	}

	call hisat2.hisat2 as hisat2 {
		input: 
		sample_id = sample_id, 
		idx = idx, 
		idx_prefix = idx_prefix, 
		Trim_R1 = fastp.Trim_R1, 
		Trim_R2 = fastp.Trim_R2,
		docker = hisat2_docker,
		cluster = hisat2_cluster,
		disk_size = disk_size,
		pen_intronlen = pen_intronlen,
		pen_cansplice = pen_cansplice,
		pen_noncansplice = pen_noncansplice,
		min_intronlen = min_intronlen,
		max_intronlen = max_intronlen,
		maxins = maxins,
		minins = minins
	}

	call samtools.samtools as samtools {
		input: 
		sample_id = sample_id, 
		sam = hisat2.sam,
		docker = samtools_docker,
		cluster = samtools_cluster,
		disk_size = disk_size,
		insert_size = insert_size
	}

	call stringtie.stringtie as stringtie {
		input: 
		sample_id = sample_id,
		gtf = gtf, 
		bam = samtools.out_bam,
		docker = stringtie_docker,
		cluster = stringtie_cluster,
		disk_size = disk_size,
		minimum_length_allowed_for_the_predicted_transcripts = minimum_length_allowed_for_the_predicted_transcripts,
		Junctions_no_spliced_reads = Junctions_no_spliced_reads,
		minimum_isoform_abundance = minimum_isoform_abundance,
		maximum_fraction_of_muliplelocationmapped_reads = maximum_fraction_of_muliplelocationmapped_reads
	}

	call ballgown.ballgown as ballgown {
		input: 
		sample_id = sample_id,
		docker = ballgown_docker,
		cluster = ballgown_cluster,
		ballgown = stringtie.ballgown,
		gene_abundance = stringtie.gene_abundance,
		disk_size = disk_size
	} 
			
	call qualimapBAMqc.qualimapBAMqc as qualimapBAMqc {
		input:
		bam = samtools.out_percent,
		docker = qualimapBAMqc_docker,
		cluster_config = qualimapBAMqc_cluster_config,
		disk_size = qualimapBAMqc_disk_size
	}

	call qualimapRNAseq.qualimapRNAseq as qualimapRNAseq {
		input:
		bam = samtools.out_percent,
		docker = qualimapRNAseq_docker,
		cluster_config = qualimapRNAseq_cluster_config,
		disk_size = qualimapRNAseq_disk_size,
		gtf = gtf
	}

}