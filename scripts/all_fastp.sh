qsub -N ABLATION scripts/fastp_SE.sh rawdata/srr_ablation_single rawdata/ablation clean/ablation
qsub -hold_jid ABLATION -N AMPUTATED scripts/fastp_PE.sh rawdata/srr_amputated_paired rawdata/amputated clean/amputated
qsub -hold_jid AMPUTATED -N CRYO scripts/fastp_SE.sh rawdata/srr_cryoinjury_single rawdata/cryoinjury clean/cryoinjury
qsub -hold_jid CRYO -N SHAM_PE scripts/fastp_PE.sh rawdata/srr_sham_paired rawdata/sham/PE clean/sham/PE
qsub -hold_jid SHAM_PE -N SHAM_SE scripts/fastp_SE.sh rawdata/srr_sham_single rawdata/sham/SE clean/sham/SE
qsub -hold_jid SHAM_SE -N UNINJURED_PE scripts/fastp_PE.sh rawdata/srr_uninjured_paired rawdata/uninjured/PE clean/uninjured/PE
qsub -hold_jid UNINJURED_PE -N UNINJURED_SE scripts/fastp_SE.sh rawdata/srr_uninjured_single rawdata/uninjured/SE clean/uninjured/SE