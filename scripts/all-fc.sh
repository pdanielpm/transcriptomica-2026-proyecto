qsub -N fc-ablation scripts/fc-SE.sh alligned/ablation featurecounts/ablation_counts.txt Danio_rerio.GRCz11.115.chr.gtf
qsub -hold_jid fc-ablation -N fc-amputated scripts/fc-PE.sh alligned/amputated featurecounts/amputated_counts.txt Danio_rerio.GRCz11.115.chr.gtf
qsub -hold_jid fc-amputated -N fc-cryo scripts/fc-SE.sh alligned/cryoinjury featurecounts/cryoinjury_counts.txt Danio_rerio.GRCz11.115.chr.gtf
qsub -hold_jid fc-cryo -N fc-shamPE scripts/fc-PE.sh alligned/sham/PE featurecounts/shamPE_counts.txt Danio_rerio.GRCz11.115.chr.gtf
qsub -hold_jid fc-shamPE -N fc-shamSE scripts/fc-SE.sh alligned/sham/SE featurecounts/shamSE_counts.txt Danio_rerio.GRCz11.115.chr.gtf
qsub -hold_jid fc-shamSE -N fc-uninjuredPE scripts/fc-PE.sh alligned/uninjured/PE featurecounts/uninjuredPE_counts.txt Danio_rerio.GRCz11.115.chr.gtf
qsub -hold_jid fc-uninjuredPE -N fc-uninjuredSE scripts/fc-SE.sh alligned/uninjured/SE featurecounts/uninjuredSE_counts.txt Danio_rerio.GRCz11.115.chr.gtf