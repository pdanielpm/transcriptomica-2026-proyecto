qsub -N ablated_bgi scripts/STAR_SE.sh rawdata/ablated_bgi clean/ablation alligned/ablation STAR_Indices/50genomeDirectory
qsub -hold_jid ablated_bgi -N ablated_genome scripts/STAR_SE.sh rawdata/ablated_genome clean/ablation alligned/ablation STAR_Indices/50genomeDirectory
qsub -hold_jid ablated_genome -N amputated_hiseq  scripts/STAR_PE.sh rawdata/amputated_hiseq clean/amputated alligned/amputated STAR_Indices/150genomeDirectory
qsub -hold_jid amputated_hiseq -N amputated_novaseq  scripts/STAR_PE.sh rawdata/amputated_novaseq clean/amputated alligned/amputated STAR_Indices/150genomeDirectory
qsub -hold_jid amputated_novaseq -N cryo_bgi scripts/STAR_SE.sh rawdata/cryo_bgi clean/cryoinjury alligned/cryoinjury STAR_Indices/50genomeDirectory
qsub -hold_jid cryo_bgi -N cryo_nexseq scripts/STAR_SE.sh rawdata/cryo_nexseq clean/cryoinjury alligned/cryoinjury STAR_Indices/75genomeDirectory
qsub -hold_jid cryo_nexseq -N sham_bgi scripts/STAR_SE.sh rawdata/sham_bgi clean/sham/SE alligned/sham/SE STAR_Indices/50genomeDirectory
qsub -hold_jid sham_bgi -N sham_nexseq scripts/STAR_SE.sh rawdata/sham_nexseq clean/sham/SE alligned/sham/SE STAR_Indices/75genomeDirectory
qsub -hold_jid sham_nexseq -N sham_novaseq scripts/STAR_PE.sh rawdata/sham_novaseq clean/sham/PE alligned/sham/PE STAR_Indices/150genomeDirectory
qsub -hold_jid sham_novaseq -N uninjured_bgi scripts/STAR_SE.sh rawdata/uninjured_bgi clean/uninjured/SE alligned/uninjured/SE STAR_Indices/50genomeDirectory
qsub -hold_jid uninjured_bgi -N uninjured_genome scripts/STAR_SE.sh rawdata/uninjured_genome clean/uninjured/SE alligned/uninjured/SE STAR_Indices/50genomeDirectory
qsub -hold_jid uninjured_genome -N uninjured_hiseq scripts/STAR_PE.sh rawdata/uninjured_hiseq clean/uninjured/PE alligned/uninjured/PE STAR_Indices/150genomeDirectory