# Specificity in Calcium Channels | A Comparative Study of CaV1, CaV2, and CaV3
The project aimed to investigate the phylogeny and specificity of calcium channel (CaV) subtypes CaV1, CaV2, and CaV3. We compiled a list of CaV homologs, performed multiple sequence alignment, and constructed phylogenetic trees to analyze the relationships between the different CaV subfamilies. The results demonstrate a clear evolutionary separation between the three main CaV subtypes. Furthermore, by comparing the homologous sequences, we identified specific amino acid positions that differ between the CaV subfamilies, assuming these positions might explain the biophysical differences between the channel types. 

##### find_homologs_swissprot.py & find_insecta_homologs.py
These scripts contain basic functions for reading and writing FASTA files, performing MSAs, and other sequence-processing tasks. Additionally, they include functions to filter homologous sequences based on annotation, sequence length, and other criteria.

##### compare_seq_by_type.py & sequence_specificity_analysis.py
These scripts analyze aligned positions to identify subtype-specific residues. They classify the observed mutations based on amino acid families and identify segments where several substitution mutations occur consecutively.
