This repository contains the scripts used for the analyses of the article. The study aimed to identify founder variants in Saguenay-Lac-Saint-Jean population that
are likely to modify the protein, and could be pathogenic.


The main script of the repository is Filters_RFD_carriers_MAF_families_IBDsharing_AlphaMissense_pLIscore.sh which allows to filter variants according 
to frequency and likelihood of modifying the protein. 
The other scripts are used for specific criteria, refered to in the main script: 
ibd_select_one_haplotype_per_carrier.jl run with obd_matrix_haplotype.sh is used fo IBD sharing
identify_carriers.jl run with run_identify_carriers.sh is used for the families criterion. 
variant_annotation.R is used to annotate variants and their gene
pLI_score.R is used to add pLI score according to the gene on which the variant is located


Genealogies from BALSAC population register are available upon request at https:/www.balsac.caen/data-infrastructure/, 
following the Policy on Access to BALSAC Data for Research Purposes. 
Genotyping, WGS and imputations of CaG participants are available under restricted access from CaG biobank due to the informed consent 
given by study participants (https://cartagene.qc.ca/en/researchers/ access-request.html). 
Data from SLSJ familial asthma cohort are available in the Canadian Dataverse Repository Borealis on reasonable request 
(https:// doi. org/ 10. 5683/ SP3/ QEDOPE). 
