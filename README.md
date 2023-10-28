# AutoPepVax
An automatic epitope selection method for peptide vaccine design.

PF.py Function Descriptions  
`get_NetMHCI`: Collects MHC I-restricted epitope-HLA complexes relevant to the listed mutations and set of HLA alleles and records their NetMHCpan-4.1 binding affinity in nM and rank.  
`get_NetMHCII`: Collects MHC II-restricted epitope-HLA complexes relevant to the listed mutations and set of HLA alleles and records their NetMHCIIpan-4.1 binding affinity in nM and rank.  
`get_NetMHC_stab`: Collects MHC II-restricted epitope-HLA complexes relevant to the listed mutations and set of HLA alleles and records their NetMHCstabpan stability prediction.  
`get_MHCFlurry`: Collects binding affinity in nM and rank data for the MHC I-restricted epitopes and the ratio of the wild-type epitope to mutant epitope rank for each HLA allele.  
`get_DeepImmuno`: Collects the DeepImmuno immunogenicity prediction values for MHC I-restricted epitopes.  
`get_MHCI_IEDB_immunogenicity`: Collects the IEDB immunogenicity prediction values for MHC I-restricted epitopes.  
`get_MHCII_IEDB_immunogenicity`: Collects the IEDB immunogenicity prediction values for MHC II-restricted epitopes.  
`get_antigenicity`: Collects the VaxiJen antigenicity prediction values for all epitopes.  
`get_toxicity`: Collects the ToxinPred toxicity predictions for all epitopes.  
`get_prot_parameters`: Collects half-life, isoelectric point, instability index, aliphatic index, and GRAVY score for all epitopes.  
`get_hydropathicity`: Collects hydropathicity for all parameters.  
`get_allergenicity`: Collects the AllerTop 2.0 allergenicity predictions for all epitopes.  
`get_INFgamma`: Collects the INFepitope INFgamma prediction values for all epitopes.  
`rank_peptides`: Scores and classifies the epitopes based on the likelihood that they are good candidates for a peptide vaccine. Epitope-HLA pairs with an ID of 1 were classified as likely candidates. Epitope-HLA pairs are ranked by score within their class.  
`parameter_filter`: Creates a new file that only includes peptides that are filtered by the user's choice of parameters, which may include IEDB immunogenicity, antigenicity, half-life, instability index, toxicity, allergenicity, INFgamma, and ID.  
