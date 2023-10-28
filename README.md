# AutoPepVax
An automatic epitope selection method for peptide vaccine design. The Python-based programs in this repository can produce a list of epitopes likely to induce tumor-specific T cells based on in-silico tools and a novel algorithm.

### Table of Contents
- [Dependencies](#Dependencies)
- [Installation](#Installation)
- [Quick Start Guide](#Quick-Start-Guide)
- [Output](#Output)
- [PF.py Function Descriptions](#PFpy-Function-Descriptions)

### Dependencies
The programs in this repository are easiest to use in the JupyterLab environment. The installation guide for JupyterLab can be found [here](https://jupyter.org/install).   
ChromeDriver is needed to use some of the web scraping functions. If you are using MacOS, it is recommended to download [ChromeDriver for Mac](https://formulae.brew.sh/cask/chromedriver) with [Homebrew](https://brew.sh/). ChromeDriver for Windows downloads can be found [here](https://chromedriver.chromium.org/downloads).  
The following Python dependencies must also be installed:
```bash
pip install pandas
pip install splinter
pip install selenium
pip install mhcflurry
pip install tensorflow
mhcflurry-downloads fetch
pip install scikit-learn
pip install mhcgnomes
```

### Installation
Download the repository at [https://github.com/enricobautista/AutoPepVax.git](https://github.com/enricobautista/AutoPepVax.git), or copy the following code into your command prompt at the desired directory.
```bash
git clone https://github.com/enricobautista/AutoPepVax.git
```
To access the code, run JupyterLab and navigate to the folder you downloaded which should named *AutoPepVax*. The folder may need to be unzipped before accessing it.

### Quick Start Guide
In JupyterLab, navigate to the *AutoPepVax* folder you downloaded, and open `Automate Protocol.ipynb`.  
Change the following variables in the first code block to specify the mutations and protein you would like to study:
|Variable|Description|Example|
|:---:|:---|:---|
|`cancer_name`|This should be the name of the cancer that harbors relevant mutations. A folder for the cancer_name will appear in the *Data* folder. This folder will contain any output files pertinent to that cancer's mutations.|`cancer_name = 'LUAD'`|
|`chromedriver_path`|This should be the file path to your ChromeDriver download.|`chromedriver_path = '/opt/homebrew/bin/chromedriver'`|
|`og_seq`|This will be a string representing the single-letter code for the primary sequence of the unmutated protein that is being studied.|`og_seq = 'MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLG'`|
|`mutations`|This will be a list of the missense mutations that are being studied.|`mutations = ['A10T,R23G']` |

Once the variables are correctly set, all the code blocks can be run.

### Output
Upon running the code in `Automate Protocol.ipynb`, a folder named after the `cancer_name` variable will appear in the *Data* folder. Inside this folder, the following five files will be produced:
|File|Description|
|:---|:---|
|CD4 Epitopes.csv|A list of all analyzed MHC II-restricted epitope-HLA alleles and their pertinent characteristics.|
|CD4 Filtered Epitopes.csv|A filtered list of MHC II-restricted epitope-HLA alleles that meet exclusion criteria.|
|CD8 Epitopes.csv|A list of all analyzed MHC I-restricted epitope-HLA alleles and their pertinent characteristics, including ID and score.|
|CD8 Filtered Epitopes.csv|A filtered list of MHC I-restricted epitope-HLA alleles that meet exclusion criteria.|

### PF.py Function Descriptions
The functions in `Automate Protocol.ipynb` come from `PF.py`. Here are descriptions of what each of these functions do.  
`get_NetMHCI`: Collects MHC I-restricted epitope-HLA complexes relevant to the listed mutations and set of HLA alleles and records their NetMHCpan-4.1 binding affinity in nM and rank.  
`get_NetMHCII`: Collects MHC II-restricted epitope-HLA complexes relevant to the listed mutations and set of HLA alleles and records their NetMHCIIpan-4.1 binding affinity in nM and rank.  
`get_NetMHC_stab`: Collects MHC II-restricted epitope-HLA complexes relevant to the listed mutations and set of HLA alleles and records their NetMHCstabpan stability prediction.  
`get_MHCFlurry`: Collects binding affinity in nM and rank data for the MHC I-restricted epitopes and the ratio of the wild-type epitope to mutant epitope rank for each HLA allele.    
`get_MHCI_IEDB_immunogenicity`: Collects the IEDB immunogenicity prediction values for MHC I-restricted epitopes.  
`get_MHCII_IEDB_immunogenicity`: Collects the IEDB immunogenicity prediction values for MHC II-restricted epitopes.  
`get_antigenicity`: Collects the VaxiJen antigenicity prediction values for all epitopes.  
`get_toxicity`: Collects the ToxinPred toxicity predictions for all epitopes.  
`get_prot_parameters`: Collects half-life, isoelectric point, instability index, aliphatic index, and GRAVY score for all epitopes.  
`get_hydropathicity`: Collects hydropathicity for all parameters.  
`get_allergenicity`: Collects the AllerTop 2.0 allergenicity predictions for all epitopes.  
`get_INFgamma`: Collects the INFepitope INFgamma prediction values for all epitopes.  
`rank_peptides`: Scores and classifies MHC I-restricted epitopes based on the likelihood that they are good candidates for a peptide vaccine. Epitope-HLA pairs with an ID of 1 were classified as likely candidates. Epitope-HLA pairs are ranked by score within their class.  
`parameter_filter`: Creates a new file that only includes peptides that are filtered by the user's choice of parameters, which may include IEDB immunogenicity, antigenicity, half-life, instability index, toxicity, allergenicity, INFgamma, and ID.  
