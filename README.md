# Natural variation in the sequestosome-related gene, *sqst-5*, underlies zinc homeostasis in *Caenorhabditis elegans*
#### *Data, scripts, and figures for zinc *sqst-5* manuscript*

Zinc is an essential trace element that acts as a cofactor for many enzymes and transcription factors required for cellular growth and development. Altering intracellular zinc levels can produce dramatic effects ranging from cell proliferation to cell death. To avoid such fates, cells have evolved mechanisms to handle both an excess and a deficiency of zinc. Zinc homeostasis is largely maintained via zinc transporters, permeable channels, and other zinc-binding proteins. Variation in these proteins might affect their ability to interact with zinc, leading to either increased sensitivity or resistance to natural zinc fluctuations in the environment. We can leverage the power of the roundworm nematode *Caenorhabditis elegans* as a tractable metazoan model for quantitative genetics to identify genes that could underlie variation in responses to zinc. We found that the laboratory-adapted strain (N2) is resistant and a natural isolate from Hawaii (CB4856) is sensitive to micromolar amounts of exogenous zinc supplementation. Using a panel of recombinant inbred lines, we identified two large-effect quantitative trait loci (QTL) on the left arm of chromosome III and the center of chromosome V that are associated with zinc responses. We validated and refined both QTL using near-isogenic lines (NILs) and identified a naturally occurring deletion in *sqst-5*, a sequestosome-related gene, that is associated with resistance to high exogenous zinc. We found that this deletion is relatively common across strains within the species and that variation in *sqst-5* is associated with zinc resistance. Our results offer a possible mechanism for how organisms can respond to naturally high levels of zinc in the environment and how zinc homeostasis varies among individuals.

### Manuscript
Link: [Here!](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008986)

### Usage
1. Download git repo
2. Install necessary packages (*tidyverse, linkagemapping, cegwas2, ggtree, phangorn, ape, easysorter, ggplotify*)
3. Set working directory at top of script (!!!!)
4. Run script

### Special Package Installs
- `linkagemapping` can be found [here](https://github.com/AndersenLab/linkagemapping) and installed with `devtools::install_github("AndersenLab/linkagemapping")`
- `cegwas2` can be found [here](https://github.com/AndersenLab/cegwas2) and installed with `devtools::install_github("AndersenLab/cegwas2")`
- `easysorter` can be found [here](https://github.com/AndersenLab/easysorter) and installed with `devtools::install_github("AndersenLab/easysorter")`

### Shiny
As part of this manuscript, I developed and published an R Shiny web application that can be run to analyze all the NIL phenotype data from this manuscript. A link to its github page can be found [here](https://github.com/katiesevans/finemap_NIL) and a link to the application can be found [here](https://katiesevans9.shinyapps.io/QTL_NIL/)!

### Files
- **S1 File. Dose response phenotype data.** Processed phenotype data from zinc dose response (standard HTA)
- **S2 File. Zinc response heritability.** Phenotypic values and used to calculate heritability and calculated heritabilities for all four zinc response traits (standard HTA)
- **S3 File. RIAIL phenotype data.** Phenotypic values for all 121 set 1 RIAILs, 253 set 2 RIAILs, and parent strains (N2 and CB4856) in response to zinc (standard HTA) 
- **S4 File. Linkage mapping results.** Linkage mapping LOD scores at 13,003 genomic markers for all four zinc-response traits with the set 2 RIAILs 
- **S5 File. Summary of two-dimensional genome scan.** Summary of the scan2 object containing data from the two-dimensional genome scan with animal optical density (median.EXT) in zinc 
- **S6 File. List of zinc-related genes.** List of all previously known zinc genes (zinc transporters and hits from mutant screens), their location in the genome, and if they have variation in CB4856 
- **S7 File. NIL sequence data.** VCF from the whole-genome sequencing for all the NILs in this study 
- **S8 File. NIL genotype data.** Simplified genotypes of the NILs in the study 
- **S9 File. NIL phenotype data.** Raw pruned phenotypes for the NILs on chromosomes III, IV, V, and X (standard HTA). 
- **S10 File. Statistical significance for NIL and CRISPR assays.** Pairwise statistical significance for all strains and high-throughput assays 
- **S11 File. ChrV NIL breakup phenotype data.** Raw pruned phenotypes for the NILs used to break up the QTL interval on chromosome V (standard HTA) 
- **S12 File. Modified HTA dose response phenotype data.** Raw pruned phenotypes for the parental dose response with the modified HTA  
- **S13 File. ChrIII dominance assay phenotype data.** Raw pruned phenotypes for the chromosome III dominance assay (modified HTA) 
- **S14 File. Genes in the chrIII QTL.** List of all genes in the chromosome III interval, their functional descriptions and GO annotations, and if they have variation in CB4856 
- **S15 File. Expression QTL mapping results for *sqst-5*.** Linkage mapping results for the *sqst-5* expression data obtained with the set 1 RIAILs 
- **S16 File. Mediation estimates for chrIII QTL.** Mediation estimates for the chromosome III QTL, including *sqst-5* 
- **S17 File. N2 and CB4856 *sqst-5* deletion phenotype data.** Raw pruned phenotypes for the *sqst-5* deletions in the parental backgrounds (modified HTA) 
- **S18 File. NIL *sqst-5* deletion phenotype data.** Raw pruned phenotypes for the *sqst-5* deletions in the NIL (ECA859) background (modified HTA) 
- **S19 File. Reciprocal hemizygosity for *sqst-5* phenotype data.** Raw pruned phenotypes for the reciprocal hemizygosity assay using the *sqst-5* deletions in the NIL background (modified HTA) 
- **S20 File. Sequence of *sqst-5* in N2 and CB4856.** Raw sequence of *sqst-5* for both N2 and CB4856 
- **S21 File. *sqst-5* gene alignment.** Nucleotide alignment of *sqst-5* for N2 and CB4856 
- **S22 File. SQST-5 protein alignment.** Protein alignment of SQST-5 for N2 and CB4856 
- **S23 File. Structural variation in *sqst-5* in 328 wild isolates.** List of all 328 wild isolates, their isolation location, and extent of structural variation in *sqst-5* 
- **S24 File. *sqst-5* phylogenetic tree.** Neighbor-joining tree for all 328 wild isolates using variants near *sqst-5*
- **S25 File. Tajima’s D for chrIII QTL.** Tajima’s D calculated for a sliding window across the zinc-response interval on chromosome III.
- **S26 File. Wild isolate phenotype data.** Residual phenotypic values for all 81 wild isolates in response to zinc (standard HTA) 
- **S27 File. GWA mapping results.** GWA mapping significance values for all markers across the genome for all four zinc-response traits
- **S28 File. Heavy metals linkage mapping results.** Linkage mapping LOD scores at 13,003 genomic markers for all four metal-response traits with the set 2 RIAILs for four heavy metals

