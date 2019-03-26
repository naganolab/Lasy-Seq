# The scripts for the analysis on "Lasy-Seq"

# Author
Mari Kamitani, Makoto Kashima, Ayumi Tezuka and Atsushi J. Nagano.

# Dependencies
R version 3.5.0
gridExtra version 2.3
TCC version 1.16.0
ggplot2 version 3.1.0

# Description
The scripts for analysis and preparation of the figures of the below paper.
The attribute of samples, gene description and the expression data (rpm) required for the analysis are included.

  TempResponse_A_thaliana: R scripts for analysis for temperature response of A. thaliana.

  ERCC_O_sativa: R scripts for analysis on false-assignment rate in Lasy-Seq.

  script for RSEMout: Python scripts for preparing "rawcnt" and "rpm" files from RSEM data.

# Usage
Download https://github.com/naganolab/Lasy-Seq. 
For the analysis on temperature response of A.thaliana, run R, execute script from no.1 to no.8.
For the analysis on false assignment rate, run R, execute script "190323_PrepRPM.R" and　calculate the false-assignment rate following to the method shown in Supplementary figure 2.

Comparison-conventinal-method_Lasy-Seq contains scripts and data for comparison of RNA-Seq results between a conventinal method and Lasy-Seq. For the analysis, run R and execute script "analysis.R"

# Citation
Mari Kamitani, Makoto Kashima, Ayumi Tezuka and Atsushi J. Nagano, 2018, Lasy-Seq: a high-throughput library preparation method for RNA-Seq and its application in the analysis of plant responses to fluctuating temperatures, bioRxiv, doi: https://doi.org/10.1101/463596

