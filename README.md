# LocSpiral-LocBSharpen-LocBFactor-LocOccupancy
Local computational methods to improve the interpretability and analysis of cryo-EM maps.

In this work, we propose 1) approaches to enhance high-resolution features of cryo-EM maps, while preventing map distortions and 2) methods to obtain local B-factors and electron density occupancy maps. These algorithms have as common link the use of the spiral phase transformation and are called LocSpiral, LocBSharpen, LocBFactor and LocOccupancy. Currently the code runs in Matlab while a Python package is under development.

LocSpiral, LocBSharpen, LocBFactor and LocOccupancy are developed by the Vargas lab at the Optics Department of the Universidad Complutense de Madrid.

The Matlab code is in the Code folder. Please see the Examples folder to see two complete examples showing how to use the code. All the results of the examples can be downloaded from [2]. In addition, all the results shown in the paper [1] can be downloaded from [3]. Please if this work was useful for you cite the paper [1] (now a preprint). 

For bug reports, questions or comments please contact Javier Vargas (jvargas@fis.ucm.es)

[1] S. Kaur, J. Gomez-Blanco, A. Khalifa, S. Adinarayanan, R. Sanchez-Garcia, D. Wrapp, J. S. McLellan, K. H. Bui, J. Vargas, Local computational methods to improve the interpretability and analysis of cryo-EM maps, bioRxiv https://doi.org/10.1101/2020.05.11.088013

[2] https://mcgill-my.sharepoint.com/:u:/g/personal/javier_vargasbalbuena_mcgill_ca/EViU56bKLflDjg0N2gfQzAIBCzmWG7sMDjx1wPEtabWfuw?e=FDBG8B

[3] https://mcgill-my.sharepoint.com/:u:/g/personal/javier_vargasbalbuena_mcgill_ca/ERl72DDZ0T9Nr561_2QxaoABzS6ClE4gFdNPR5RvOCGRNw?e=cKlqpO

1. System requirements:
The code provided in the 'Code' folder requires Matlab to be run. We have tested the code on Matlab's versions 2018a, 2018b and 2019a running in Windows 10 and Linux (Ubuntu 16.04). This software requires Matlab´s Parallel Computing toolbox if parallel acceleration is desired. To run the methods proposed here it is not required any non-standard hardware.

2. Installation guide:
To use the proposed approaches copy the files of the 'Code folder to your local machine at your desired destination. No compilation is required so the install time on a "normal" desktop computer is zero, taking into account that Matlab is already installed. Examples of use of the proposed methods can be found in [2].
