# A-clustering-methodology-for-European-banks-business-models
This repository contains the MATLAB functions relative to the paper 'Business models of the Banks in the Euro Area', No. 2070, ECB Working Paper (https://www.ecb.europa.eu/pub/pdf/scpwps/ecb.wp2070.en.pdf?ee58f8028aa3d7b55dd977292218b268). 

In that paper, the trimmed version of the factorial k-means algorithm by Vichi and Kiers (https://s3.amazonaws.com/academia.edu.documents/46752397/Factorial_k-means_analysis_for_two-way_d20160624-9057-19q8f9o.pdf?AWSAccessKeyId=AKIAIWOWYYGZ2Y53UL3A&Expires=1514563179&Signature=K8ahFxipvNcKrT2RnoD8s5Ireqw%3D&response-content-disposition=inline%3B%20filename%3DFactorial_k-means_analysis_for_two-way_d.pdf) is proposed.

The function 'tfkm.m' performs a version of trimmed factorial k-means where the latent rank and the number of clusters are selected maximizing Hartigan's statistics, while the trimming proportion must be set by the user.

The function 'tfkm_alpha.m' performs a version of trimmed factorial k-means where the trimming proportion is also selected maximizing Hartigan's statistics. 

Both MATLAB functions contain the detailed explanation of input and output arguments.
The function 'tfkm_data.m' provides some examples about the use of both functions.
