# ESS
An Unsupervised Approach for Eye Sclera Segmentation (ESS) 

--------------------------------------------------------------------------

This project implements an unsupervised sclera segmentation method for eye color images. The proposed approach operates on a visible spectrum RGB eye image and does not require any prior knowledge such as eyelid or iris center coordinate detection. The eye color input image is enhanced by an adaptive histogram normalization to produce a gray level image in which the sclera is highlighted. A feature extraction process is involved both in the image binarization and in the computation of scores to assign to each connected components of the foreground. The binarization process is based on clustering and adaptive thresholding. Finally, the selection of foreground components identifying the sclera is performed on the analysis of the computed scores and of the positions between the foreground components. The proposed method was ranked 2nd in the Sclera Segmentation and Eye Recognition Benchmarking Competition (SSRBC 2017), providing satisfactory performance in terms of precision. 

--------------------------------------------------------------------------

Matlab Implementation by Daniel Riccio, April 05, 2017. 
Copyright (C) 2017 Daniel Riccio (dnl.riccio@gmail.com)
