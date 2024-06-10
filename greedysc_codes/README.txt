Before running the matlab codes, please read the following instructions.

1. Download the matlab codes for SSC, LRR, SCC, TSC from

SSC : http://www.cis.jhu.edu/~ehsan/Codes/SSC_ADMM_v1.1.zip
LRR : https://sites.google.com/site/guangcanliu/lrr%28motion_face%29.zip?attredirects=0
SCC : http://www.math.duke.edu/~glchen/scc.zip
TSC : http://www.nari.ee.ethz.ch/commth/research/downloads/supp_material_RSCT.zip


2. Unzip the files and place the folders in the current folder.

SSC_ADMM_v1.1, code2, TSC_supp_material_RSCT/include, SCC (make this folder and place the files there)


3. For motion segmentation and face clustering, you need to modify SSC_ADMM_v1.1/SSC.m so that it returns the estimated labels.

function [missrate,CMat] = SSC(X,r,affine,alpha,outlier,rho,s) to 
--> function [missrate,CMat,grps] = SSC(X,r,affine,alpha,outlier,rho,s)

And removed the line

missrate = Misclassification(grps,s);

since it doesn't need to be computed.


4. The following codes generates the figures and numbers in the paper.

run_sc_synthetic.m : Figures 1 and 2
run_sc_speedcomparison.m : Figure 3
run_sc_motionseg.m : Table 1
run_sc_faceimages.m : Table 2

--

If you have any questions or comments, please send an email to [dhpark at utexas dot edu]
