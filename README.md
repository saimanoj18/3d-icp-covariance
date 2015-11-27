Covariance of ICP with 3D Point to Point and Point to Plane Error Metrics  

site: https://sites.google.com/site/icpcovariance/  

/***************************************  
 * This file contains a class that has functions to create SHOT descriptors for feature matching  
 * Feature matching is used to bring the two point clouds close so that ICP can converge to a global minima  
 *  
 *  
 *  
 * Generally it is observed that the ICP covaraince is lower at global minima than at a local minima.  
 * And my initial experimental results show that  
 * the ICP's covariance is lower when it converges to a global minima as compared to the ICP's covariance when it converges to a local minima ( when ICP fails).  
 *  
 * If by any chance, you come across a case where ICP converges to a global minima and yet offers higher covariance than the covariance of local minima,  
 * please report it to me so that I can remove these statements and clarify myself.  
 *
 * My email is  
 * saimanoj18@yahoo.com,  
 * saimanoj001@e.ntu.edu.sg  
 * https://sites.google.com/site/saimanoj18/  
 *  
 *  
 * **************************************/  

There are two parts in it. The first one is MatLab scripts while the second one is C++  

You need to work with MatLab scripts if you want to change the objective function for ICP. We provide for both 3D point to point and 3D point to plane error metrics  


Required Packages:  

PCL (Point Cloud Library)  

To RUN :  

(cd to C++/build folder)

cd code/C++/build

cmake ..

make

./cov_point_to_point

./cov_point_to_plane


Sample Point Clouds with the package  

data_bunny.pcd  

model_bunny.pcd  

// Important things to Note  
There are comments inline of the code!  

Finally...if this work is useful for you...Please cite :)  



@INPROCEEDINGS{3d_icp_cov,   
author={Prakhya, S.M. and Liu Bingbing and Yan Rui and Weisi Lin},  
booktitle={Machine Vision Applications (MVA), 2015 14th IAPR International Conference on},  
title={A closed-form estimate of 3D ICP covariance},  
year={2015},  
pages={526-529},  
doi={10.1109/MVA.2015.7153246},  
month={May},}  


Please contact me at  
saimanoj18@yahoo.com  
saimanoj001@e.ntu.edu.sg  


