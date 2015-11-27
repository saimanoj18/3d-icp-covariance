Covariance of ICP with 3D Point to Point and Point to Plane Error Metrics

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
 * please report it to me so that I can remove these statements and clarify meself.
 *
 * My email is
 * saimanoj18@yahoo.com,
 * saimanoj18@gmail.com
 * saimanoj001@e.ntu.edu.sg
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


Please contact me at
saimanoj18@yahoo.com
saimanoj18@gmail.com
saimanoj001@e.ntu.edu.sg


