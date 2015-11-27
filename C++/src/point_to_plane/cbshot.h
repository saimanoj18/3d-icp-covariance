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


#include <cmath>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/transforms.h>
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/correspondence_estimation.h>
#include <pcl/registration/correspondence_rejection_sample_consensus.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/kdtree/flann.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/impl/kdtree_flann.hpp>
#include <pcl/registration/transformation_estimation_svd.h>
#include <pcl/registration/transformation_estimation_lm.h>
#include <pcl/registration/transformation_estimation_point_to_plane.h>
#include <pcl/registration/transformation_estimation_point_to_plane_lls.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/shot_omp.h>
#include <pcl/features/fpfh_omp.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/registration/ia_ransac.h>
#include <pcl/registration/icp.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d_omp.h>

#include "cov_func_point_to_plane.h"



class cbshot
{

public :

    pcl::PointCloud<pcl::PointXYZ> cloud1, cloud2;

    pcl::PointCloud<pcl::Normal> cloud1_normals, cloud2_normals;
    pcl::PointCloud<pcl::PointXYZ> cloud1_keypoints, cloud2_keypoints;

    pcl::PointCloud<pcl::SHOT352> cloud1_shot, cloud2_shot;
    pcl::PointCloud<pcl::FPFHSignature33> cloud1_fpfh, cloud2_fpfh;



    
    void calculate_normals ( float radius )
    {
        // Estimate the normals.
        pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> normalEstimation;
        normalEstimation.setRadiusSearch(radius);
        normalEstimation.setNumberOfThreads(12);
        pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
        normalEstimation.setSearchMethod(kdtree);

        normalEstimation.setInputCloud(cloud1.makeShared());
        normalEstimation.compute(cloud1_normals);

        normalEstimation.setInputCloud(cloud2.makeShared());
        normalEstimation.compute(cloud2_normals);
    }




    void  calculate_voxel_grid_keypoints ( float leaf_size )
    {
        // Find Keypoints on the input cloud
        pcl::VoxelGrid<pcl::PointXYZ> voxel_grid;
        voxel_grid.setLeafSize(leaf_size, leaf_size, leaf_size);

        voxel_grid.setInputCloud(cloud1.makeShared());
        voxel_grid.filter(cloud1_keypoints);

        voxel_grid.setInputCloud(cloud2.makeShared());
        voxel_grid.filter(cloud2_keypoints);
    }





    void calculate_SHOT ( float radius )
    {

        // SHOT estimation object.
        pcl::SHOTEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::SHOT352> shot;
        shot.setRadiusSearch(radius);
        shot.setNumberOfThreads(12);
        pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
        shot.setSearchMethod(kdtree);

        shot.setInputCloud(cloud1_keypoints.makeShared());
        shot.setSearchSurface(cloud1.makeShared());
        shot.setInputNormals(cloud1_normals.makeShared());
        shot.compute(cloud1_shot);

        shot.setInputCloud(cloud2_keypoints.makeShared());
        shot.setSearchSurface(cloud2.makeShared());
        shot.setInputNormals(cloud2_normals.makeShared());
        shot.compute(cloud2_shot);

    }




    void calculate_FPFH ( float radius )
    {

        // SHOT estimation object.
        pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh;
        fpfh.setRadiusSearch(radius);
        fpfh.setNumberOfThreads(12);
        pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
        fpfh.setSearchMethod(kdtree);

        fpfh.setInputCloud(cloud1_keypoints.makeShared());
        fpfh.setSearchSurface(cloud1.makeShared());
        fpfh.setInputNormals(cloud1_normals.makeShared());
        fpfh.compute(cloud1_fpfh);


        fpfh.setInputCloud(cloud2_keypoints.makeShared());
        fpfh.setSearchSurface(cloud2.makeShared());
        fpfh.setInputNormals(cloud2_normals.makeShared());
        fpfh.setNumberOfThreads(12);
        fpfh.compute(cloud2_fpfh);


    }


};


