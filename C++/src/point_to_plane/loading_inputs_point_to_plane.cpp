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

#include "cbshot.h"


int main(int argc, char** argv)
{
    // Two Input Point Clouds
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out (new pcl::PointCloud<pcl::PointXYZ>);


    /************************************/
    if (pcl::io::loadPCDFile<pcl::PointXYZ> ("../data_bunny.pcd", *cloud_in) == -1) //* load the file
    {
        PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
        return (-1);
    }

    if (pcl::io::loadPCDFile<pcl::PointXYZ> ("../model_bunny.pcd", *cloud_out) == -1) //* load the file
    {
        PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
        return (-1);
    }
    /**************************************/


     std::cout << "Status : Loaded both the point clouds and started to compute SHOT descriptors" << endl;



    /***************************************************/
    // Model's Normals should be calculated..checked via (RPi-Qi + T).Ni
    // Here Ni is the Normal of Qi
    //hence we need to calculate the normal of model

    //    Pi is data -> cloud_in
    //    Qi is model -> cloud_out
    /***************************************************/




     /******************************************************
      *
      *
      *
      * Below, we bring the two input point clouds close to each other, as ICP can converge to a global minimum
      * We do it by matching keypoints by feature descriptors, then estimating the consensual transformation via RANSAC
      * Then transform the point clouds and pass them for ICP refinement
      *
      *
      *
      * ****************************************************/




    cbshot cb;

    cb.cloud1 = *cloud_in;
    cb.cloud2 = *cloud_out;

    cb.calculate_voxel_grid_keypoints(0.05);
    cb.calculate_normals(0.02);
    cb.calculate_SHOT(0.30);



    pcl::Correspondences corresp;
    pcl::registration::CorrespondenceEstimation<pcl::SHOT352, pcl::SHOT352> shot_corr;
    shot_corr.setInputSource(cb.cloud1_shot.makeShared());
    shot_corr.setInputTarget(cb.cloud2_shot.makeShared());
    std::cout << "Status : Calculated SHOT descriptors and finding the correspondences" << endl;
    cout << "May take a while...depends on the number of feature descriptors and its support size" << endl;
    shot_corr.determineReciprocalCorrespondences(corresp);


    pcl::CorrespondencesConstPtr correspond = boost::make_shared< pcl::Correspondences >(corresp);

    pcl::Correspondences corr;
    pcl::registration::CorrespondenceRejectorSampleConsensus< pcl::PointXYZ > Ransac_based_Rejection;
    Ransac_based_Rejection.setInputSource(cb.cloud1_keypoints.makeShared());
    Ransac_based_Rejection.setInputTarget(cb.cloud2_keypoints.makeShared());
    Ransac_based_Rejection.setInputCorrespondences(correspond);
    Ransac_based_Rejection.getCorrespondences(corr);

    Eigen::Matrix4f ransac_mat = Ransac_based_Rejection.getBestTransformation();
    cout << "RANSAC based Transformation Matrix : \n" << ransac_mat << endl;

    pcl::transformPointCloud(*cloud_in, *cloud_in, ransac_mat);


    std::cout << "Now a window with the established correspondences between source and target keypoints pops up!" << endl;
    cout << "Based on these correspondences, the source and target point clouds are moved closer and then ICP works on them" << endl;
    cout << "If the established correspondences are not good, then the estimated transformation is not reliable!" << endl;
    cout << "It means that there is high probability that ICP may converge to a local minima !!!" << endl;



    /******************************************
     * If the ICP covariance values do not make sense
     * or
     * the model and data are not registered properly
     * then
     * look at the correspondences established by feature matching after RANSAC
     * It will show if the transformation can make two input point clouds close to each other
     *
     * Please remember that ICP works well for point clouds that are close enough
     * or else
     * ICP may converge to local minama..other way of saying it is that ICP fails
     * ***************************************/



    /************************************************************************************/

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer1 (new pcl::visualization::PCLVisualizer ("Feature Correspondences after RANSAC"));
    viewer1->setBackgroundColor (0, 0, 0);

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color3(cb.cloud1_keypoints.makeShared(), 255, 0, 0);
    viewer1->addPointCloud<pcl::PointXYZ> (cb.cloud1_keypoints.makeShared(), single_color3, "sample cloud3");
    viewer1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "sample cloud3");
    viewer1->initCameraParameters ();

    Eigen::Matrix4f mat1; mat1 << 1,0,0,3,0,1,0,0,0,0,1,0,0,0,0,1;
    pcl::PointCloud<pcl::PointXYZ> temp_here;
    pcl::transformPointCloud(cb.cloud2_keypoints,temp_here,mat1);

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color4(temp_here.makeShared(), 0, 255, 0);
    viewer1->addPointCloud<pcl::PointXYZ> (temp_here.makeShared(), single_color4, "sample cloud4");
    viewer1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "sample cloud4");

    viewer1->addCorrespondences<pcl::PointXYZ>(cb.cloud1_keypoints.makeShared(), temp_here.makeShared(), corr, "correspondences");

    /**************************************************************************************/



    cout << "Started Point to Plane ICP" << endl;

    // Now do the Point to Plane ICP

    pcl::PointCloud<pcl::PointNormal> cld1, cld2;
    cbshot temp;
    temp.cloud1 = *cloud_in;
    temp.cloud2 = *cloud_out;
    temp.calculate_normals(0.02);
    pcl::copyPointCloud(temp.cloud1, cld1);// Take care that you pass the transformed clouds after the feature matching and RANSAC and not the input ones as they may not be close for global convergence
    pcl::copyPointCloud(temp.cloud2,cld2);
    pcl::copyPointCloud(temp.cloud1_normals, cld1);
    pcl::copyPointCloud(temp.cloud2_normals, cld2);

    Eigen::Matrix4f final_transformation;
    final_transformation << 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1;

    pcl::Correspondences correspondeces_reciprocal_shot;

    int iterations = 50;// We just hard code the number of ICP iterations to be 50
    for (int i = 0; i < iterations; i++)
    {

        pcl::registration::CorrespondenceEstimation<pcl::PointNormal, pcl::PointNormal> corr_est;
        corr_est.setInputSource(cld1.makeShared()); // + setIndices(...)
        corr_est.setInputTarget(cld2.makeShared());
        corr_est.determineReciprocalCorrespondences(correspondeces_reciprocal_shot);
        //cout << "No. of Reciprocal Correspondences : " << correspondeces_reciprocal_shot.size() << endl;


        Eigen::Matrix4f transform_eigen;
        pcl::registration::TransformationEstimationPointToPlaneLLS<pcl::PointNormal, pcl::PointNormal> pp_lls;
        pp_lls.estimateRigidTransformation (cld1, cld2, correspondeces_reciprocal_shot, transform_eigen);


        // rotate/transform data based on this estimated transformation
        pcl::transformPointCloud(cld1, cld1, transform_eigen);


        // accumulate incremental tf
        final_transformation = transform_eigen * final_transformation;


    }


    cout <<  "Final transformation (feature matching * ICP alone)  :\n" << ransac_mat * final_transformation << endl;

    cout << "*************\n!!! Here we estimate the ICP's covariance of transformation between closely moved point clouds !!! \n*************"<< endl;

    //cout << "No. of Reciprocal Correspondences : " << correspondeces_reciprocal_shot.size() << endl;


    std::vector<int> data_idx;
    std::vector<int> model_idx;

    pcl::Correspondence temp1;
    for (int i = 0; i < correspondeces_reciprocal_shot.size(); i++)
    {
        temp1 = correspondeces_reciprocal_shot[i];
        data_idx.push_back(temp1.index_query);
        model_idx.push_back(temp1.index_match);
    }


    pcl::PointCloud<pcl::PointNormal> data_pi; //Put all the pi in this cloud and its size will be equal to number of correspondences
    pcl::PointCloud<pcl::PointNormal> model_qi;// Put all the qi in this cloud and its size will be equal to number of correspondences

    pcl::copyPointCloud(cld1,data_idx,data_pi);
    pcl::copyPointCloud(cld2,model_idx,model_qi);

    Eigen::MatrixXd ICP_COV(6,6);
    ICP_COV = Eigen::MatrixXd::Zero(6,6);


    clock_t start, end;
    double cpu_time_used;
    start = clock();

    calculate_ICP_COV(data_pi, model_qi, final_transformation, ICP_COV);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "\nComputational Time for ICP's covariance computation alone (in sec) = " << cpu_time_used;

    std::cout << "\n\n" << "ICP_COV Trace = \n\n" << ICP_COV.trace() <<"\n\n\n"<< std::endl;






    /***************************************************/

    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("Registered Point Clouds after ICP"));
    viewer->setBackgroundColor (0, 0, 0);

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointNormal> single_color1(cld1.makeShared(), 255, 0, 0);
    viewer->addPointCloud<pcl::PointNormal> (cld1.makeShared(), single_color1, "sample cloud1");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "sample cloud1");
    viewer->initCameraParameters ();

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointNormal> single_color2(cld2.makeShared(), 0, 255, 0);
    viewer->addPointCloud<pcl::PointNormal> (cld2.makeShared(), single_color2, "sample cloud2");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "sample cloud2");

    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }


    /***************************************************/



    return 0;
}

