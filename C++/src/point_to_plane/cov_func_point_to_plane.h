/*************************************
 * The inputs for this function are
 * data_pi --> the corresponding points in the data_cloud or source cloud or Pi
 * model_qi --> the correspondences in the model or target or Qi
 * RT --> the transformation matrix as returned by ICP
 * Output is a 6x6 matrix which is the COVARIANCE of ICP in [x,y,z,a,b,c]
 * [a b c] are Yaw, Pitch and Roll respectively
***************************************/




void calculate_ICP_COV(pcl::PointCloud<pcl::PointNormal>& data_pi, pcl::PointCloud<pcl::PointNormal>& model_qi, Eigen::Matrix4f& transform, Eigen::MatrixXd& ICP_COV)
{

    double Tx = transform(0,3);
    double Ty = transform(1,3);
    double Tz = transform(2,3);
    double roll  = atan2f(transform(2,1), transform(2,2));
    double pitch = asinf(-transform(2,0));
    double yaw   = atan2f(transform(1,0), transform(0,0));

    double x, y, z, a, b, c;
    x = Tx; y = Ty; z = Tz;
    a = yaw; b = pitch; c = roll;// important // According to the rotation matrix I used and after verification, it is Yaw Pitch ROLL = [a,b,c]== [R] matrix used in the MatLab also :)

    /* Flushing out in the form of XYZ ABC */
    std::cout << "\nPrinting out [x, y, z, a, b, c] =  " <<x<<"    "<<y<<"    "<<z<<"    "<<a<<"    "<<b<<"    "<<c<<std::endl;

    //Matrix initialization
    Eigen::MatrixXd d2J_dX2(6,6);
    d2J_dX2 = Eigen::MatrixXd::Zero(6,6);

    /****  Calculating d2J_dX2  ****/
    for (size_t s = 0; s < data_pi.points.size(); ++s )
    {
        double pix = data_pi.points[s].x;
        double piy = data_pi.points[s].y;
        double piz = data_pi.points[s].z;
        double qix = model_qi.points[s].x;
        double qiy = model_qi.points[s].y;
        double qiz = model_qi.points[s].z;

        double nix = model_qi[s].normal_x;
        double niy = model_qi[s].normal_y;
        double niz = model_qi[s].normal_z;

        if (niz!=niz)// for nan removal in the input point cloud data:)
            continue;

        /***********************************************************************

d2J_dX2 -- X is the [R|T] in the form of [x, y, z, a, b, c]
x, y, z is the translation part 
a, b, c is the rotation part in Euler format
[x, y, z, a, b, c] is acquired from the Transformation Matrix returned by ICP.

Now d2J_dX2 is a 6x6 matrix of the form

d2J_dx2                   
d2J_dxdy    d2J_dy2
d2J_dxdz    d2J_dydz    d2J_dz2
d2J_dxda    d2J_dyda    d2J_dzda   d2J_da2
d2J_dxdb    d2J_dydb    d2J_dzdb   d2J_dadb   d2J_db2
d2J_dxdc    d2J_dydc    d2J_dzdc   d2J_dadc   d2J_dbdc   d2J_dc2

*************************************************************************/



        double 	d2J_dx2,     d2J_dydx,	  d2J_dzdx,   d2J_dadx,   d2J_dbdx,     d2J_dcdx,
                d2J_dxdy,    d2J_dy2,	  d2J_dzdy,   d2J_dady,   d2J_dbdy,     d2J_dcdy,
                d2J_dxdz,    d2J_dydz,    d2J_dz2,    d2J_dadz,   d2J_dbdz,     d2J_dcdz,
                d2J_dxda,    d2J_dyda,    d2J_dzda,   d2J_da2,	  d2J_dbda,     d2J_dcda,
                d2J_dxdb,    d2J_dydb,    d2J_dzdb,   d2J_dadb,   d2J_db2,      d2J_dcdb,
                d2J_dxdc,    d2J_dydc,    d2J_dzdc,   d2J_dadc,   d2J_dbdc,     d2J_dc2;

        // These terms are generated from the provided Matlab scipts. We just have to copy
        // the expressions from the matlab output with two very simple changes.
        // The first one being the the sqaure of a number 'a' is shown as a^2 in matlab,
        // which is converted to pow(a,2) in the below expressions.
        // The second change is to add ';' at the end of each expression :)
        // In this way, matlab can be used to generate these terms for various objective functions of ICP
        // and they can simply be copied to the C++ files and with appropriate changes to ICP estimation,
        // its covariance can be easily estimated.


        d2J_dx2 =

                2*pow(nix,2);


        d2J_dy2 =

                2*pow(niy,2);


        d2J_dz2 =

                2*pow(niz,2);


        d2J_dydx =

                2*nix*niy;


        d2J_dxdy =

                2*nix*niy;


        d2J_dzdx =

                2*nix*niz;


        d2J_dxdz =

                2*nix*niz;


        d2J_dydz =

                2*niy*niz;


        d2J_dzdy =

                2*niy*niz;


        d2J_da2 =

                (niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a))) - (2*nix*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) + 2*niy*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_db2 =

                (niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c))) - (2*niy*(pix*cos(b)*sin(a) + piz*cos(c)*sin(a)*sin(b) + piy*sin(a)*sin(b)*sin(c)) + 2*niz*(piz*cos(b)*cos(c) - pix*sin(b) + piy*cos(b)*sin(c)) + 2*nix*(pix*cos(a)*cos(b) + piz*cos(a)*cos(c)*sin(b) + piy*cos(a)*sin(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dc2 =

                (nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c))) - (2*niy*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b))) - 2*nix*(piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b))) + 2*niz*(piz*cos(b)*cos(c) + piy*cos(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));

        d2J_dxda =

                nix*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dadx =

                2*nix*(niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dyda =

                niy*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dady =

                2*niy*(niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dzda =

                niz*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dadz =

                2*niz*(niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dxdb =

                nix*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dbdx =

                2*nix*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dydb =

                niy*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dbdy =

                2*niy*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));

        d2J_dzdb =

                niz*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dbdz =

                2*niz*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dxdc =

                nix*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dcdx =

                2*nix*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dydc =

                niy*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));

        d2J_dcdy =

                2*niy*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dzdc =

                niz*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dcdz =

                2*niz*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dadb =

                (niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c))) - (2*nix*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niy*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dbda =

                (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c))) - (2*nix*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niy*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dbdc =

                (2*nix*(piy*cos(a)*cos(b)*cos(c) - piz*cos(a)*cos(b)*sin(c)) - 2*niz*(piy*cos(c)*sin(b) - piz*sin(b)*sin(c)) + 2*niy*(piy*cos(b)*cos(c)*sin(a) - piz*cos(b)*sin(a)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c))) + (niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dcdb =

                (2*nix*(piy*cos(a)*cos(b)*cos(c) - piz*cos(a)*cos(b)*sin(c)) - 2*niz*(piy*cos(c)*sin(b) - piz*sin(b)*sin(c)) + 2*niy*(piy*cos(b)*cos(c)*sin(a) - piz*cos(b)*sin(a)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c))) + (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dcda =

                (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c))) + (2*nix*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niy*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dadc =

                (niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c))) + (2*nix*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niy*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));





        Eigen::MatrixXd d2J_dX2_temp(6,6);

        d2J_dX2_temp << d2J_dx2,     d2J_dydx,	  d2J_dzdx,   d2J_dadx,   d2J_dbdx,     d2J_dcdx,
                        d2J_dxdy,    d2J_dy2,	  d2J_dzdy,   d2J_dady,   d2J_dbdy,     d2J_dcdy,
                        d2J_dxdz,    d2J_dydz,    d2J_dz2,    d2J_dadz,   d2J_dbdz,     d2J_dcdz,
                        d2J_dxda,    d2J_dyda,    d2J_dzda,   d2J_da2,	  d2J_dbda,     d2J_dcda,
                        d2J_dxdb,    d2J_dydb,    d2J_dzdb,   d2J_dadb,   d2J_db2,      d2J_dcdb,
                        d2J_dxdc,    d2J_dydc,    d2J_dzdc,   d2J_dadc,   d2J_dbdc,     d2J_dc2;


        d2J_dX2 = d2J_dX2 + d2J_dX2_temp;

    }// End of the FOR loop!!!

    std::cout << "\n**************\n Successfully Computed d2J_dX2 \n**************\n" << std::endl;



    // Now its time to calculate d2J_dZdX , where Z are the measurements Pi and Qi, X = [x,y,z,a,b,c]

    // n is the number of correspondences
    int n = data_pi.points.size();





    /*  Here we check if the number of correspondences between the source and the target point clouds are greater than 200.
      if yes, we only take the first 200 correspondences to calculate the covariance matrix
      You can try increasing it but if its too high, the system may run out of memory and give an exception saying

     terminate called after throwing an instance of 'std::bad_alloc'
     what():  std::bad_alloc
     Aborted (core dumped)

  */
    if (n > 200) n = 200;////////////****************************IMPORTANT CHANGE***** but may not affect********************/////////////////////////////////////////

    std::cout << "\nNumber of Correspondences used for ICP's covariance estimation = " << n << std::endl;

    Eigen::MatrixXd d2J_dZdX(6,6*n);

    for (int k = 0; k < n ; ++k) // row
    {
        //here the current correspondences are loaded into Pi and Qi
        double pix = data_pi.points[k].x;
        double piy = data_pi.points[k].y;
        double piz = data_pi.points[k].z;
        double qix = model_qi.points[k].x;
        double qiy = model_qi.points[k].y;
        double qiz = model_qi.points[k].z;

        double nix = model_qi[k].normal_x;
        double niy = model_qi[k].normal_y;
        double niz = model_qi[k].normal_z;

        if (niz!=niz) // for nan removal in input point cloud data
            continue;

        Eigen::MatrixXd d2J_dZdX_temp(6,6);


        double 	d2J_dpix_dx,    d2J_dpiy_dx,	d2J_dpiz_dx,  	   d2J_dqix_dx,    d2J_dqiy_dx,	   d2J_dqiz_dx,
                d2J_dpix_dy,    d2J_dpiy_dy,	d2J_dpiz_dy,   	   d2J_dqix_dy,    d2J_dqiy_dy,	   d2J_dqiz_dy,
                d2J_dpix_dz,    d2J_dpiy_dz,    d2J_dpiz_dz,       d2J_dqix_dz,    d2J_dqiy_dz,    d2J_dqiz_dz,
                d2J_dpix_da,    d2J_dpiy_da,    d2J_dpiz_da,       d2J_dqix_da,    d2J_dqiy_da,    d2J_dqiz_da,
                d2J_dpix_db,    d2J_dpiy_db,    d2J_dpiz_db,       d2J_dqix_db,    d2J_dqiy_db,    d2J_dqiz_db,
                d2J_dpix_dc,    d2J_dpiy_dc,    d2J_dpiz_dc,       d2J_dqix_dc,    d2J_dqiy_dc,    d2J_dqiz_dc;


        d2J_dpix_dx =

                2*nix*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));


        d2J_dpix_dy =

                2*niy*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));


        d2J_dpix_dz =

                2*niz*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));


        d2J_dpix_da =

                (2*niy*cos(a)*cos(b) - 2*nix*cos(b)*sin(a))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c))) + (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));


        d2J_dpix_db =

                (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a)) - (2*niz*cos(b) + 2*nix*cos(a)*sin(b) + 2*niy*sin(a)*sin(b))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dpix_dc =

                (2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));


        d2J_dpiy_dx =

                2*nix*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c));


        d2J_dpiy_dy =

                2*niy*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c));


        d2J_dpiy_dz =

                2*niz*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c));


        d2J_dpiy_da =

                (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c)) - (2*nix*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + 2*niy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dpiy_db =

                (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c)) + (2*nix*cos(a)*cos(b)*sin(c) - 2*niz*sin(b)*sin(c) + 2*niy*cos(b)*sin(a)*sin(c))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dpiy_dc =

                (2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c)) + (2*nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - 2*niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2*niz*cos(b)*cos(c))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dpiz_dx =

                2*nix*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c));


        d2J_dpiz_dy =

                2*niy*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c));


        d2J_dpiz_dz =

                2*niz*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c));


        d2J_dpiz_da =

                (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c)) + (2*nix*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2*niy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dpiz_db =

                (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c)) + (2*nix*cos(a)*cos(b)*cos(c) - 2*niz*cos(c)*sin(b) + 2*niy*cos(b)*cos(c)*sin(a))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dpiz_dc =

                (2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c)) - (2*niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - 2*nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + 2*niz*cos(b)*sin(c))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));


        d2J_dqix_dx =

                -2*pow(nix,2);


        d2J_dqix_dy =

                -2*nix*niy;


        d2J_dqix_dz =

                -2*nix*niz;


        d2J_dqix_da =

                -nix*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dqix_db =

                -nix*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dqix_dc =

                -nix*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dqiy_dx =

                -2*nix*niy;


        d2J_dqiy_dy =

                -2*pow(niy,2);


        d2J_dqiy_dz =

                -2*niy*niz;


        d2J_dqiy_da =

                -niy*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dqiy_db =

                -niy*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dqiy_dc =

                -niy*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));


        d2J_dqiz_dx =

                -2*nix*niz;


        d2J_dqiz_dy =

                -2*niy*niz;


        d2J_dqiz_dz =

                -2*pow(niz,2);


        d2J_dqiz_da =

                -niz*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));


        d2J_dqiz_db =

                -niz*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));


        d2J_dqiz_dc =

                -niz*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));



        d2J_dZdX_temp <<    d2J_dpix_dx,    d2J_dpiy_dx,	d2J_dpiz_dx,  	   d2J_dqix_dx,    d2J_dqiy_dx,	   d2J_dqiz_dx,
                            d2J_dpix_dy,    d2J_dpiy_dy,	d2J_dpiz_dy,   	   d2J_dqix_dy,    d2J_dqiy_dy,	   d2J_dqiz_dy,
                            d2J_dpix_dz,    d2J_dpiy_dz,        d2J_dpiz_dz,       d2J_dqix_dz,    d2J_dqiy_dz,    d2J_dqiz_dz,
                            d2J_dpix_da,    d2J_dpiy_da,        d2J_dpiz_da,       d2J_dqix_da,    d2J_dqiy_da,    d2J_dqiz_da,
                            d2J_dpix_db,    d2J_dpiy_db,        d2J_dpiz_db,       d2J_dqix_db,    d2J_dqiy_db,    d2J_dqiz_db,
                            d2J_dpix_dc,    d2J_dpiy_dc,        d2J_dpiz_dc,       d2J_dqix_dc,    d2J_dqiy_dc,    d2J_dqiz_dc;


        d2J_dZdX.block<6,6>(0,6*k) = d2J_dZdX_temp;

    }


    //By reaching here both the matrices d2J_dX2 and d2J_dZdX are calculated and lets print those values out;

    //std::cout << "\n Finally here are the two matrices \n\n" << "d2J_dX2 = \n " << d2J_dX2 <<std::endl;
    //std::cout << "\n\n\n" << "d2J_dZdX = \n " << d2J_dZdX <<std::endl;




    /**************************************
     *
     * Here we create the matrix cov(z) as mentioned in Section 3.3 in the paper, "Covariance of ICP with 3D Point to Point and Point to Plane Error Metrics"
     *
     * ************************************/
    Eigen::MatrixXd cov_z(6*n,6*n);
    cov_z = 0.01 * Eigen::MatrixXd::Identity(6*n,6*n);



    ICP_COV =  d2J_dX2.inverse() * d2J_dZdX * cov_z * d2J_dZdX.transpose() * d2J_dX2.inverse();


    std::cout << "\n\n********************** \n\n" << "ICP_COV = \n" << ICP_COV <<"\n*******************\n\n"<< std::endl;

    std::cout << "\nSuccessfully Computed the ICP's Covariance !!!\n" << std::endl;

}














