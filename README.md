# ParallaxBA

This is the MATLAB code comparing Standard Bundle Adjustment(SBA) and Parallax Bundle Adjustment(PBA) for our paper **"Comparison Between MATLAB Bundle Adjustment Function and Parallax Bundle Adjustment".**

SBA MATLAB function from the Computer Vision toolbox is used, and PBA MATLAB is implemented by ourselves. 
For the comparison between the two BA methods, the "MALAGA PARKING-6L" [1] and "Starry Night" [2] datasets are used.

To compare two BA with "MALAGA PARKING-6L" and "Starry Night" datasets in MATLAB, **'Malaga_BundleAdjustment.m'** and **'Starry_BundleAdjustment.m'** can be executed, respectively.

Visual-Inertial SLAM can be executed with **"BAwIMU.m"**, which use "Starry Night" dataset.   

In the case of the "Starry Might" dataset, the type of data can be changed between the data with 40, 80, 100, and 500 features. The type of initial value and the number of images also can be chosen to run MATLAB with different initial inputs. 
   
                
                                
----         
### References
1. L. Zhao, S. Huang, L. Yan, J. J. Wang, G. Hu, and G. Dissanayake, “Large-scale monocular SLAM by local bundle adjustment and map joining”, in 2010 11th International Conference on Control Automation Robotics Vision, 7-10 Dec. 2010 2010, pp. 431-436, doi: 10.1109/ICARCV.2010.5707820.
2. L. E. Clement, V. Peretroukhin, J. Lambert, and J. Kelly, “The Battle for Filter Supremacy: A Comparative Study of the Multi-State Constraint Kalman Filter and the Sliding Window Filter”, 2015 12th Conference on Computer and Robot Vision, 3-5 June 2015, pp. 23-30, doi: 10.1109/CRV.2015.11.               
----
