VIPSS - Zhiyang Huang (2019)

------------------------------------

This code implements the algorithm described in

  **Variational Implicit Point Set Surface**  
   Zhiyang Huang, Nathan Carr, Tao Ju.  
   *ACM Transactions on Graphics (Proc. ACM Siggraph 2019)*  

The primary function of this program is to predict the normal as well as the underlying surface of a given set of unoriented points.

Currently, the code is only tested on Mac OS X 10.13 or above, it should work at other platforms with minor modification.


BUILDING
======================================================================================================


The code has only two dependencies: 1)Armadillo,   2)NLOPT

1) http://arma.sourceforge.net/

2) https://nlopt.readthedocs.io/en/latest/

You can download & install them by yourself, or run the env.sh script which will install homebrew first, and then install the two dependencies via homebrew.

$source env.sh  

Then go to the vipss folder, build the Cmake file and make:  
$cd vipss  
$cmake .  
$make  

In the vipss directory, there should be an executable called "vipss" (or "vipss.exe" on Windows if it is successfully built).


RUNNING
======================================================================================================

To run the code from the command line, type:

$./vipss -i input_file_name [-l user_lambda] [-s number_voxel_per_line] [-o output_file_path]

where:
1. -i: followed by the path of the input file. input_file_name is a path to the input file. currently, support file format includes ".xyz". The format of .xyz is:  
x1, y1, z1  
x2, y2, z2  
.....  
xn, yn, zn  

2. -l: optional argument. Followed by a float number indicating the lambda which balances the energy (see the paper for details). Default 0 (exact interpolation), you should set and tune this number according to your inputs.

3. -s: optional argument. Followed by a unsigned integer number indicating the number of voxels in each dimension for the implicit surfacing. Only If -s is included in the command line, the program would output the surface ([input file name]_surface.ply). We recomment using 100 for a default value, and you should set this according to your inputs and the precision of the output. Notices that the surfacing algorithm takes quite a long time for surfacing the zero-level set, and it depends on the resolution and the shape of the zero-level set.

4. -o: optional argument. Followed by the path of the output path. output_file_path is a path to the folder for generating output files. Default the folder of the input file.

5. -t: optional argument. when it is activated, the program will create a txt file ([input file name]_time.txt) which records the timing information in this run.

Some examples have been placed at data folder for testing:
1. $./vipss -i ../data/hand_ok/input.xyz -l 0 -s 200
2. $./vipss -i ../data/walrus/input.xyz -l 0.003 -s 100

The program will generate the predicted normal in [input file name]_normal.ply.
If -s is included in the command line, the program will generate the surface as the zero-level set of the solved implicit function ([input file name]_surface.ply).

:bell: To generate all the example in the paper, please run the makefigure.sh script in the vipss folder:  
$source makefigure.sh  
The result will be generated into the data folder respectively.

:bell: Notice: The current surface tracker does NOT producing multi-component surface. We refer the user to use CGAL implicit mesher if needed. We will update the program to it soon.

:mega: For further questions about the code and the paper, please contact Zhiyang Huang at adshhzy@gmail.com or zhiyang.huang@wustl.edu (might be invalid after he graduated). You can also contact Prof. Tao Ju at taoju@wustl.edu.



