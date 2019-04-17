fcm - Zhiyang Huang (2018)

------------------------------------

This code implements the algorithm described in

Huang, Zhi Yang, et al. "Variational Implicit Point Set Surface" Computer Graphics Forum. Vol. 37. No. 2. 2018.

The primary function of this program is to predict the normal as well as the underlying surface of a given set of unoriented points.

Currently, the code is only tested on Mac OS X 10.10 or above, it should work at other platforms with minor modification.


BUILDING
======================================================================================================


The code has only two dependencies: 1)Armadillo,   2)NLOPT,  3)Suite-Sparse

1) http://arma.sourceforge.net/

2) https://nlopt.readthedocs.io/en/latest/

After download/install the three dependencies, please go to folder vipss and modify the relevant path according to your installation in the CMakeList.txt, e.g.

Armadillo include file path:    SET(ARMADILLO_INCLUDE_DIRS "/opt/local/include/")
Armadillo lib file path:    SET(ARMADILLO_LIB_DIRS "/opt/local/lib/")

NLOPT include file and lib path:
SET(NLOPT_INCLUDE_DIRS "/opt/local/include/")
SET(NLOPT_LIB_DIR "/opt/local/lib/")


Then build the Cmake file and make:
$cmake .
$make

In the vipss directory, there should be an executable called "vipss" (or "vipss.exe" on Windows if it is successfully built).


RUNNING
======================================================================================================

To run the code from the command line, type:

$./vipss -i input_file_name -o output_file_path [-s] [-v number_voxel_per_line] [-l user_lambda] 

where:
1. -i: followed by the path of the input file. input_file_name is a path to the input file. currently, support file format includes ".xyz". The format of .xyz is:
	x1, y1, z1,
	x2, y2, z2,
	.....
	xn, yn, zn
	
2. -o: followed by the path of the output path. output_file_path is a path to the folder for generating output files. Default the folder of the input file.

3. -s: optional argument. If -w is included in the command line, the program output the surface (predicted_surface.obj).

4. -v: optional argument. Followed by a unsigned integer number indicating the number of voxels in each dimension for the implicit surfacing. Default 100, you should set this according to your inputs, you can check the triangulated result in predicted_surface.obj by turn on -s option

5. -l: optional argument. Followed by a float number indicating the lambda which balances the energy (see the paper for details). Default 0 (exact interpolation), you should set and tune this number according to your inputs.


An example is placed at data folder for testing:
1. $./vipss -i ../data/hand_ok/input.xyz -l 0 -s -v 200

The program will generate the predict normal in predicted_normal.ply 


We are still reconstructing and optimizing the code from the original research code, thus some details in the paper has not been fully reproduced. A more compact and friendly version with visualization code and more examples will be updated soon.



