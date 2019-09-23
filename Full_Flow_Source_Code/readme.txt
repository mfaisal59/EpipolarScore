Please run and read "demo.m" for demonstration. The program was tested in Linux.

Here are the settings to our program:
-"im1" and "im2" are two input images.
-"maxDisp" is set to 100 in "demo.m" for fast processing. The optimzation takes about 5 seconds when maxDisp=100. Notice that we use maxDisp=242 in our experiments in the paper.
-"opt" includes the parameters for fullflow and Epicflow. The parameters used in the Sintel and Kitti evaluation are given in "demo.m".
-To use your own data term and regularization term, please modify "fullflow.m" accordingly.
-The solver "trws_l1" in "fullflow.m" only works for the L1 regurization. If the regularization is not L1 but convex, please use "trws_convex".
-If you want to recompile the mex solver "trws_l1", run "trws_l1/linux/mex_command" or "trws_l1/windows/mex_file.bat". Repeat the same steps for "trws_convex".

If you have any question, please email me at cqf@stanford.edu.
