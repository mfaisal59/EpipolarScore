# EpipolarScore

This repository is part of public implementation of our "Exploiting Geometric Constraints on Dense Trajectories for Motion Saliency" paper. This repository contains Optical Flow & Epipolar Score Computation code. You can check the main repository [here](https://github.com/mfaisal59/EpONet).

### Installations:

This source code is based on MATLAB framework and tested on Ubuntu 16.04 with MATLAB 2016b.

### Instructions:

##### 1) Clone the repository
	
```
git clone https://github.com/mfaisal59/EpipolarScore.git
```

##### 2) Download Dataset

Download and unpack the DAVIS 2016 dataset and as well as the evaluatio code from https://davischallenge.org/davis2016/code.html

##### 3) Compute Optical Flow

The optical flow is based on Full Flow Method (https://cqf.io/fullflow/). To compute the optical flow for DAVIS Dataset, run the following script:

```
cd ./Full_Flow_Source_Code/
run davisBatch.m file
#modify the path to DAVIS dataset directory
```

##### 4) Compute Epipolar Score

To compute the Epipolar Score, modify the paths in 'testDAVIS.m' file and run:
```
cd ./EpipolarScoreMain/
run testDAVIS.m script
#modify the path to DAVIS dataset, forward and backward optical flow directory
```

##### 5) Convert Flow to X-Y Displacement Images

```
cd ./EpipolarScoreMain/
run flow2Displacement.m script
#modify the paths to DAVIS dataset, forward and backward optical flow directory
```

##### 6) Generate Motion Images

```
cd ./EpipolarScoreMain/
run generateMotionImages.m script
#modify the paths to DAVIS dataset and Optical Flow directory
```

