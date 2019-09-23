g++ -m64 -O3 -fopenmp -msse4 -shared -w -I"D:/Matlab_R2015b/extern/include" -DMATLAB_MEX_FILE -o ../../trws_l1.mexw64 *.cpp -L"D:/Matlab_R2015b/bin/win64" -lmex -lmx -leng -lmat
pause


