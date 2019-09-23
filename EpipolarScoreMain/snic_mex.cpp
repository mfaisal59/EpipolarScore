//=================================================================================
//  snic_mex.cpp
//
//
//  AUTORIGHTS
//  Copyright (C) 2016 Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland.
//
//  Created by Radhakrishna Achanta on 05/November/16 (firstname.lastname@epfl.ch)
//
//
//  Code released for research purposes only. For commercial purposes, please
//  contact the author.
//=================================================================================
/*Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met
 
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of EPFL nor the names of its contributors may
 be used to endorse or promote products derived from this software
 without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <mex.h>
#include <cmath>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <queue>
#include <string.h>

using namespace std;

void rgbtolab(int* rin, int* gin, int* bin, int sz, double* lvec, double* avec, double* bvec)
{
    int i; int sR, sG, sB;
    double R,G,B;
    double X,Y,Z;
    double r, g, b;
    const double epsilon = 0.008856;	//actual CIE standard
    const double kappa   = 903.3;		//actual CIE standard
    
    const double Xr = 0.950456;	//reference white
    const double Yr = 1.0;		//reference white
    const double Zr = 1.088754;	//reference white
    double xr,yr,zr;
    double fx, fy, fz;
    double lval,aval,bval;
    
    for(i = 0; i < sz; i++)
    {
        sR = rin[i]; sG = gin[i]; sB = bin[i];
        R = sR/255.0;
        G = sG/255.0;
        B = sB/255.0;
        
        if(R <= 0.04045)	r = R/12.92;
        else				r = pow((R+0.055)/1.055,2.4);
        if(G <= 0.04045)	g = G/12.92;
        else				g = pow((G+0.055)/1.055,2.4);
        if(B <= 0.04045)	b = B/12.92;
        else				b = pow((B+0.055)/1.055,2.4);
        
        X = r*0.4124564 + g*0.3575761 + b*0.1804375;
        Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
        Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
        
        //------------------------
        // XYZ to LAB conversion
        //------------------------
        xr = X/Xr;
        yr = Y/Yr;
        zr = Z/Zr;
        
        if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
        else				fx = (kappa*xr + 16.0)/116.0;
        if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
        else				fy = (kappa*yr + 16.0)/116.0;
        if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
        else				fz = (kappa*zr + 16.0)/116.0;
        
        lval = 116.0*fy-16.0;
        aval = 500.0*(fx-fy);
        bval = 200.0*(fy-fz);
        
        lvec[i] = lval; avec[i] = aval; bvec[i] = bval;
    }
}

//===========================================================================
/// FindSeeds
//===========================================================================
void FindSeeds(const int width, const int height, int& numk, vector<int>& kx, vector<int>& ky)
{
    const int sz = width*height;
    int gridstep = sqrt(double(sz)/double(numk)) + 0.5;
    int halfstep = gridstep/2;
    double h = height; double w = width;
    
    int xsteps = int(width/gridstep);
    int ysteps = int(height/gridstep);
    int err1 = abs(xsteps*ysteps - numk);
    int err2 = abs(int(width/(gridstep-1))*int(height/(gridstep-1)) - numk);
    if(err2 < err1)
    {
        gridstep -= 1.0;
        xsteps = width /(gridstep);
        ysteps = height/(gridstep);
    }
    
    numk = (xsteps*ysteps);
    kx.resize(numk); ky.resize(numk);
    int n = 0;
    for(int y = halfstep, rowstep = 0; y < height && n < numk; y += gridstep, rowstep++)
    {
        for(int x = halfstep; x < width && n < numk; x += gridstep)
        {
            if( y <= h-halfstep && x <= w-halfstep)
            {
                kx[n] = x;
                ky[n] = y;
                n++;
            }
        }
    }
}


//===========================================================================
/// runSNIC
///
/// Runs the priority queue base Simple Non-Iterative Clustering (SNIC) algo.
//===========================================================================
void runSNIC(
             double*		lv,
             double*		av,
             double*		bv,
             const int					width,
             const int					height,
             int*                       labels,
             int*						outnumk,
             const int                  innumk,
             const double               compactness)
{
    const int w = width;
    const int h = height;
    const int sz = w*h;
    const int dx8[8] = {-1,  0, 1, 0, -1,  1, 1, -1};//for 4 or 8 connectivity
    const int dy8[8] = { 0, -1, 0, 1, -1, -1, 1,  1};//for 4 or 8 connectivity
    const int dn8[8] = {-1, -w, 1, w, -1-w,1-w,1+w,-1+w};
    
    struct NODE
    {
        unsigned int i; // the x and y values packed into one
        unsigned int k; // the label
        double d;       // the distance
    };
    struct compare
    {
        bool operator()(const NODE& one, const NODE& two)
        {
            return one.d > two.d;//for increasing order of distances
        }
    };
    //-------------
    // Find seeds
    //-------------
    vector<int> cx(0),cy(0);
    int numk = innumk;
    FindSeeds(width,height,numk,cx,cy);//the function may modify numk from its initial value
    //-------------
    // Initialize
    //-------------
    NODE tempnode;
    priority_queue<NODE, vector<NODE>, compare> pq;
    memset(labels,-1,sz*sizeof(int));
    for(int k = 0; k < numk; k++)
    {
        NODE tempnode;
        tempnode.i = cx[k] << 16 | cy[k];
        tempnode.k = k;
        tempnode.d = 0;
        pq.push(tempnode);
    }
    vector<double> kl(numk, 0), ka(numk, 0), kb(numk, 0);
    vector<double> kx(numk,0),ky(numk,0);
    vector<double> ksize(numk,0);
    
    const int CONNECTIVITY = 4;//values can be 4 or 8
    const double M = compactness;//10.0;
    const double invwt = (M*M*numk)/double(sz);
    
    int qlength = pq.size();
    int pixelcount = 0;
    int xx(0),yy(0),ii(0);
    double ldiff(0),adiff(0),bdiff(0),xdiff(0),ydiff(0),colordist(0),xydist(0),slicdist(0);
    //-------------
    // Run main loop
    //-------------
    while(qlength > 0) //while(nodevec.size() > 0)
    {
        NODE node = pq.top(); pq.pop(); qlength--;
        const int k = node.k;
        const int x = node.i >> 16 & 0xffff;
        const int y = node.i & 0xffff;
        const int i = y*width+x;
        
        if(labels[i] < 0)
        {
            labels[i] = k; pixelcount++;
            kl[k] += lv[i];
            ka[k] += av[i];
            kb[k] += bv[i];
            kx[k] += x;
            ky[k] += y;
            ksize[k] += 1.0;
            
            for(int p = 0; p < CONNECTIVITY; p++)
            {
                xx = x + dx8[p];
                yy = y + dy8[p];
                if(!(xx < 0 || xx >= w || yy < 0 || yy >= h))
                {
                    ii = i + dn8[p];
                    if(labels[ii] < 0)//create new nodes
                    {
                        ldiff = kl[k] - lv[ii]*ksize[k];
                        adiff = ka[k] - av[ii]*ksize[k];
                        bdiff = kb[k] - bv[ii]*ksize[k];
                        xdiff = kx[k] - xx*ksize[k];
                        ydiff = ky[k] - yy*ksize[k];
                        
                        colordist   = ldiff*ldiff + adiff*adiff + bdiff*bdiff;
                        xydist      = xdiff*xdiff + ydiff*ydiff;
                        slicdist    = (colordist + xydist*invwt)/(ksize[k]*ksize[k]);//late normalization by ksize[k], to have only one division operation
            
                        tempnode.i = xx << 16 | yy;
                        tempnode.k = k;
                        tempnode.d = slicdist;
                        pq.push(tempnode); qlength++;
                        
                    }
                }
            }
        }
    }
    *outnumk = numk;
    //---------------------------------------------
    // Label the rarely occuring unlabelled pixels
    //---------------------------------------------
    if(labels[0] < 0) labels[0] = 0;
    for(int y = 1; y < height; y++)
    {
        for(int x = 1; x < width; x++)
        {
            int i = y*width+x;
            if(labels[i] < 0)//find an adjacent label
            {
                if(labels[i-1] >= 0) labels[i] = labels[i-1];
                else if(labels[i-width] >= 0) labels[i] = labels[i-width];
            }//if labels[i] < 0 ends
        }
    }
    
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1) {
        mexErrMsgTxt("At least one argument is required.") ;
    } else if(nrhs > 3) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("SLIC:nlhs","Two outputs required, a labels and the number of labels, i.e superpixels.");
    }
    //---------------------------
    int numelements   = mxGetNumberOfElements(prhs[0]) ;
    mwSize numdims = mxGetNumberOfDimensions(prhs[0]) ;
    const mwSize* dims  = mxGetDimensions(prhs[0]) ;
    unsigned char* imgbytes  = (unsigned char*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it
    int width = dims[1];
    int height = dims[0];//Note: first dimension provided is height and second is width
    int sz = width*height;
    //---------------------------
    const int numk = mxGetScalar(prhs[1]);
    double compactness  = mxGetScalar(prhs[2]);
    //---------------------------
    // Allocate memory
    //---------------------------
    int* rin    = (int*)mxMalloc( sizeof(int)      * sz ) ;
    int* gin    = (int*)mxMalloc( sizeof(int)      * sz ) ;
    int* bin    = (int*)mxMalloc( sizeof(int)      * sz ) ;
    double* lvec    = (double*)mxMalloc( sizeof(double)      * sz ) ;
    double* avec    = (double*)mxMalloc( sizeof(double)      * sz ) ;
    double* bvec    = (double*)mxMalloc( sizeof(double)      * sz ) ;
    int* klabels = (int*)mxMalloc( sizeof(int)         * sz );//original k-means labels
    //---------------------------
    // Perform color conversion
    //---------------------------
    //if(2 == numdims)
    if(numelements/sz == 1)//if it is a grayscale image, copy the values directly into the lab vectors
    {
        for(int x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(int y = 0; y < height; y++)
            {
                int i = y*width+x;
                lvec[i] = imgbytes[ii];
                avec[i] = imgbytes[ii];
                bvec[i] = imgbytes[ii];
                ii++;
            }
        }
    }
    else//else covert from rgb to lab
    {
        if(1)//convert from rgb to cielab space
        {
            for(int x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
            {
                for(int y = 0; y < height; y++)
                {
                    int i = y*width+x;
                    rin[i] = imgbytes[ii];
                    gin[i] = imgbytes[ii+sz];
                    bin[i] = imgbytes[ii+sz+sz];
                    ii++;
                }
            }
            rgbtolab(rin,gin,bin,sz,lvec,avec,bvec);
        }
        else//else use rgb values directly
        {
            for(int x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
            {
                for(int y = 0; y < height; y++)
                {
                    int i = y*width+x;
                    lvec[i] = imgbytes[ii];
                    avec[i] = imgbytes[ii+sz];
                    bvec[i] = imgbytes[ii+sz+sz];
                    ii++;
                }
            }
        }
    }
    //---------------------------
    // Compute superpixels
    //---------------------------
    int numklabels = 0;
    runSNIC(lvec,avec,bvec,width,height,klabels,&numklabels,numk,compactness);
    //---------------------------
    // Assign output labels
    //---------------------------
    plhs[0] = mxCreateNumericMatrix(height,width,mxINT32_CLASS,mxREAL);
    int* outlabels = (int*)mxGetData(plhs[0]);
    for(int x = 0, ii = 0; x < width; x++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(int y = 0; y < height; y++)
        {
            int i = y*width+x;
            outlabels[ii] = klabels[i];
            ii++;
        }
    }
    //---------------------------
    // Assign number of labels/seeds
    //---------------------------
    plhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    int* outputNumSuperpixels = (int*)mxGetData(plhs[1]);//gives a void*, cast it to int*
    *outputNumSuperpixels = numklabels;
    //---------------------------
    // Deallocate memory
    //---------------------------
    mxFree(rin);
    mxFree(gin);
    mxFree(bin);
    mxFree(lvec);
    mxFree(avec);
    mxFree(bvec);
    mxFree(klabels);
}
