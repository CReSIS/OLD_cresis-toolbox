// MCMC-based ice layer finding
// (c) 2014 Stefan Lee, Jerome Mitchell, David Crandall
//
// Please see README for instructions on how to use, 
// LICENSE for licensing information.
//
// For details on the technique, please see:
// Stefan Lee Jerome Mitchell David J. Crandall Geoffrey C. Fox, 
//    "ESTIMATING BEDROCK AND SURFACE LAYER BOUNDARIES AND CONFIDENCE 
//     INTERVALS IN ICE SHEET RADAR IMAGERY USING MCMC", in International 
//     Conference on Pattern Recognition 2014.
//
// Also see our project website: vision.soic.indiana.edu/icelayers
//
// This work was supported in part by the National Science Foundation (CNS-0723054, 
//     OCI-0636361, IIS-1253549, ANT-0424589) and by a NASA Earth and Space Science 
//     Fellowship (NNX13AN82H). Thanks to the Center for the Remote Sensing of 
//     Ice Sheets (CReSIS) for providing datasets.
//
//

#include "SImage.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <limits>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "ranker.h"
#include <utility>
#include <map>

using namespace std;
typedef vector< pair<int, int> > PathType;
/////////////////////////////////////////////////////////////////////////////////////////////////
//Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////

#define ITERS 30000
#define BURN_IN 20000

//Small Images
//#define OPT_SUB_WINDOW 30
//#define NON_MAX_WIDTH 30
//#define HORZ_STDEV 10
//#define VERT_THRESH 5
//#define SPAN 3


#define NON_MAX_WIDTH 50
#define HORZ_STDEV 75
#define OPT_SUB_WINDOW 20
#define VERT_THRESH 25
#define OFFSET 3
#define SPAN 5

#define VERT_THRESH_VALUE 0.1
#define SMALL_CONSTANT 0.00001

/////////////////////////////////////////////////////////////////////////////////////////////////

//Utility
double **cacheV, **cacheH; 
double normal_pdf(double x, double m, double s);    
double unifRand();
void setIJParam(int val,int i, int j, int *parameters,SDoublePlane &input);
int getIJParam(int i, int j, int *parameters,const SDoublePlane &input);
void precomputePairwiseVert(int rows);
void precomputePairwiseHort(int rows);

//Pairwise Constraints
double pairwiseVertical(int r1, int r2);
double pairwiseHorizontal(int r1, int r2);

//Gibbs Sampler Core
int scene_labeling(SDoublePlane &input, int init_k, int **parameters,int **lower,int **upper, const PathType &pts1, const PathType &pts2);
void gibbsUpdate(int *parameters, int k, SDoublePlane &input, int *fixed);
int gibbsSampleIJ(int i, int j, int k, int *parameters, SDoublePlane &input, int *fixed);
int optGibbsSampleIJ(int i, int j, int k, int *parameters, SDoublePlane &input, int *fixed);
void conditionHorizontal(int r, int rows, double* dist, int offset);
void conditionVertical(int r, int rows, double* dist, int top, int offset);

//Reversible Jumps
void birth(int *parameters,int *k, SDoublePlane &input);
void death(int *parameters, int *k, SDoublePlane &input);
double computePosteriorProbability(int *parameters, int k, int n, SDoublePlane &input);

//Image IO 
SDoublePlane resample(const SDoublePlane &input, int newRow, int newCol);
void write_results(const SDoublePlane &input, int *parameters, const char* filename, int k);
void write_plane(const SDoublePlane &unary_top, const char *fname);

/////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef STAND_ALONE
#include <SImageIO.h>
int main(int argc, char *argv[])
{
    srand (time(NULL));
    string input_filename = argv[1], output_filename = argv[2];
    string basename_g = output_filename.c_str();


    PathType a,b;
    //SDoublePlane input = resample(SImageIO::read_png_file(input_filename.c_str()), 140,180);  
    SDoublePlane input = SImageIO::read_png_file(input_filename.c_str());
    int *parameters, *upper, *lower;
    int k = scene_labeling(input, 2, &parameters, &lower, &upper,a,b);




    SDoublePlane R, G, B;
    R = input;  G = input;  B = input;
    for (int c=input.cols()-1; c>=0; c--){
        for(int r = 0; r < k; r++){
            int u = getIJParam(r, c, upper,input);			
            int t = getIJParam(r, c, parameters,input);		
            int l = getIJParam(r, c, lower, input);
            if(r%2 == 0){
                R[u][c] = 255.0;G[u][c] = 125; B[u][c]=125;
                R[t][c] = 255.0;G[t][c] = 0; B[t][c]=0;
                R[l][c] = 255.0;G[l][c] = 125; B[l][c]=125;
            }
            else{
                R[u][c] = 125;G[u][c] = 125; B[u][c]=255.0;
                R[t][c] = 0;G[t][c] = 0; B[t][c]=255.0;
                R[l][c] = 125;G[l][c] = 125; B[l][c]=255.0;
            }
        }
    }

    SImageIO::write_png_file(output_filename.c_str(), R, G, B);

    //Compare to ground truth
    if(argc == 4){
        double totalError = 0;
        int numPoints = 0;	
        int surface[input.cols()];
        int bedrock[input.cols()];
        SImageIO::read_png_file(argv[3],R,G,B);
        for(int c = 0; c < input.cols(); c++){
            surface[c] = -1;
            bedrock[c] = -1;			
            for(int r = 0; r < input.rows(); r++)
                if(bedrock[c] == -1 && R[r][c] == 255 && B[r][c] == 0 && G[r][c]==0){
                    bedrock[c] = r;
                    numPoints++;
                }
                else if(surface[c] == -1 && R[r][c] == 255 && B[r][c] == 255 & G[r][c]==0){
                    surface[c] = r;
                    numPoints++;
                }
        }

        for(int c = 0; c< input.cols(); c++){
            if(surface[c] >=0)
                totalError += abs(surface[c]-getIJParam(0,c,parameters,input));
            if(bedrock[c] >=0)
                totalError += abs(bedrock[c]-getIJParam(1,c,parameters,input));
        }
        cout << argv[3] << ", " << numPoints << ", " << totalError << ", "<< totalError/(double)numPoints << endl;	
        cout.flush();
    }



    return 0;
}

#else



#include "mex.h" 
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    double *input = mxGetPr(prhs[0]);
    int rows= mxGetM(prhs[0]), cols = mxGetN(prhs[0]);
    SDoublePlane P(rows, cols);

    for(int j=0, cp=0; j<cols; j++){
        for(int i=0; i<rows; i++, cp++){
            P[i][j] = input[cp];
        }
    }

    PathType pts1, pts2;
    if(nrhs > 1)
    {
        int m = mxGetN(prhs[1]);
        double *t = mxGetPr(prhs[1]);

        for(int i=0; i<m; i++){
            pts1.push_back( pair<int, int>((int)t[i*2], (int)t[i*2+1]) );
        }
        printf("%s\n", "IN if1");
        m = mxGetN(prhs[2]);
        t = mxGetPr(prhs[2]);
        printf("%s\n", "IN if2");
        for(int i=0; i<m; i++){
            pts2.push_back( pair<int, int>((int)t[i*2], (int)t[i*2+1]) );}
        printf("%s\n", "IN if3");
    }

    int *parameters, *upper, *lower;
    int k = scene_labeling(P, 2, &parameters, &lower, &upper, pts1, pts2);

    plhs[0] = mxCreateDoubleMatrix(2, cols, mxREAL);
    double *P_ptr = mxGetPr(plhs[0]);
    for(int i=0; i<cols; i++){
        P_ptr[i*2] = getIJParam(0, i, parameters,P);
        P_ptr[i*2+1] = getIJParam(1,i,parameters,P);
    }
    plhs[1] = mxCreateDoubleMatrix(2, cols, mxREAL);
    P_ptr = mxGetPr(plhs[1]);
    for(int i=0; i<cols; i++){
        P_ptr[i*2] = getIJParam(0, i, upper,P);
        P_ptr[i*2+1] = getIJParam(1,i,upper,P);
    }

    plhs[2] = mxCreateDoubleMatrix(2, cols, mxREAL);
    P_ptr = mxGetPr(plhs[2]);
    for(int i=0; i<cols; i++){
        P_ptr[i*2] = getIJParam(0, i, lower,P);
        P_ptr[i*2+1] = getIJParam(1,i,lower,P);
    }

    delete[] parameters;
    delete[] upper;
    delete[] lower;
}

#endif


int maxIndex(double *r, int rows, int nonMax){
    double maxV = -1;
    int maxI = -1;
    for(int i = 0; i < rows; i++)
        if(r[i] > maxV){
            maxI = i;
            maxV = r[i];
        }
    int upper = 0;//maxI-nonMax>0 ? maxI-nonMax : 0;
    int lower = maxI+nonMax<rows ? maxI+nonMax : rows-1;
    for(int i = upper; i <= lower; i++)
        r[i] = -1;

    return maxI;
}


int scene_labeling(SDoublePlane &input, int init_k, int **parameters, int **lower, int **upper, const PathType &pts1, const PathType &pts2){
    int rows = input.rows(), cols = input.cols();

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Caching
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //Precompute Pairwise Probabilities
    precomputePairwiseVert(rows);
    precomputePairwiseHort(rows);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Prior Calculation
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //Compute Gradient Feature
    SDoublePlane gradients=SDoublePlane(input.rows(),input.cols());
    for (int i=0;i<rows;i++){
        int low_i = i-SPAN>=1 ? i-SPAN: 1;
        int high_i = i+SPAN<rows?i+SPAN:rows-1;
        for (int k=0;k<cols;k++){
            double gradient = 0.0;
            for (int j=low_i;j<high_i;j++){
                gradient+=fabs(input[j][k]-input[high_i-(j-low_i)][k]);
            }
            int low_j = k-SPAN>=0 ? k-SPAN:1;
            int high_j = k+SPAN<cols ? k+SPAN:cols-1;
            for (int j=low_j;j<=high_j;j++){
                gradient+=0.5*fabs(input[i][j]-input[i][high_j-(j-low_j)]);
            }
            gradients[i][k]=(gradient);
        }
    }

    SDoublePlane prior=SDoublePlane(rows,cols);
    for(int i = 0; i < cols; i++){
        double total = 0;
        double totalc = 0;//
        for(int j = 0; j < rows; j++){
            total+= gradients[j][i];
            totalc+= 255-input[j][i];//
        }
        for(int j = 0; j < rows; j++){
            gradients[j][i] = gradients[j][i]/total;
            prior[j][i] = max(SMALL_CONSTANT,gradients[j][i]*max(0.01,(255-input[j][i])/totalc));
        }	
    }

    //Normalize Cols
    for(int i = 0; i < cols; i++){
        double total = 0;
        for(int j = 0; j < rows; j++){
            total+= prior[j][i];
        }
        for(int j = 0; j < rows; j++){
            prior[j][i] = prior[j][i]/total;
        }	
    }

    //Output prior specification image
#ifdef DEBUG
    cout << "Outputting prior..";	
    cout.flush();
    write_plane(prior,"run/prior.png");
    cout << "done." << endl;
    cout.flush(); 
#endif


    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Line Location Initialization
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //Initialize lines to row maximums of the prior with non-max suppression
    int *fixed;
    try{	
        *parameters = new int[cols*init_k];	
        fixed = new int[cols*init_k];
        for(int i = 0; i < cols*init_k; i++)
            fixed[i] = 0;
    }catch(std::bad_alloc& exc)
    {
        cout << "Failed to malloc... terminating." << endl;
        cout.flush();
        exit(0);
    }	

    //Compute row sum of prior probability	
    double rowLike[rows];
    for(int j = 0; j < rows; j++){
        rowLike[j]=0;		
        for(int i = 0; i < cols; i++)
            rowLike[j] += prior[j][i];
    }

    SDoublePlane temp(prior);
    double big, bigger;
    int big_c, bigger_c;
    for(int c = 0; c < cols; c++){
        big = -2;bigger = -1;		

        //find biggest and supress		
        for(int r = 0; r < rows; r++){
            if(temp[r][c] > bigger){				
                bigger = temp[r][c];
                bigger_c = r;
            }
        }
        int upper = 0;
        int lower = bigger_c+NON_MAX_WIDTH<rows ? bigger_c+NON_MAX_WIDTH : rows-1;
        for(int i = upper; i <= lower; i++)
            temp[i][c] = -1;


        //find biggest and supress		
        for(int r = 0; r < rows; r++){
            if(temp[r][c] >= big){				
                big = temp[r][c];
                big_c = r;
            }
        }

        if(bigger_c < big_c){
            setIJParam(bigger_c, 0, c, *parameters, prior);
            setIJParam(big_c, 1, c, *parameters, prior);
        }
        else{
            setIJParam(bigger_c, 1, c, *parameters, prior);
            setIJParam(big_c, 0, c, *parameters, prior);
        }

    }


    /*
    //Get k maxima with supression
    double maxima[init_k];
    for(int k = 0; k < init_k; k++)
    maxima[k] = maxIndex(rowLike, rows,NON_MAX_WIDTH);


    //Get in order of lines (i.e. highest goes first line!)
    for(int k = 0; k < init_k; k++){
    int index; double minV = rows+1;
    for(int t = 0; t < init_k; t++)
    if(maxima[t] < minV){
    index = t;
    minV = maxima[t];
    }
    maxima[index] = rows+1;
    for(int j = 0; j <  cols; j++)
    setIJParam(minV,k,j,*parameters, prior);	}

*/
    for(int i = 0; i < pts1.size(); i++){
        cerr << pts1[i].second-OFFSET << " " << pts1[i].first-1 << endl;
        cerr.flush();
        setIJParam(max(0,pts1[i].second-OFFSET-1), 0, max(0,pts1[i].first-1), *parameters, prior);
        setIJParam(1, 0, max(0,pts1[i].first-1), fixed, prior);
    }
    for(int i = 0; i < pts2.size(); i++){
        cerr << pts2[i].second-OFFSET << " " << pts2[i].first-1 << endl;
        cerr.flush();		
        setIJParam(max(0,pts2[i].second-OFFSET-1), 0, max(0,pts2[i].first-1), *parameters, prior);
        setIJParam(1, 1, max(0,pts2[i].first-1), fixed, prior);		
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Gibbs Sampler
    /////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
    cout << "Starting sampler..."<<endl;
    cout.flush();
#endif

    //Allocate memory for storing distribution after BURN_IN period
    double *mean, **dist;
    try{	
        mean = new double[cols*init_k];
        bzero(mean, sizeof(double)*cols*init_k);
        dist = new double*[init_k*cols];
        for(int i = 0; i < init_k*cols; ++i)
            dist[i] = new double[ITERS-BURN_IN];	
    }catch(std::bad_alloc& exc)
    {
        cout << "Failed to malloc... terminating." << endl;
        cout.flush();
        exit(0);
    }	


    //Run Sampler
    int k = init_k;
    for(int i = 0; i < ITERS; i++){
#ifdef DEBUG		
        if(i%1000 == 0){
            char buffer[10];		
            sprintf(buffer,"run/out_%d.png",i);
            cout << "Writing state " << i << "...";				
            write_results(prior,*parameters, buffer, init_k);
            cout << "done." << endl;				
        }
#endif
        //cerr << i << endl;
        gibbsUpdate(*parameters, init_k, prior, fixed);
        if(i >= BURN_IN){	
            for(int j = 0; j < cols*k; j++){	
                dist[j][i-BURN_IN] = (*parameters)[j];			
                mean[j] += (*parameters)[j];
            }
        }		
    }

    //Take expectation
    for(int j = 0; j < input.cols()*k; j++){
        (*parameters)[j] = min((double)rows-1,floor(mean[j]/(ITERS-BURN_IN))+OFFSET);	
    }

    //Find confidence interval	
    *lower = new int[k*cols];
    *upper = new int[k*cols];
    for(int j = 0; j < input.cols()*k; j++){
        (*lower)[j] = max(0,(int)quantile(dist[j], ITERS-BURN_IN, 0.025)+OFFSET);
        (*parameters)[j] = max(0,min(rows-1,(int)floor(mean[j]/(ITERS-BURN_IN))+OFFSET));
        (*upper)[j] = min(rows-1,(int)quantile(dist[j], ITERS-BURN_IN, 0.975)+OFFSET);
    }
    delete[] fixed;
    delete[] mean;
    for(int i = 0; i < init_k*cols; i++)
        delete[] dist[i];
    delete [] dist;

    return k;
}


//TODO needed for birth/death steps
double computePosteriorProbability(int *parameters, int k, int n, SDoublePlane &input){


}

//TODO will likely need to adjust the way distribution storage is handled
void birth(int *parameters,int *k, SDoublePlane &input){


}

//TODO will likely need to adjust the way distribution storage is handled
void death(int *parameters, int *k, SDoublePlane &input){


}


void gibbsUpdate(int *parameters, int k, SDoublePlane &input, int *fixed){
    int new_param;	
    for(int i = 0; i < k; i++){
        for(int j = 0; j < input.cols(); j++){
#ifdef OPT			
            new_param = optGibbsSampleIJ(i, j, k, parameters, input, fixed);
#else
            new_param = gibbsSampleIJ(i, j, k, parameters, input, fixed);
#endif
            setIJParam(new_param, i, j, parameters, input);
        }
    }
}



void inline setIJParam(int val,int i, int j, int *parameters,SDoublePlane &input){
    parameters[i*input.cols()+j] = max(0,min(input.rows()-1,val));
}

int inline getIJParam(int i, int j, int *parameters,const SDoublePlane &input){
    return parameters[i*input.cols()+j];
}



//Resample the column assignment for the Jth node in the Ith line
int inline gibbsSampleIJ(int i, int j, int k, int *parameters, SDoublePlane &input, int *fixed){
    double dist[input.rows()];
    int n = input.cols();


    if(getIJParam(i,j, fixed, input)==1)
        return getIJParam(i,j,parameters,input);

    //Initialize unary
    for(int m = 0; m < input.rows(); m++){
        dist[m] = input[m][j];
    }

    //Condition on neighbors
    if(i+1 < k){
        //Condition on (i+1, j)
        conditionVertical(getIJParam(i+1, j, parameters, input), input.rows(), dist,0,0);
    }

    if(i-1 >= 0){
        //Condition on (i-1, j)
        conditionVertical(getIJParam(i-1, j, parameters, input), input.rows(), dist,1,0);
    }

    if(j+1 < input.cols()){
        //Condition on (i, j+1)
        conditionHorizontal(getIJParam(i, j+1, parameters, input), input.rows(), dist,0);
    }

    if(j-1 >= 0){
        //Condition on (i, j-1)
        conditionHorizontal(getIJParam(i, j-1, parameters, input), input.rows(), dist,0);
    }

    double total = 0;
    for(int m = 0; m < input.rows(); m++){
        total += dist[m];		
    }

    //Sample
    double u = unifRand()*total;
    total = 0;
    int m;
    for(m = 0; m < input.rows(); m++){
        total += dist[m];	
        if(total >= u){
            return m;
        }
    }

    return m-1;
} 



int inline optGibbsSampleIJ(int i, int j, int k, int *parameters, SDoublePlane &input, int *fixed){

    if(getIJParam(i,j, fixed, input) == 1){

        return getIJParam(i,j,parameters,input);

    }
    int n = input.rows();
    int upper, lower, leftU,leftL, rightU,rightL;
    leftU = j+1<n ? min(n, getIJParam(i,j+1,parameters,input)+OPT_SUB_WINDOW) : 0;
    leftL = j+1<n ? max(0, getIJParam(i,j+1,parameters,input)-OPT_SUB_WINDOW) : n;
    rightU = j-1>=0 ? min(n, getIJParam(i,j-1,parameters,input)+OPT_SUB_WINDOW) : 0;
    rightL = j-1>=0 ? max(0, getIJParam(i,j-1,parameters,input)-OPT_SUB_WINDOW) : n;

    upper = max(leftU,rightU);
    lower = min(leftL, rightL);

    n = upper-lower;
    double dist[n];

    //Initialize unary
    for(int m = 0; m < n; m++){
        dist[m] = input[lower+m][j];
    }

    //Condition on neighbors
    if(i+1 < k){
        //Condition on (i+1, j)
        conditionVertical(getIJParam(i+1, j, parameters, input), n, dist,0,lower);
    }

    if(i-1 >= 0){
        //Condition on (i-1, j)
        conditionVertical(getIJParam(i-1, j, parameters, input), n, dist,1,lower);
    }

    if(j+1 < input.cols()){
        //Condition on (i, j+1)
        conditionHorizontal(getIJParam(i, j+1, parameters, input), n, dist,lower);
    }

    if(j-1 >= 0){
        //Condition on (i, j-1)
        conditionHorizontal(getIJParam(i, j-1, parameters, input), n, dist,lower);
    }

    //Compute normalization constant	
    double total = 0;
    for(int m = 0; m < n; m++)
        total += dist[m];		


    //Sample from conditional
    double u = unifRand()*total;
    total = 0;
    int m;
    for(m = 0; m < n; m++){
        total += dist[m];	
        if(total >= u){
            return lower+m;
        }
    }
    return lower+m-1;
} 


void inline conditionHorizontal(int r, int rows, double* dist, int offset){	
    double *c = cacheH[r];
    for(int i = 0; i < rows; i++){
        dist[i] *= c[offset+i];
    }
}

void inline conditionVertical(int r, int rows, double* dist, int top, int offset){
    if(top){
        double *c = cacheV[r];
        for(int i = 0; i < rows; i++){
            dist[i] *= c[offset+i];
        }
    }else
        for(int i = 0; i < rows; i++){
            dist[i] *= cacheV[offset+i][r];
        }
}

void inline precomputePairwiseVert(int rows){
    cacheV = new double*[rows];
    for(int i = 0; i < rows; ++i)
        cacheV[i] = new double[rows];	

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < rows; j++)
            cacheV[i][j] = (pairwiseVertical(i, j));
}

void inline precomputePairwiseHort(int rows){
    cacheH = new double*[rows];
    for(int i = 0; i < rows; ++i)
        cacheH[i] = new double[rows];	

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < rows; j++)
            cacheH[i][j] = (pairwiseHorizontal(i, j));
}



double pairwiseVertical(int r1, int r2){ return r2-r1 < VERT_THRESH ? VERT_THRESH_VALUE : 1; }
double pairwiseHorizontal(int r1, int r2){	return max(SMALL_CONSTANT,normal_pdf(r1-r2, 0, sqrt(HORZ_STDEV)));}

double normal_pdf(double x, double m, double s)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * exp(-0.5f * a * a);
}
double unifRand(){ return rand() / double(RAND_MAX);}

#ifdef STAND_ALONE
void write_plane(const SDoublePlane &unary_top, const char *fname)
{


    SDoublePlane g = unary_top;
    double x=-INFINITY, n=INFINITY;
    for(int i=0; i<g.rows(); i++)
        for(int j=0; j<g.cols(); j++)
        {
            if(!isinf(g[i][j]))
                x = max(x, g[i][j]), n=min(n, g[i][j]);
        }
    for(int i=0; i<g.rows(); i++)
        for(int j=0; j<g.cols(); j++)
        {
            if(!isinf(g[i][j]))
                g[i][j] = (g[i][j] - n) / (x-n)  * 255.0;
            else if(isinf(g[i][j]) > 0)
                g[i][j] = 255;
            else
                g[i][j] = 0;
        }

    SImageIO::write_png_file(fname, g, g, g);

}


void write_results(const SDoublePlane &input, int *parameters, const char* filename, int k){

    SDoublePlane g = input;
    double x=-INFINITY, n=INFINITY;
    for(int i=0; i<g.rows(); i++)
        for(int j=0; j<g.cols(); j++)
        {
            if(!isinf(g[i][j]))
                x = max(x, g[i][j]), n=min(n, g[i][j]);
        }
    for(int i=0; i<g.rows(); i++)
        for(int j=0; j<g.cols(); j++)
        {
            if(!isinf(g[i][j]))
                g[i][j] = (g[i][j] - n) / (x-n)  * 255.0;
            else if(isinf(g[i][j]) > 0)
                g[i][j] = 255;
            else
                g[i][j] = 0;
        }
    SDoublePlane R,G,B;
    R = g; G=g; B=g;
    for (int c=input.cols()-1; c>=0; c--){ 
        for(int r = 0; r < k; r++){			
            int t = getIJParam(r, c, parameters,input);

            if(r%2 == 0){
                R[t][c] = 255.0;G[t][c] = 0; B[t][c]=0;}
            else{
                R[t][c] = 0;G[t][c] = 0; B[t][c]=255.0;}
        }
    }
    SImageIO::write_png_file(filename, R, G, B);
}

#endif

//Changes size of an image
SDoublePlane resample(const SDoublePlane &input, int newRow, int newCol)
{
    SDoublePlane output(newRow, newCol);
    int oldRow = input.rows();
    int oldCol = input.cols();
    for (int i=0;i<newRow;i++)
    {
        for (int j=0;j<newCol;j++)
        {
            double src_i = ((double)i)*oldRow/newRow;
            double src_j = ((double)j)*oldCol/newCol;
            int fi0 = int(src_i);
            int fi1 = fi0+1;
            int fj0 = int(src_j);
            int fj1 = fj0+1;
            double a = src_i - int(src_i);
            double b = src_j - int(src_j);
            output[i][j] = (1-b)*(1-a)*input[fi0][fj0]+(1-b)*a*input[fi1][fj0]+b*(1-a)*input[fi0][fj1]+b*a*input[fi1][fj1];
        }
    }
    return output;
}



