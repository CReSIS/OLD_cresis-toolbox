#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//


// Rotate an image about its center by the specified angle, using a nearest-neighbor approach
//  (i.e., no interpolation...)
//

double real_abs(double x)
{
  return x<0?-x:x;
}

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

SDoublePlane rotate_nn(const SDoublePlane &input, double angle)
{
  SDoublePlane output(input.rows(), input.cols());

  // rotation matrix = [ rot00 rot01; rot10 rot11 ]
  double rot00 = cos(angle),  rot01 = -sin(angle);
  double rot10 = sin(angle), rot11 = cos(angle);

  // invert rotation matrix
  double irot00 = rot00,  irot01 = rot10;
  double irot10 = rot01,  irot11 = rot11;

  for(int i=0; i<input.rows(); i++)
    for(int j=0; j<input.cols(); j++)
      {
	// move origin to center of image
	int dest_i = i - input.rows()/2, dest_j = j - input.cols()/2;

	int src_i = int( irot00 * dest_i + irot01 * dest_j );
	int src_j = int( irot10 * dest_i + irot11 * dest_j );

	src_i += input.rows()/2;   src_j += input.cols()/2;
	dest_i += input.rows()/2;  dest_j += input.cols()/2;

	// handle case where corresponding pixel in input image is 
	//  outside of the image boundaries...
	if(src_i >=0 && src_i < input.rows() && src_j >= 0 && src_j < input.cols())
	  output[dest_i][dest_j] = input[src_i][src_j];
	else
	  output[i][j] = 0;
      }

  return output;
}


// Rotate an image about its center by the specified angle, using bilinear interpolation
//
SDoublePlane rotate_bilinear(const SDoublePlane &input, double angle)
{
  SDoublePlane output(input.rows(), input.cols());

  // rotation matrix = [ rot00 rot01; rot10 rot11 ]
  double rot00 = cos(angle),  rot01 = -sin(angle);
  double rot10 = sin(angle), rot11 = cos(angle);

  // invert rotation matrix
  double irot00 = rot00,  irot01 = rot10;
  double irot10 = rot01,  irot11 = rot11;

  for(int i=0; i<input.rows(); i++)
    for(int j=0; j<input.cols(); j++)
      {
	// move origin to center of image
	int dest_i = i - input.rows()/2, dest_j = j - input.cols()/2;

	double src_i = irot00 * dest_i + irot01 * dest_j;
	double src_j = irot10 * dest_i + irot11 * dest_j;
	src_i += input.rows()/2;
	src_j += input.cols()/2;
	int fi0 = int(src_i);
	int fi1 = fi0+1;
	int fj0 = int(src_j);
	int fj1 = fj0+1;
	double a = src_i - int(src_i);
	double b = src_j - int(src_j);

	// handle case where corresponding pixel in input image is 
	//  outside of the image boundaries...
	if(fi0 >=0 && fi1 < input.rows() && fj0 >= 0 && fj1 < input.cols())
	  output[i][j] = (1-b)*(1-a)*input[fi0][fj0]+(1-b)*a*input[fi1][fj0]+b*(1-a)*input[fi0][fj1]+b*a*input[fi1][fj1];
	else
	  output[i][j] = 0;
      }

  return output;
}


// an SDoubleMatrix is exactly the same as an SDoublePlane; we'll just give it an
//  alias to make it clear when we're talking about images and when we're talking
//  about matrices.
//
typedef SDoublePlane SDoubleMatrix;

// Apply a given 3x3 projective transformation matrix (aka a homography) to an image
//  
SDoublePlane projective_transformation(const SDoublePlane &input, const SDoubleMatrix &homography)
{
  SDoublePlane output(input.rows(), input.cols());

  // homography matrix = [a b c; d e f; g h k;]
  double a = homography[0][0],b = homography[0][1], c = homography[0][2], d = homography[1][0], e = homography[1][1], f = homography[1][2], g = homography[2][0], h = homography[2][1], k = homography[2][2];
  double det_homography = a*e*k+b*f*g+d*h*c-a*f*h-b*d*k-c*e*g;

  double ia = (e*k-f*h)/det_homography;
  double ib = (f*g-k*d)/det_homography;
  double ic = (d*h-e*g)/det_homography;
  double id = (c*h-k*b)/det_homography;
  double ie = (a*k-c*g)/det_homography;
  double i_f = (b*g-a*h)/det_homography;
  double ig = (b*f-c*e)/det_homography;  
  double ih = (c*d-f*a)/det_homography;
  double ik = (a*e-d*b)/det_homography;

  for(int i=0; i<input.rows(); i++)
    for(int j=0; j<input.cols(); j++)
      {

	double src_i = ib * j + ie * i + ih;
	double src_j = ia * j + id * i + ig;
	double src_w = ic * j + i_f * i + ik;
    src_i /= src_w;
    src_j /= src_w;
	int fi0 = int(src_i);
	int fi1 = fi0+1;
	int fj0 = int(src_j);
	int fj1 = fj0+1;
	double a = src_i - int(src_i);
	double b = src_j - int(src_j);

	// handle case where corresponding pixel in input image is 
	//  outside of the image boundaries...
	if(fi0 >=0 && fi1 < input.rows() && fj0 >= 0 && fj1 < input.cols())
	  output[i][j] = (1-b)*(1-a)*input[fi0][fj0]+(1-b)*a*input[fi1][fj0]+b*(1-a)*input[fi0][fj1]+b*a*input[fi1][fj1];
	else
	  output[i][j] = 0;
      }

  return output;
}

SDoublePlane disparity_map(const SDoublePlane &input1, const SDoublePlane &input2)
{
  // implement this in step 4...
  //  this placeholder just returns a random disparity map
  SDoublePlane result(input1.rows(), input1.cols());
  int w = 3;
  int dmax = 50;
  for(int i=0; i<input1.rows(); i++)
    for(int j=0; j<input1.cols(); j++)
    {
      double mincost = 1e32;
      int d=0;
      int minD=0;
      for (d=0;d<=dmax;d++)
      {
        double cost = 0;
        int xmin = i-w>=0?i-w:0;
        int xmax = i+w<input1.rows()?i+w:input1.rows()-1;
        int ymin = j-w>=0?j-w:0;
        int ymax = j+w<input1.cols()?j+w:input1.cols()-1;
        for (int x = xmin; x<=xmax; x++)
        {
          for (int y = ymin; y<=ymax; y++)
          {
            if (y+d<input1.cols())
            {
              double l = input1[x][y];
              double r = input2[x][y+d];
              cost += (l-r)*(l-r);
            }
            else
            {
              double l = input1[x][y];
              double r = 0;
              cost += (l-r)*(l-r);
            }
          }
        }
        if (cost<mincost)
        {
          mincost = cost;
          minD = d;
        }
      }

//      result[i][j] = minD;
//		3x version
			result[i][j] = minD*3;	
    }
  return result;
}
//-------------------------hmm
double abs(double x)
{
	return x<0?-x:x;
}

void getDisparity(int i, int j, int dmax, int w, double* vectd, const SDoublePlane &input1, const SDoublePlane &input2)
{
	int d=0;
	for (d=0;d<dmax;d++)
	{
		double cost = 0;
		int xmin = i-w>=0?i-w:0;
		int xmax = i+w<input1.rows()?i+w:input1.rows()-1;
		int ymin = j-w>=0?j-w:0;
		int ymax = j+w<input1.cols()?j+w:input1.cols()-1;
		for (int x = xmin; x<=xmax; x++)
		{
			for (int y = ymin; y<=ymax; y++)
			{
				if (y+d<input1.cols())
				{
					double l = input1[x][y];
					double r = input2[x][y+d];
					cost += (l-r)*(l-r);
				}
				else
				{
					double l = input1[x][y];
					double r = 0;
					cost += (l-r)*(l-r);
				}
			}
		}
		vectd[d]=cost;
//		printf("%lf, ", cost);
	}
//	getchar();
}

double potts(int x, int y)
{
	return abs(x-y)<=3?0:10000;
}

SDoublePlane hmm_stereo(const SDoublePlane &input1, const SDoublePlane &input2)
{
  // implement this in step 5...
  //  this placeholder just returns the result of disparity_map();

  SDoublePlane result(input1.rows(), input1.cols());
  int w = 3;
  int dmax = 50;
  SDoublePlane mind(input1.cols(), dmax);
  SDoublePlane prevS(input1.cols(), dmax);
  double tmpd[dmax];
  int d;
  for(int i=0; i<input1.rows(); i++)
  {
		getDisparity(i, 0, dmax, w, tmpd, input1, input2);
		for (d=0;d<dmax;d++)
		{
			mind[0][d]=tmpd[d];
			prevS[0][d]=-1;
		}
		for (int j=1; j<input1.cols(); j++)
		{
			getDisparity(i, j, dmax, w, tmpd, input1, input2);
			for (d=0;d<dmax;d++)
			{
				double minCost = 1e32;
				int minPrevS=-1;
				for (int k=0;k<dmax;k++)
				{
					//double v = abs(d-k)*10000;
					double v=potts(d,k);
					double cost = tmpd[d]+v+mind[j-1][k];
					if (cost<minCost)
					{
						minCost = cost;
						minPrevS = k;
					}
				}
				mind[j][d]=minCost;
				prevS[j][d]=minPrevS;
//				printf("min cost=%lf , previous step=%d",minCost,minPrevS);
//				getchar();
			}
		}
		double minC = 1e32;
		int bestLast = -1;
		int lastCol = input1.cols()-1;
		for (d=0;d<dmax;d++)
		{
			if (mind[lastCol][d]<minC)
			{
				minC = mind[lastCol][d];
				bestLast = d;
			}
		}
		int prev = bestLast;
		for (int j=input1.cols()-1;j>=0;j--)
		{
//			printf("%d, ",prev);
			result[i][j] = prev*3;
			prev = (int)prevS[j][prev];
		}
//		printf("\n");
//		getchar();
  }
  return result;
}

//get result from matrix[i][j-1]
double get_left(const SDoublePlane &matrix, int i, int j)
{
  if (j-1>=0)
    return matrix[i][j-1];
  else
    return 0;
}

double get_right(const SDoublePlane &matrix, int i, int j)
{
  if (j+1<matrix.cols())
    return matrix[i][j+1];
  else
    return 0;
}

double get_up(const SDoublePlane &matrix, int i, int j)
{
  if (i-1>=0)
    return matrix[i-1][j];
  else
    return 0;
}

double get_down(const SDoublePlane &matrix, int i, int j)
{
  if (i+1<matrix.rows())
    return matrix[i+1][j];
  else
    return 0;
}

double min(double a, double b)
{
  return a>b?b:a;
}

SDoublePlane mrf_stereo(const SDoublePlane &inputleft, const SDoublePlane &inputright)
{
  int w = 3;
  int dmax = 30;
  int i,j,k;
  int d=0;
  double dist = 2.5e4;
  SDoublePlane input1 = resample(inputleft, inputleft.rows()/2, inputleft.cols()/2);
  SDoublePlane input2 = resample(inputright, inputright.rows()/2, inputright.cols()/2);
  SDoublePlane result(input1.rows(), input1.cols());
  SDoublePlane dMatrix(input1.rows()*input1.cols(),dmax);
  int max_iter = input1.cols()+input1.rows();
  //up[i][j] means message from node (i,j) to (i-1,j)
  SDoublePlane up[dmax];
  SDoublePlane tmp_up[dmax];
  printf("start initializing up...\n");
  for (i=0;i<dmax;i++)
  {
    up[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_up[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  //down[i][j] means message from node (i,j) to (i+1,j)
  SDoublePlane down[dmax];
  SDoublePlane tmp_down[dmax];
  printf("start initializing down...\n");
  for (i=0;i<dmax;i++)
  {
    down[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_down[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  //left[i][j] means message from node (i,j) to (i,j-1)
  SDoublePlane left[dmax];
  SDoublePlane tmp_left[dmax];
  printf("start initializing left...\n");
  for (i=0;i<dmax;i++)
  {
    left[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_left[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  //right[i][j] means message from node (i,j) to (i,j+1)
  SDoublePlane right[dmax];
  SDoublePlane tmp_right[dmax];
  printf("start initializing right...\n");
  for (i=0;i<dmax;i++)
  {
    right[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_right[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  printf("start initializing d-matrix...\n");
  for (i=0;i<input1.rows();i++)
    for (j=0;j<input1.cols();j++)
    {
	    for (d=0;d<dmax;d++)
	    {
		    double cost = 0;
		    int xmin = i-w>=0?i-w:0;
		    int xmax = i+w<input1.rows()?i+w:input1.rows()-1;
		    int ymin = j-w>=0?j-w:0;
		    int ymax = j+w<input1.cols()?j+w:input1.cols()-1;
		    for (int x = xmin; x<=xmax; x++)
		    {
			    for (int y = ymin; y<=ymax; y++)
			    {
				    if (y+d<input1.cols())
				    {
					    double l = input1[x][y];
					    double r = input2[x][y+d];
					    cost += (l-r)*(l-r);
				    }
				    else
				    {
					    double l = input1[x][y];
					    double r = 0;
					    cost += (l-r)*(l-r);
				    }
			    }
		    }
		    dMatrix[i*input1.cols()+j][d]=cost;
		    //printf("d-m[%d][%d]=%f\n",i,j,cost);
	    }
	  }

  double *up_h = new double[dmax], *down_h=new double[dmax], *left_h=new double[dmax],*right_h=new double[dmax];
  printf("Start calculating...\n");
//  for (int iter = 0;iter<dmax;iter++)
  for (int iter = 0;iter<max_iter;iter++)
  {
    if(iter%50==0)
      printf("Iteration #%d\n",iter);
    for (i=0;i<input1.rows();i++)
    {
      for (j=0;j<input1.cols();j++)
      {
        double up_min_h = 1e32,down_min_h = 1e32,left_min_h = 1e32,right_min_h = 1e32;
//        printf("\tLoop over p...\n");
        for (d=0;d<dmax;d++) //loops over p to get min_h
        {
          double tmp_d = dMatrix[i*input1.cols()+j][d]/2.0;
          up_h[d] = tmp_d +(get_left(right[d], i, j)+get_down(up[d], i, j)+get_right(left[d], i, j))/6.0;
          down_h[d] = tmp_d + (get_left(right[d], i, j)+get_up(down[d], i, j)+get_right(left[d], i, j))/6.0;
          left_h[d] = tmp_d + (get_down(up[d], i, j)+get_up(down[d], i, j)+get_right(left[d], i, j))/6.0;
          right_h[d] = tmp_d + (get_left(right[d], i, j)+get_down(up[d], i, j)+get_up(down[d], i, j))/6.0;
          if (up_h[d]<up_min_h)
          {
            up_min_h = up_h[d];
          }
          if (down_h[d]<down_min_h)
          {
            down_min_h = down_h[d];
          }
          if (left_h[d]<left_min_h)
          {
            left_min_h = left_h[d];
          }
          if (right_h[d]<right_min_h)
          {
            right_min_h = right_h[d];
          }
        }
//        printf("\tLoop over q...\n");
        for (k=0;k<dmax;k++)//loops over q for m
        {
          tmp_up[k][i][j] = min(up_h[k], up_min_h+dist);
          tmp_down[k][i][j] = min(down_h[k], down_min_h+dist);
          tmp_left[k][i][j] = min(left_h[k], left_min_h+dist);
          tmp_right[k][i][j] = min(right_h[k], right_min_h+dist);
        }
      }
    }
    for (i=0;i<input1.rows();i++)
    {
      for (j=0;j<input1.cols();j++)
      {
        for (k=0;k<dmax;k++)
        {
          up[k][i][j] = tmp_up[k][i][j];
          down[k][i][j] = tmp_down[k][i][j];
          left[k][i][j] = tmp_left[k][i][j];
          right[k][i][j] = tmp_right[k][i][j];
        }
      }
    }
  }
  for (i=0;i<input1.rows();i++)
    for (j=0;j<input1.cols();j++)
    {
      double min_cost=1e32;
      int min_d=-1;
      for (k=0;k<dmax;k++)
      {
        double tmp_cost = dMatrix[i*input1.cols()+j][k];
        tmp_cost += (get_up(down[k],i,j)+get_left(right[k], i, j)+get_down(up[k], i, j)+get_right(left[k], i, j))/4.0;
        if (tmp_cost<min_cost)
        {
          min_cost = tmp_cost;
          min_d = k;
        }
      }
      result[i][j] = min_d*3;
    }
  return result;
}

//convert a (i,j) label to index in storing matrix
int ij2index(int i, int j, int n)
{
  return (2*n-i+1)/2+j-i;
}

//convert storing matrix index back to a (i,j) label
void index2ij(int index, int n, int &i, int &j)
{
  i=0;
  j=0;
  int col_num=n;
  while (index - col_num >= 0)
  {
    index-=col_num;
    col_num--;
    i++;
  }
  j = index+i;
  return;
}

double calc_u(SDoublePlane image, int col, int i, int j, double top, double bottom)
{
  double cost = 0.0;
  double mean_top = 0.0, mean_bottom=0.0;
  double scale=10;
  for (int row=0;row<=i;row++)
  {
    cost += real_abs(image[row][col]-top)*scale;
    mean_top +=image[row][col];
  }
//  printf("i=%d,j=%d\ntop assign cost = %f, ",i,j,cost);
  mean_top/=(i+1);
  double sigma_top=0.0;
  for (int row=0;row<=i;row++)
  {
    sigma_top += (image[row][col]-mean_top)*(image[row][col]-mean_top);
  }
  for (int row=i+1; row<j;row++)
  {
    cost += 20.0*scale;
  }
  for (int row=j;row<image.rows();row++)
  {
    cost += real_abs(image[row][col]-bottom)*scale;
    mean_bottom +=image[row][col];
  }
  mean_bottom/=(image.rows()-j);
  double sigma_bottom=0.0;
  for (int row=j;row<image.rows();row++)
  {
    sigma_bottom += (image[row][col]-mean_bottom)*(image[row][col]-mean_bottom);
  }
//  printf("label cost = %f, sig_top = %f, sig_bottom = %f\n",cost, sigma_top, sigma_bottom);
//  getchar();
  return cost+sigma_top+sigma_bottom;
}

double calc_h(int rows, int i1, int j1, int i2, int j2)
{
  double cost = 0.0;
  double label_diff_1 = 200.0, label_diff_2 = 1000.0;
/*  for (int row = 0;row<rows;row++)
  {
    if(row<i1 && row>=i2 && row<j2)
      cost += label_diff_1;
    else if(row<i1 && row >=j2)
      cost += label_diff_2;
    else if(row>=i1 && row<j1 && row < i2)
      cost += label_diff_1;
    else if (row>=i1 && row<j1 && row >= j2)
      cost += label_diff_1;
    else if (row>=j1 && row<i2)
      cost += label_diff_2;
    else if (row>=j1 && row>=i2 && row<j2)
      cost += label_diff_1;
  }
*/
  if (i1<i2)
  {
    if(j1<i2)
    {
      cost = (j1-i1+j2-i2)*label_diff_1+(i2-j1)*label_diff_2;
    }
    else
    {
      if(j1<j2)
      {
        cost = (i2-i1+j2-j1)*label_diff_1;
      }
      else
      {
        cost = (i2-i1+j1-j2)*label_diff_1;
      }
    }
  }
  else
  {
    if(i1>j2)
    {
      cost = (i2-j2+j1-i1)*label_diff_1+(i1-j2)*label_diff_2;
    }
    else
    {
      if(j1<j2)
      {
        cost = (i1-i2+j2-j1)*label_diff_1;
      }
      else
      {
        cost = (i1-i2+j1-j2)*label_diff_1;
      }
    }
  }
//  printf("i1=%d,i2=%d,j1=%d,j2=%d, h smooth cost = %f\n",i1,i2,j1,j2,cost);
//  getchar();
  return cost;
}

void scene_labeling(const SDoublePlane &image, SDoublePlane &R, SDoublePlane &G, SDoublePlane &B)
{
  printf("Initialization...\n");
  printf("Resampling...\n");
  SDoublePlane input = resample(image, image.rows()/3, image.cols()/3);
  R = SDoublePlane(input.rows(),input.cols());
  G = SDoublePlane(input.rows(),input.cols());
  B = SDoublePlane(input.rows(),input.cols());
	int t_length = input.rows()/5+1;
	int *t1 = new int[t_length];
	int *t2 = new int[t_length];
	int *t3 = new int[t_length];
	int *t4 = new int[t_length];
	int *t5 = new int[t_length];
	int *t6 = new int[t_length];
	for (int i=0;i<t_length;i++)
	{
		t1[i]=t2[i]=t3[i]=t4[i]=t5[i]=t6[i]=-1;
	}
//  SDoublePlane result(input.rows(),input.cols());
  int num_label = (input.rows()+1)*input.rows()/2;
  SDoublePlane Uk(input.cols(),num_label);
  double l_top=0.0, l_bottom=0.0;
  int max_line = 5;
  for (int i=0;i<max_line;i++)
  {
    for (int j=0;j<input.cols();j++)
    {
      l_top+=input[i][j];
      l_bottom+=input[input.rows()-i-1][j];
    } 
  }
  l_top/=input.cols()*max_line;
  l_bottom/=input.cols()*max_line;
  printf("Calculating Uk...\n");
  for (int k=0;k<input.cols();k++)
  {
    int index = 0;
    for (int i=0;i<input.rows();i++)
    {
      for (int j=i;j<input.rows();j++)
      {
        Uk[k][index]=calc_u(input, k, i, j, l_top, l_bottom);
        index++;
      }
    }
    if(k%50==0)
      printf("Finished column %d\n",k);
  }
  printf("Calculating Hk...\n");
  SDoublePlane Hk(num_label,num_label);
  for (int i=0;i<num_label;i++)
  {
    for (int j=0;j<num_label;j++)
    {
      int i1,j1,i2,j2;
      index2ij(i, input.rows(),i1,j1);
      index2ij(j, input.rows(),i2,j2);
      Hk[i][j]=calc_h(input.rows(),i1,j1,i2,j2);
    }
    if(i%500==0)
      printf("Finished label %d\n",i);
  }
  //viterbi
  printf("Calculating Ek...\n");
  SDoublePlane Ek(input.cols(),num_label);
  SDoublePlane Path(input.cols(),num_label);
  for (int i=0;i<num_label;i++)
  {
    Ek[0][i] = Uk[0][i];
    Path[0][i] = -1;
  }
  
  for (int k=1; k< input.cols();k++)
  {
    for (int s = 0;s<num_label;s++)
    {
      double min_cost = 1e32;
      int min_prev_step = -1;
      for (int sp = 0;sp<num_label;sp++)
      {
        double tmp_cost = Ek[k-1][sp];
        tmp_cost += Hk[sp][s];
        if(tmp_cost<min_cost)
        {
          min_cost = tmp_cost;
          min_prev_step = sp;
        }
      }
      Ek[k][s]=Uk[k][s]+min_cost;
      Path[k][s]=min_prev_step;
    }
    if(k%50==0)
      printf("Finished column %d\n",k);
  }
  for (int i=0;i<input.rows();i++)
  {
    for (int j=0;j<input.cols();j++)
    {
      R[i][j] = input[i][j];
      G[i][j] = input[i][j];
      B[i][j] = input[i][j];
    }
  }

  int min_s = -1;
  double min_cost = 1e32;
  for (int s = 0;s<num_label;s++)
  {
    if (Ek[input.cols()-1][s]<min_cost)
    {
      min_cost=Ek[input.cols()-1][s];
      min_s = s;
    }
  }
  printf("Generating path...\n");
  for (int k=input.cols()-1;k>=0;k--)
  {
    int c1, c2;
    index2ij(min_s, input.rows(), c1, c2);
    R[c1][k]=255.0;
    G[c1][k]=0;
    B[c1][k]=0;
    R[c2][k]=0;
    B[c2][k]=0;
    G[c2][k]=255.0;
    min_s = (int)Path[k][min_s];
  }
  printf("Finished!\n");
  return;
}

SDoublePlane mrf_linear_stereo(const SDoublePlane &inputleft, const SDoublePlane &inputright)
{
  int w = 3;
  int dmax = 30;
  int i,j,k;
  int d=0;
  double dist = 2.5e4;
  double c = 1e3;
  SDoublePlane input1 = resample(inputleft, inputleft.rows()/2, inputleft.cols()/2);
  SDoublePlane input2 = resample(inputright, inputright.rows()/2, inputright.cols()/2);
  SDoublePlane result(input1.rows(), input1.cols());
  SDoublePlane dMatrix(input1.rows()*input1.cols(),dmax);
  int max_iter = input1.cols()+input1.rows();
  //up[i][j] means message from node (i,j) to (i-1,j)
  SDoublePlane up[dmax];
  SDoublePlane tmp_up[dmax];
  printf("start initializing up...\n");
  for (i=0;i<dmax;i++)
  {
    up[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_up[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  //down[i][j] means message from node (i,j) to (i+1,j)
  SDoublePlane down[dmax];
  SDoublePlane tmp_down[dmax];
  printf("start initializing down...\n");
  for (i=0;i<dmax;i++)
  {
    down[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_down[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  //left[i][j] means message from node (i,j) to (i,j-1)
  SDoublePlane left[dmax];
  SDoublePlane tmp_left[dmax];
  printf("start initializing left...\n");
  for (i=0;i<dmax;i++)
  {
    left[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_left[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  //right[i][j] means message from node (i,j) to (i,j+1)
  SDoublePlane right[dmax];
  SDoublePlane tmp_right[dmax];
  printf("start initializing right...\n");
  for (i=0;i<dmax;i++)
  {
    right[i] = SDoublePlane(input1.rows(), input1.cols());
    tmp_right[i] = SDoublePlane(input1.rows(), input1.cols());
  }
  printf("start initializing d-matrix...\n");
  for (i=0;i<input1.rows();i++)
    for (j=0;j<input1.cols();j++)
    {
	    for (d=0;d<dmax;d++)
	    {
		    double cost = 0;
		    int xmin = i-w>=0?i-w:0;
		    int xmax = i+w<input1.rows()?i+w:input1.rows()-1;
		    int ymin = j-w>=0?j-w:0;
		    int ymax = j+w<input1.cols()?j+w:input1.cols()-1;
		    for (int x = xmin; x<=xmax; x++)
		    {
			    for (int y = ymin; y<=ymax; y++)
			    {
				    if (y+d<input1.cols())
				    {
					    double l = input1[x][y];
					    double r = input2[x][y+d];
					    cost += (l-r)*(l-r);
				    }
				    else
				    {
					    double l = input1[x][y];
					    double r = 0;
					    cost += (l-r)*(l-r);
				    }
			    }
		    }
		    dMatrix[i*input1.cols()+j][d]=cost;
		    //printf("d-m[%d][%d]=%f\n",i,j,cost);
	    }
	  }

  double *up_h = new double[dmax], *down_h=new double[dmax], *left_h=new double[dmax],*right_h=new double[dmax];
//  for (int iter = 0;iter<dmax;iter++)
  printf("Start calculating...\n");
  for (int iter = 0;iter<max_iter;iter++)
  {
    if(iter%50==0)
      printf("Iteration #%d\n",iter);
    for (i=0;i<input1.rows();i++)
    {
      for (j=0;j<input1.cols();j++)
      {
        double up_min_h = 1e32,down_min_h = 1e32,left_min_h = 1e32,right_min_h = 1e32;
//        printf("\tLoop over p...\n");
        for (d=0;d<dmax;d++) //loops over p to get min_h
        {
          double tmp_d = dMatrix[i*input1.cols()+j][d]/2.0;
          up_h[d] = tmp_d +(get_left(right[d], i, j)+get_down(up[d], i, j)+get_right(left[d], i, j))/6.0;
          down_h[d] = tmp_d + (get_left(right[d], i, j)+get_up(down[d], i, j)+get_right(left[d], i, j))/6.0;
          left_h[d] = tmp_d + (get_down(up[d], i, j)+get_up(down[d], i, j)+get_right(left[d], i, j))/6.0;
          right_h[d] = tmp_d + (get_left(right[d], i, j)+get_down(up[d], i, j)+get_up(down[d], i, j))/6.0;
        }
        //foward step
        for (d=1;d<dmax;d++)
        {
          up_h[d] = min(up_h[d], up_h[d-1]+c);
          down_h[d] = min(down_h[d], down_h[d-1]+c);
          left_h[d] = min(left_h[d], left_h[d-1]+c);
          right_h[d] = min(right_h[d], right_h[d-1]+c);
        }
//        printf("\tLoop over q...\n");
        //backward step
        for (k=dmax-2;k>=0;k--)//loops over q for m
        {
          tmp_up[k][i][j] = min(up_h[k], up_h[k+1]+c);
          tmp_down[k][i][j] = min(down_h[k], down_h[k+1]+c);
          tmp_left[k][i][j] = min(left_h[k], left_h[k+1]+c);
          tmp_right[k][i][j] = min(right_h[k], right_h[k+1]+c);
        }
      }
    }
    for (i=0;i<input1.rows();i++)
    {
      for (j=0;j<input1.cols();j++)
      {
        for (k=0;k<dmax;k++)
        {
          up[k][i][j] = tmp_up[k][i][j];
          down[k][i][j] = tmp_down[k][i][j];
          left[k][i][j] = tmp_left[k][i][j];
          right[k][i][j] = tmp_right[k][i][j];
        }
      }
    }
  }
  for (i=0;i<input1.rows();i++)
    for (j=0;j<input1.cols();j++)
    {
      double min_cost=1e32;
      int min_d=-1;
      for (k=0;k<dmax;k++)
      {
        double tmp_cost = dMatrix[i*input1.cols()+j][k];
        tmp_cost += (get_up(down[k],i,j)+get_left(right[k], i, j)+get_down(up[k], i, j)+get_right(left[k], i, j))/4.0;
        if (tmp_cost<min_cost)
        {
          min_cost = tmp_cost;
          min_d = k;
        }
      }
      result[i][j] = min_d*3;
    }
  return result;
}

int main(int argc, char *argv[])
{
  if(argc != 4 && argc != 3)
    {
      cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file]" << endl;
      return 1;
    }

  string input_filename1 = argv[1], input_filename2 = argv[2];
  string gt_filename;
  if(argc == 4)
    gt_filename = argv[3];

  // read in images and gt
  SDoublePlane image1 = SImageIO::read_png_file(input_filename1.c_str());
  SDoublePlane image2 = SImageIO::read_png_file(input_filename2.c_str());
  SDoublePlane gt;
  if(gt_filename != "")
    gt = SImageIO::read_png_file(gt_filename.c_str());

/*
  // rotate image1 using nearest neighbor interpolation...
  const int angle_degrees = 72;
  SDoublePlane rotated_nn = rotate_nn(image1, angle_degrees/180.0*M_PI);
  SImageIO::write_png_file("rotated_nn.png", rotated_nn, rotated_nn, rotated_nn);

  // rotate image1 using bilinear interpolation...
  SDoublePlane rotated_bilin = rotate_bilinear(image1, angle_degrees/180.0*M_PI);
  SImageIO::write_png_file("rotated_bilin.png", rotated_bilin, rotated_bilin, rotated_bilin);


  // apply a projection matrix
  // as a placeholder, just use an identity matrix
  SDoubleMatrix proj_matrix(3,3);
  proj_matrix[0][0] = 0.907;
  proj_matrix[0][1] = 0.258;
  proj_matrix[0][2] = -182;
  proj_matrix[1][0] = -0.153;
  proj_matrix[1][1] = 1.44;
  proj_matrix[1][2] = 58;
  proj_matrix[2][0] = -0.000306;
  proj_matrix[2][1] = 0.000731;
  proj_matrix[2][2] = 1;

  SDoublePlane projected = projective_transformation(image1, proj_matrix);
  SImageIO::write_png_file("proj.png", projected, projected, projected);


  // do stereo using simple technique
  SDoubleMatrix disp = disparity_map(image1, image2);
  SImageIO::write_png_file("disp_simple.png", disp, disp, disp);

  // do stereo using hmm
  SDoubleMatrix disp2 = hmm_stereo(image1, image2);
  SImageIO::write_png_file("disp_hmm.png", disp2, disp2, disp2);
*/
  // do stereo using mrf
//  SDoubleMatrix disp3 = mrf_stereo(image1, image2);
//  SImageIO::write_png_file("disp_mrf.png", disp3, disp3, disp3);

//  SDoubleMatrix disp4 = mrf_linear_stereo(image1, image2);
//  SImageIO::write_png_file("disp_mrf_linear.png", disp4, disp4, disp4);
  // Measure error with respect to ground truth, if we have it...
  SDoublePlane r,g,b;
  scene_labeling(image1, r, g, b);
  SImageIO::write_png_file("scene_labeling.png", r, g, b);
  /*
  if(gt_filename != "")
    {
      double err=0;
      for(int i=0; i<gt.rows(); i++)
	for(int j=0; j<gt.cols(); j++)
	  err += sqrt((disp5[i][j] - gt[i][j])*(disp5[i][j] - gt[i][j]));

      cout << "Simple stereo technique mean error = " << err/gt.rows()/gt.cols() << endl;

      err=0;
      for(int i=0; i<gt.rows(); i++)
	for(int j=0; j<gt.cols(); j++)
	  err += sqrt((disp2[i][j] - gt[i][j])*(disp2[i][j] - gt[i][j]));

      cout << "HMM stereo technique mean error = " << err/gt.rows()/gt.cols() << endl;

    }
*/
  return 0;
}
