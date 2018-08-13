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

SDoublePlane mrf_stereo(const SDoublePlane &input1, const SDoublePlane &input2)
{
  int w = 3;
  int dmax = 30;
  int i,j,k;
  double tmpd[dmax];
  int max_iter = input1.cols()+input1.rows();
  SDoublePlane result(input1.rows(), input1.cols());
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

//  for (int iter = 0;iter<max_iter;iter++)
  for (int iter = 0;iter<dmax;iter++)
  {
    printf("Iteration #%d\n",iter);
    for (i=0;i<input1.rows();i++)
      for (j=0;j<input1.cols();j++)
      {
        getDisparity(i, j, dmax, w, tmpd, input1, input2);
        for (k=0;k<dmax;k++)//k loops over qj
        {
          double up_min_cost = 1e32,down_min_cost = 1e32,left_min_cost = 1e32,right_min_cost = 1e32;
          int up_min_d = -1,down_min_d = -1,left_min_d = -1,right_min_d = -1;
          for (int curr_d=0;curr_d<dmax;curr_d++)//cd loops over qi
          {
            double tmp_cost = 0;
            double up_cost = 0, down_cost=0, left_cost=0,right_cost=0;
            tmp_cost += tmpd[curr_d];
            tmp_cost += potts(k, curr_d);
            tmp_cost /= 2.0;
            up_cost = tmp_cost + (get_left(right[curr_d], i, j)+get_down(up[curr_d], i, j)+get_right(left[curr_d], i, j))/6.0;
            down_cost = tmp_cost + (get_left(right[curr_d], i, j)+get_up(down[curr_d], i, j)+get_right(left[curr_d], i, j))/6.0;
            left_cost = tmp_cost + (get_down(up[curr_d], i, j)+get_up(down[curr_d], i, j)+get_right(left[curr_d], i, j))/6.0;
            right_cost = tmp_cost + (get_left(right[curr_d], i, j)+get_down(up[curr_d], i, j)+get_up(down[curr_d], i, j))/6.0;
            if (up_cost<up_min_cost)
            {
              up_min_cost = up_cost;
              up_min_d = curr_d;
            }
            if (down_cost<down_min_cost)
            {
              down_min_cost = down_cost;
              down_min_d = curr_d;
            }
            if (left_cost<left_min_cost)
            {
              left_min_cost = left_cost;
              left_min_d = curr_d;
            }
            if (right_cost<right_min_cost)
            {
              right_min_cost = right_cost;
              right_min_d = curr_d;
            }
          }
          tmp_up[k][i][j] = up_min_cost;
          tmp_down[k][i][j] = down_min_cost;
          tmp_left[k][i][j] = left_min_cost;
          tmp_right[k][i][j] = right_min_cost;
        }
        /*
        tmp_up[i][j] = get_left(right, i, j)+get_down(up, i, j)+get_right(left, i, j);
        tmp_down[i][j] = get_left(right, i, j)+get_up(down, i, j)+get_right(left, i, j);
        tmp_left[i][j] = get_down(up, i, j)+get_up(down, i, j)+get_right(left, i, j);
        tmp_right[i][j] = get_left(right, i, j)+get_down(up, i, j)+get_up(down, i, j);
        */
      }
    for (i=0;i<input1.rows();i++)
      for (j=0;j<input1.cols();j++)
        for (k=0;k<dmax;k++)
        {
          up[k][i][j] = tmp_up[k][i][j];
          down[k][i][j] = tmp_down[k][i][j];
          left[k][i][j] = tmp_left[k][i][j];
          right[k][i][j] = tmp_right[k][i][j];
        }
  }
  for (i=0;i<input1.rows();i++)
    for (j=0;j<input1.cols();j++)
    {
      getDisparity(i, j, dmax, w, tmpd, input1, input2);
      double min_cost=1e32;
      int min_d=-1;
      for (k=0;k<dmax;k++)
      {
        double tmp_cost = tmpd[k];
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
  SDoubleMatrix disp3 = mrf_stereo(image1, image2);
  SImageIO::write_png_file("disp_mrf.png", disp3, disp3, disp3);

  // Measure error with respect to ground truth, if we have it...
  /*
  if(gt_filename != "")
    {
      double err=0;
      for(int i=0; i<gt.rows(); i++)
	for(int j=0; j<gt.cols(); j++)
	  err += sqrt((disp2[i][j] - gt[i][j])*(disp2[i][j] - gt[i][j]));

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
