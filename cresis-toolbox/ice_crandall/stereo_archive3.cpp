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

inline double real_abs(double x)
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
// an SDoubleMatrix is exactly the same as an SDoublePlane; we'll just give it an
//  alias to make it clear when we're talking about images and when we're talking
//  about matrices.
//
typedef SDoublePlane SDoubleMatrix;

// Apply a given 3x3 projective transformation matrix (aka a homography) to an image
//  
double abs(double x)
{
	return x<0?-x:x;
}

double potts(int x, int y)
{
	return abs(x-y)<=3?0:10000;
}


double min(double a, double b)
{
  return a>b?b:a;
}

//convert a (i,j) label to index in storing matrix
int ij2index(int i, int j, int n)
{
  return (2*n-i+1)*i/2+j-i;
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

double calc_u(int col, int i, int j, SDoublePlane &gradients, double min_cost, double scale)
{
  double cost = 20000.0;
  cost -= gradients[i][col]+gradients[j][col];
  //cost += 50000/(j-i+1);
  double vert_penal = 10000.0-real_abs((j-i)*scale);
  if (vert_penal<min_cost)
  {
    vert_penal = min_cost;
  }
  cost += vert_penal;
  return cost;
}
/*
double calc_u(int col, int i, int j, SDoublePlane &gradients, double cost, double gradScale)
{
  //double cost = 20000.0;
  cost -= gradients[i][col]+gradients[j][col];
  //cost += 50000/(j-i+1);
  cost += gradScale/(j-i+1);
  return cost;
}
*/
/*
double calc_u(int col, int i, int j, SDoublePlane &gradients, SDoublePlane &cumulativeGradients)
{
  double cost = 20000.0;
  double middle_scale=20.0;
  cost -= gradients[i][col]+gradients[j][col];
  cost += 50000/(j-i+1);
  cost += cumulativeGradients[j][col];
  return cost;
}
*/
double calc_h(int rows, int i1, int j1, int i2, int j2)
{
  double cost = 0.0;
  double label_diff_1 = 200.0, label_diff_2 = 1000.0;
  if (i1<i2)
  {
    if(j1<i2)
    {
      cost = (j1-i1+j2-i2)*label_diff_1+(i2-j1)*label_diff_2;
    }
    else
    {
      if(j1<j2)
      {//t=1
        cost = (i2-i1+j2-j1)*label_diff_1;
      }
      else
      {//t=2
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
      {//t=3
        cost = (i1-i2+j2-j1)*label_diff_1;
      }
      else
      {//t=4
        cost = (i1-i2+j1-j2)*label_diff_1;
      }
    }
  }
//  printf("i1=%d,i2=%d,j1=%d,j2=%d, h smooth cost = %f\n",i1,i2,j1,j2,cost);
//  getchar();
  return cost;
}

void scene_labeling(const SDoublePlane &image, SDoublePlane &R, SDoublePlane &G, SDoublePlane &B, double cost, double gradScale)
{
  printf("Initialization...\n");
  printf("Resampling...\n");
  SDoublePlane input = resample(image, image.rows()/3, image.cols()/3);
//  SDoublePlane input = image;
  R = SDoublePlane(input.rows(),input.cols());
  G = SDoublePlane(input.rows(),input.cols());
  B = SDoublePlane(input.rows(),input.cols());
  int rows = input.rows();
  double label_diff_1 = 200.0, label_diff_2 = 1000.0;
  int num_label = (input.rows()+1)*input.rows()/2;
  SDoublePlane Uk(input.cols(),num_label);
  printf("Calculating Uk...\n");
  SDoublePlane gradients(input.rows(),input.cols());
  //SDoublePlane cumulativeGradients(input.rows(),input.cols());
  double scale=100;
  int grade = 1;
  for (int i=0;i<input.rows();i++)
  {
    int low_i = i;
    int high_i = i+grade<rows?i+grade:rows-1;
    for (int k=0;k<input.cols();k++)
    {
      double gradient = 0.0;
      for (int j=low_i+1;j<=high_i;j++)
      {
        gradient+=abs(input[j-1][k]-input[j][k])*scale;
      }
      gradients[i][k]=gradient;
      /*
      if (i==0)
        cumulativeGradients[i][k]=gradient;
      else
        cumulativeGradients[i][k]=gradient+cumulativeGradients[i-1][k];
        */
    }
  }
  /*
  for (int i=0;i<input.rows();i++)
  {
    for (int k=0;k<input.cols();k++)
    {
      cumulativeGradients[i][k]=cumulativeGradients[input.rows()-1][k]-cumulativeGradients[i][k];
      //printf("row=%d,col=%d,grad=%f\n",i,k,cumulativeGradients[i][k]);
    }
//    getchar();
  }*/
  int index = 0;
  for (int i=0;i<input.rows();i++)
  {
    for (int j=i;j<input.rows();j++)
    {
      for (int k=0;k<input.cols();k++)
      {
        //Uk[k][index]=calc_u(k, i, j, gradients, cumulativeGradients);
        Uk[k][index]=calc_u(k, i, j, gradients, cost, gradScale);
      }
      index++;
    }
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
    //t=1, i1<=i2<=j1<=j2
  double **Ep1 = new double*[rows];
  double **F1 = new double*[rows];
  //t=2, i1<=i2<=j2<=i1
  double **Ep2 = new double*[rows];
  double **F2 = new double*[rows];
  //t=3, i2<=i1<=j1<=j2
  double **Ep3 = new double*[rows];
  double **F3 = new double*[rows];
  //t=4, i2<=i1<=j2<=j1
  double **Ep4 = new double*[rows];
  double **F4 = new double*[rows];
  for (int k=1; k< input.cols();k++)
  {
//    printf("calculating t=1\n");
    int **min_ip_t1 = new int*[rows];
//    printf("calculating Ep1\n");
    for (int jp=0;jp<rows;jp++)// for each jp
		{
		  double *g = new double[jp+1];
			Ep1[jp] = new double[jp+1];
			min_ip_t1[jp] = new int[jp+1];
			for (int i = 0;i<=jp;i++)
			{
			  g[i] = Ek[k-1][ij2index(i,jp,rows)]-i*label_diff_1;
				if (i==0)
				{
				  Ep1[jp][i] = g[i];
  			  min_ip_t1[jp][i] = 0;
				}
				else if(g[i]<Ep1[jp][i-1])
				{
				  Ep1[jp][i] = g[i];
				  min_ip_t1[jp][i]=i;
				}
				else
				{
				  Ep1[jp][i]=Ep1[jp][i-1];
				  min_ip_t1[jp][i]=min_ip_t1[jp][i-1];
				}
			}
			delete[] g;
	  }
//	  printf("calculating F1\n");
	  int **min_jp_t1 = new int*[rows];
	  for (int i=0;i<rows;i++)
	  {
	    min_jp_t1[i] = new int[rows-i+1];
	    double *g = new double[rows-i+1];
	    F1[i] = new double[rows-i+1];
	    double *h = new double[rows-i+1];
	    for (int j=i;j<rows;j++)
	    {
	      g[j-i]=-j*label_diff_1+Ep1[j][i];
	      if(j-i==0)
	      {
	        h[0]=g[0];
	        min_jp_t1[i][j-i] = j;
	      }
	      else if(g[j-i]<h[j-i-1])
	      {
	        h[j-i]=g[j-i];
	        min_jp_t1[i][j-i] = j;
	      }
	      else
	      {
	        h[j-i] = h[j-i-1];
	        min_jp_t1[i][j-i] = min_jp_t1[i][j-i-1];
	      }
	      F1[i][j-i] = (i+j)*label_diff_1+h[j-i];
	    }
	    delete[] g;
	    delete[] h;
	  }
	  
//    printf("calculating t=2\n");
    int **min_ip_t2 = new int*[rows];
//    printf("calculating Ep2\n");
    for (int jp=0;jp<rows;jp++)// for each jp
		{
		  double *g = new double[jp+1];
			Ep2[jp] = new double[jp+1];
			min_ip_t2[jp] = new int[jp+1];
			for (int i = 0;i<=jp;i++)
			{
			  g[i] = Ek[k-1][ij2index(i,jp,rows)]-i*label_diff_1;
				if (i==0)
				{
				  Ep2[jp][i] = g[i];
  			  min_ip_t2[jp][i] = 0;
				}
				else if(g[i]<Ep2[jp][i-1])
				{
				  Ep2[jp][i] = g[i];
				  min_ip_t2[jp][i]=i;
				}
				else
				{
				  Ep2[jp][i]=Ep2[jp][i-1];
				  min_ip_t2[jp][i]=min_ip_t2[jp][i-1];
				}
			}
			delete[] g;
	  }
//	  printf("calculating F2\n");
	  int **min_jp_t2 = new int*[rows];
	  for (int i=0;i<rows;i++)
	  {
	    min_jp_t2[i] = new int[rows-i+1];
	    double *g = new double[rows-i+1];
	    F2[i] = new double[rows-i+1];
	    double *h = new double[rows-i+1];
	    for (int j=rows-1;j>=i;j--)
	    {
	      g[j-i]=j*label_diff_1+Ep2[j][i];
	      if(j==rows-1)
	      {
	        h[j-i]=g[j-i];
	        min_jp_t2[i][j-i] = j;
	      }
	      else if(g[j-i]<h[j-i+1])
	      {
	        h[j-i]=g[j-i];
	        min_jp_t2[i][j-i] = j;
	      }
	      else
	      {
	        h[j-i] = h[j-i+1];
	        min_jp_t2[i][j-i] = min_jp_t2[i][j-i+1];
	      }
	      F2[i][j-i] = (i-j)*label_diff_1+h[j-i];
	    }
	    delete[] g;
	    delete[] h;
	  }
	  
//	  printf("calculating t=3\n");
    int **min_ip_t3 = new int*[rows];
//    printf("calculating Ep3\n");
    for (int jp=0;jp<rows;jp++)// for each jp
		{
		  double *g = new double[jp+1];
			Ep3[jp] = new double[jp+1];
			min_ip_t3[jp] = new int[jp+1];
			for (int i = jp;i>=0;i--)
			{
			  g[i] = Ek[k-1][ij2index(i,jp,rows)]+i*label_diff_1;
				if (i==jp)
				{
				  Ep3[jp][i] = g[i];
  			  min_ip_t3[jp][i] = jp;
				}
				else if(g[i]<Ep3[jp][i+1])
				{
				  Ep3[jp][i] = g[i];
				  min_ip_t3[jp][i]=i;
				}
				else
				{
				  Ep3[jp][i]=Ep3[jp][i+1];
				  min_ip_t3[jp][i]=min_ip_t3[jp][i+1];
				}
			}
			delete[] g;
	  }
//	  printf("calculating F3\n");
	  int **min_jp_t3 = new int*[rows];
	  for (int i=0;i<rows;i++)
	  {
	    min_jp_t3[i] = new int[rows-i+1];
	    double *g = new double[rows-i+1];
	    F3[i] = new double[rows-i+1];
	    double *h = new double[rows-i+1];
	    for (int j=i;j<rows;j++)
	    {
	      g[j-i]=-j*label_diff_1+Ep3[j][i];
	      if(j-i==0)
	      {
	        h[0]=g[0];
	        min_jp_t3[i][j-i] = j;
	      }
	      else if(g[j-i]<h[j-i-1])
	      {
	        h[j-i]=g[j-i];
	        min_jp_t3[i][j-i] = j;
	      }
	      else
	      {
	        h[j-i] = h[j-i-1];
	        min_jp_t3[i][j-i] = min_jp_t3[i][j-i-1];
	      }
	      F3[i][j-i] = (j-i)*label_diff_1+h[j-i];
	    }
	    delete[] g;
	    delete[] h;
	  }
	  
//	  printf("calculating t=4\n");
    int **min_ip_t4 = new int*[rows];
//    printf("calculating Ep4\n");
    for (int jp=0;jp<rows;jp++)// for each jp
		{
		  double *g = new double[jp+1];
			Ep4[jp] = new double[jp+1];
			min_ip_t4[jp] = new int[jp+1];
			for (int i = jp;i>=0;i--)
			{
			  g[i] = Ek[k-1][ij2index(i,jp,rows)]+i*label_diff_1;
				if (i==jp)
				{
				  Ep4[jp][i] = g[i];
  			  min_ip_t4[jp][i] = jp;
				}
				else if(g[i]<Ep4[jp][i+1])
				{
				  Ep4[jp][i] = g[i];
				  min_ip_t4[jp][i]=i;
				}
				else
				{
				  Ep4[jp][i]=Ep4[jp][i+1];
				  min_ip_t4[jp][i]=min_ip_t4[jp][i+1];
				}
			}
			delete[] g;
	  }
//	  printf("calculating F4\n");
	  int **min_jp_t4 = new int*[rows];
	  for (int i=0;i<rows;i++)
	  {
	    min_jp_t4[i] = new int[rows-i+1];
	    double *g = new double[rows-i+1];
	    F4[i] = new double[rows-i+1];
	    double *h = new double[rows-i+1];
	    for (int j=rows-1;j>=i;j--)
	    {
	      g[j-i]=j*label_diff_1+Ep4[j][i];
	      if(j==rows-1)
	      {
	        h[j-i]=g[j-i];
	        min_jp_t4[i][j-i] = j;
	      }
	      else if(g[j-i]<h[j-i+1])
	      {
	        h[j-i]=g[j-i];
	        min_jp_t4[i][j-i] = j;
	      }
	      else
	      {
	        h[j-i] = h[j-i+1];
	        min_jp_t4[i][j-i] = min_jp_t4[i][j-i+1];
	      }
	      F4[i][j-i] = (-i-j)*label_diff_1+h[j-i];
	    }
	    delete[] g;
	    delete[] h;
	  }
	  
//	  printf("calculating minimum ip jp\n");
	  for (int index=0;index<num_label;index++)
	  {	    
	    int i,j, ip, jp;
	    index2ij(index,rows,i,j);
	    int min_t = 1;
	    double min_F = F1[i][j-i];
	    if (F2[i][j-i]<min_F)
	    {
	      min_F = F2[i][j-i];
	      min_t = 2;
	    }
	    if (F3[i][j-i]<min_F)
	    {
	      min_F = F3[i][j-i];
	      min_t = 3;
	    }
	    if (F4[i][j-i]<min_F)
	    {
	      min_F = F4[i][j-i];
	      min_t = 4;
	    }
	    switch(min_t)
	    {
	    case 1:
	      jp = min_jp_t1[i][j-i];
	      ip = min_ip_t1[jp][i];
	      break;
	    case 2:
	    	jp = min_jp_t2[i][j-i];
	      ip = min_ip_t2[jp][i];
	      break;
	    case 3:
	    	jp = min_jp_t3[i][j-i];
	      ip = min_ip_t3[jp][i];
	      break;
	    case 4:
	    	jp = min_jp_t4[i][j-i];
	      ip = min_ip_t4[jp][i];
	      break;
	    }
	    Ek[k][index]=Uk[k][index]+min_F;
	    Path[k][index]=ij2index(ip,jp,rows);
    }
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
    index2ij(min_s, rows, c1, c2);
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

int main(int argc, char *argv[])
{
  if(argc != 5)
    {
      cerr << "usage: " << argv[0] << " input_image_file output_image_file, cost, gradScale" << endl;
      return 1;
    }

  string input_filename = argv[1], output_filename = argv[2];
  double cost = atof(argv[3]), gradScale = atof(argv[4]);
  // read in images and gt
  SDoublePlane input_image = SImageIO::read_png_file(input_filename.c_str());

  SDoublePlane r,g,b;
  scene_labeling(input_image, r, g, b, cost, gradScale);
  SImageIO::write_png_file(output_filename.c_str(), r, g, b);
  return 0;
}
