#include "SImage.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <limits>

#ifdef STANDALONE
#include <SImageIO.h>
#endif

//#define STANDALONE

//#define INFINITY numeric_limits<double>::infinity()

using namespace std;
typedef vector< pair<int, int> > PathType;
const char *basename_g = 0;
double top_smoothness_g = 1e3, bottom_smoothness_g = 1e3;
double top_peak_g = 0.5, bottom_peak_g = 0.5;
double repulse_g = 10;

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

inline double potts(int x, int y)
{
  return fabs(x-y)<=3?0:10000;
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
/*
  double calc_u(int col, int i, int j, SDoublePlane &gradients)
  {
  double cost = 20000.0;
  double min_cost=1000.0;
  double scale = 1000.0;	
  cost -= gradients[i][col]+gradients[j][col];
  //cost += 50000/(j-i+1);
  double vert_penal = 10000.0-fabs((j-i)*scale);
  if (vert_penal<min_cost)
  {
  vert_penal = min_cost;
  }
  cost += vert_penal;
  return cost;
  }
*/

double calc_u(int col, int i, int j, SDoublePlane &gradients, double cost, double gradScale)
{
  //double cost = 20000.0;
  cost -= gradients[i][col]+gradients[j][col];
  //cost += 50000/(j-i+1);
  cost += gradScale/(j-i+1);
  return cost;
}

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
double calc_h(int rows, int i1, int j1, int i2, int j2, double label_diff_1, double label_diff_2)
{
  double cost = 0.0;
  //double label_diff_1 = 200.0, label_diff_2 = 1000.0;
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

SDoublePlane calc_uk(const SDoublePlane &input, int num_label, double cost, double gradScale, SDoublePlane &gradients)
{  
  int rows = input.rows();
  SDoublePlane Uk(input.cols(), num_label);
  printf("Calculating Uk...\n");
  gradients=SDoublePlane(input.rows(),input.cols());
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
	      gradient+=fabs(input[j-1][k]-input[j][k])*scale;
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


  return Uk;
}

inline double sqr(double x) { return x * x; }


// dt with quadratic distance
void dt(const double *src, double *dst, double *dst_ind, int s1, int s2, int d1, int d2, double scale, int off=0)
{
  int d = (d1+d2) >> 1;
  int s = s1;
  for (int p = s1; p <= s2; p++)
    if (src[s] + sqr(s-d-off) * scale> src[p] + sqr(p-d-off) * scale)
      s = p;
  dst[d] = src[s] + sqr(s-d-off) * scale;
  dst_ind[d] = s;
  
  if(d-1 >= d1)
    dt(src, dst, dst_ind, s1, s, d1, d-1, scale, off);
  if(d2>=d+1)
    dt(src, dst, dst_ind, s, s2, d+1, d2, scale, off);
}


void dt_1d(const double *f, double scale, double *result, double *dst_ind, int beg, int end, int off=0) 
{
  dt(f, result, dst_ind, beg, end-1, beg, end-1, scale, off);
}

PathType find_single_path(const SDoublePlane &unary, double lambda)
{
  //  SDoublePlane Ek1(unary.rows(), unary.cols()), path(unary.rows(), unary.cols());
  SDoublePlane Ek1(unary.cols(), unary.rows()), path(unary.cols(), unary.rows());
  for (int i=0;i<unary.rows();i++)
    Ek1[0][i] = (unary[i][0]);
  
  /*
  for (int k=1; k<unary.cols(); k++)
    {
      for(int new_i=0; new_i < unary.rows(); new_i++)
	{
	  int min_i = -1;
	  double min_cost = INFINITY;
	  double my_unary = unary[new_i][k];
	  for(int old_i=0; old_i < unary.rows(); old_i++)
	    {
	      double cost = lambda * sqr(old_i-new_i);
	      cost += Ek1[old_i][k-1] + my_unary;
	      if(cost < min_cost)
  		min_i = old_i, min_cost = cost;
	    }
	  Ek1[new_i][k] = min_cost;
	  path[new_i][k] = min_i;
	}
    }
  */

  // use DT
  //  D(j) = v(j) + min_i d(i) + v(i,j)
  for (int k=1; k<unary.cols(); k++)
    {
      //      double unary_old[unary.rows()], dt_result[unary.rows()];
      //      int dt_result_ind[unary.rows()];
      //      for(int i=0; i<unary.rows(); i++)
      //	unary_old[i] = Ek1[i][k-1];
      double *Ek1_k = Ek1[k];

      dt_1d(Ek1[k-1], lambda, Ek1_k, path[k], 0, unary.rows()-1);

      for(int i=0; i<unary.rows(); i++)
	Ek1_k[i] += unary[i][k];
      //      for(int i=0; i<unary.rows(); i++)
      //	{
      //	  Ek1[k][i] = dt_result[i] + unary[i][k];
      //	  path[k] = dt_result_ind[i];
      //	}
    }

  int min_s = -1;
  double min_cost = INFINITY;
  double *last_ek1 = Ek1[unary.cols()-1];
  for (int s = 0 ; s<unary.rows();s++)
    {
      if (last_ek1[s]<min_cost)
	{
	  min_cost=last_ek1[s];
	  min_s = s;
	}
    }

  cout << "Cost of path = " << min_cost << endl;
  cout << "Final state = " << min_s << endl;
  //  min_s=173*2;

  PathType Path(unary.cols());
  for (int k=unary.cols()-1; k>=0; k--)
    {
      int c1, c2;
      //      index2ij(min_s, rows, c1, c2);
      Path[k].first = min_s; Path[k].second = min_s;
      min_s = (int)path[k][min_s];
    }

  return Path;
}

void write_plane(const SDoublePlane &unary_top, const char *fname, const char *bn)
{
  return;

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
#ifdef STANDALONE
  SImageIO::write_png_file((string(bn) + "-" + fname).c_str(), g, g, g);
#endif
}

inline void set_both_to_max(double &A, double &B)
{
  if(A > B) B=A; else A=B;
}


PathType scene_labeling_1(const SDoublePlane &input, double cost, double gradeScale, const PathType &pts1, const PathType &pts2)
{
  cerr << "computing unary terms..." << endl;
  int rows = input.rows();
  int num_label = (input.rows()+1)*input.rows()/2;
  
  SDoublePlane gradients(input.rows(), input.cols());

    const int kk=5;
  /*  for(int i=0; i<input.rows(); i++)
    for(int j=0; j<input.cols(); j++)
      if(i >= kk && i < input.rows()-kk && j >= kk && j < input.cols()-kk)
	 //	gradients[i][j] = fabs(input[i+kk][j]-input[i-kk][j]) + fabs(input[i][j+kk] - input[i][j-kk]);
	 	 gradients[i][j] = sqrt(sqr(input[i+kk][j]-input[i-kk][j]) + sqr(input[i][j+kk] - input[i][j-kk]));
	 //	 gradients[i][j] = max(0.0, input[i-kk][j]-input[i+kk][j]); // + fabs(input[i][j+kk] - input[i][j-kk]);
      else
	gradients[i][j] = 0;
  */
  
    
  for(int i=0; i<input.rows(); i++)
    for(int j=0; j<input.cols(); j++)
      if(i >= kk && i < input.rows()-kk) 
	{
	  double acc = 0;
	  for(int k=-kk; k<=-1; k++)
	    acc += sqr(255-input[i+k][j]);
	  for(int k=1; k<=kk; k++)
	    acc += sqr(input[i+k][j]);
	  //	  gradients[i][j] = -sqr(acc); 
	  gradients[i][j] = -(acc); 
	}
      else
	gradients[i][j] = -INFINITY;

  /*	 
  const int kk=2;
  for(int i=0; i<input.rows(); i++)
    for(int j=0; j<input.cols(); j++)
      if(i >= kk && i < input.rows()-kk) 
	{
	  double acc = 0;
	  for(int k=-kk; k<=-1; k++)
	    acc += sqr(input[i+k][j]);
	  for(int k=1; k<=kk; k++)
	    acc += sqr(input[i+k][j]);
	  //	    acc += sqr(input[i+k][j]);
	  gradients[i][j] = -(acc); 
	}
      else
	gradients[i][j] = -INFINITY;
  */

  cerr << "1" << endl;
  // for gradients, higher values are good
  
  SDoublePlane cumulative_fromtop_grad = gradients, cumulative_frombottom_grad = gradients;
  for(int j=0; j<input.cols(); j++)
    {
      double top_acc=-INFINITY, bottom_acc=-INFINITY;
      for(int i=0; i<input.rows(); i++)
	set_both_to_max(cumulative_fromtop_grad[i][j], top_acc);

      for(int i=input.rows()-1; i >= 0; i--)
	set_both_to_max(cumulative_fromtop_grad[i][j], bottom_acc);
    }

  // for cumulative gradients, higher values are bad
  cerr << "2" << endl;
  SDoublePlane unary_bottom(input.rows(), input.cols());
  SDoublePlane unary_top(input.rows(), input.cols());
  for(int i=0; i<input.rows(); i++)
    for(int j=0; j<input.cols(); j++)
      {
	if(i < input.rows() - 5*2) 
	  unary_bottom[i][j] = -gradients[i][j] + cumulative_frombottom_grad[i+5][j] * bottom_peak_g;
	else
	  unary_bottom[i][j] = INFINITY;

	if(i >= 5*2) 
	  unary_top[i][j] = -gradients[i][j] + cumulative_fromtop_grad[i-5][j] * top_peak_g;
	else
	  unary_top[i][j] = INFINITY;

      }
  cerr << "3" << endl;
  

  if(pts2.size() > 0)
    {
      for(int i=0; i<pts2.size(); i++)
	unary_bottom[pts2[i].second][pts2[i].first] = -1e10;
    }


  if(pts1.size() > 0)
    {
      for(int i=0; i<pts1.size(); i++)
	unary_top[pts1[i].second][pts1[i].first] = -1e10;
    }

  write_plane(unary_top, "unary_top.png", basename_g); 
  write_plane(unary_bottom, "unary_bottom.png", basename_g); 
  write_plane(gradients, "grad.png", basename_g); 
  write_plane(cumulative_frombottom_grad, "cum_bottom.png", basename_g);
  write_plane(cumulative_fromtop_grad, "cum_top.png", basename_g);



  printf("Doing inference...\n");
  //  SDoublePlane Ek(input.cols(),num_label);
  //  SDoublePlane Path(input.cols(),num_label);

  PathType top_path = find_single_path(unary_top, top_smoothness_g);

  // make bottom path be far away from top path
  for(int i=0; i<input.cols(); i++)
    for(int j=0; j<top_path[i].first+repulse_g; j++)
      if(j < unary_bottom.rows())
  	unary_bottom[j][i] = INFINITY;

  PathType bottom_path = find_single_path(unary_bottom, bottom_smoothness_g);
  PathType final_path(input.cols());
  for(int j=0; j<input.cols(); j++)
    final_path[j] = make_pair(top_path[j].first, bottom_path[j].first);

  return final_path;
}


PathType scene_labeling(const SDoublePlane &input, double cost, double gradeScale)
{
  int rows = input.rows();
  double label_diff_1 = 300.0, label_diff_2 = 2000.0;
  int num_label = (input.rows()+1)*input.rows()/2;

  SDoublePlane grad;
  SDoublePlane Uk = calc_uk(input, num_label, cost, gradeScale, grad);

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
      printf("calculating t=%d %d\n", k, input.cols());
      int **min_ip_t1 = new int*[rows];

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

  PathType path(input.cols());
  for (int k=input.cols()-1; k>=0; k--)
    {
      int c1, c2;
      index2ij(min_s, rows, c1, c2);
      path[k].first = c1; path[k].second = c2;
      min_s = (int)Path[k][min_s];
    }

  return path;
}


#ifdef STANDALONE
int main(int argc, char *argv[])
{
  if(argc != 5)
    {
      cerr << "usage: " << argv[0] << " input_image_file output_image_file, cost, gradScale" << endl;
      return 1;
    }

  printf("Reading image...\n");
  string input_filename = argv[1], output_filename = argv[2];
  basename_g = output_filename.c_str();
  double cost = atof(argv[3]), gradScale = atof(argv[4]);

  SDoublePlane input = SImageIO::read_png_file(input_filename.c_str());

  PathType path = scene_labeling_1(input, cost, gradScale, PathType(), PathType());

  printf("Writing output...\n");
  SDoublePlane R, G, B;
  R = input;  G = input;  B = input;
  for (int k=input.cols()-1; k>=0; k--)
    {
      int c1 = path[k].first, c2=path[k].second;
      R[c1][k]=255.0;      G[c1][k]=0;      B[c1][k]=0;
      R[c2][k]=0;      B[c2][k]=0;      G[c2][k]=255.0;
    }

  printf("png...\n");
  SImageIO::write_png_file(output_filename.c_str(), R, G, B);
  return 0;
}
#else

#include "mex.h" 
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {

  if(nlhs != 2 || nrhs < 3)
    {
      cerr << nlhs << " " << nrhs << endl;
      mexErrMsgTxt("input or output variable problem");
    }

  if(mxGetM(prhs[0]) != 1 ||  mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("first param must be a scalar");

  double *input = mxGetPr(prhs[1]);
  int rows= mxGetM(prhs[1]), cols = mxGetN(prhs[1]);
  cerr << rows << " " << cols << endl;
  SDoublePlane P(rows, cols);
  cerr << "init" << endl;
  for(int j=0, cp=0; j<cols; j++)
    for(int i=0; i<rows; i++, cp++)
      {
	P[i][j] = input[cp];
      }

  // get parameters
  double *params = mxGetPr(prhs[2]);
  if(mxGetM(prhs[2])*mxGetN(prhs[2]) != 5)
    mexErrMsgTxt("need two params");

  top_smoothness_g = params[0];
  bottom_smoothness_g = params[1];
  top_peak_g = params[2];
  bottom_peak_g = params[3];
  repulse_g = params[4];

  PathType pts1, pts2;
  if(nrhs > 3)
    {
      int m = mxGetN(prhs[3]);
      double *t = mxGetPr(prhs[3]);
      for(int i=0; i<m; i++)
	//	if(t[i*2] != -1)
	  pts1.push_back( pair<int, int>((int)t[i*2], (int)t[i*2+1]) );

      m = mxGetN(prhs[4]);
      t = mxGetPr(prhs[4]);
      for(int i=0; i<m; i++)
	//	if(t[i*2] != -1)
	  pts2.push_back( pair<int, int>((int)t[i*2], (int)t[i*2+1]) );

      cerr << pts1.size() << " " << pts2.size() << endl;
    }

  cerr << "scene labeling" << endl;
  PathType path = scene_labeling_1(P, 1, 1, pts1, pts2);

  /*  printf("Writing output...\n");
  SDoublePlane R, G, B;
  R = P;  G = P;  B = P;
  for (int k=P.cols()-1; k>=0; k--)
    {
      int c1 = path[k].first, c2=path[k].second;
      R[c1][k]=255.0;      G[c1][k]=0;      B[c1][k]=0;
      R[c2][k]=0;      B[c2][k]=0;      G[c2][k]=255.0;
    }
  */
  //  printf("png...\n");
  //  SImageIO::write_png_file("aa.png", R, G, B);
  
  mwSize dims[] = {rows, cols, 3};
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  /*  double *P_ptr = mxGetPr(plhs[0]);
  int cp=0;
  for(int j=0; j<cols; j++)
    for(int i=0; i<rows; i++, cp++)
      P_ptr[cp] = R[i][j];
  for(int j=0; j<cols; j++)
    for(int i=0; i<rows; i++, cp++)
      P_ptr[cp] = G[i][j];
  for(int j=0; j<cols; j++)
    for(int i=0; i<rows; i++, cp++)
      P_ptr[cp] = B[i][j];
  */

  plhs[1] = mxCreateDoubleMatrix(2, cols, mxREAL);
  double *P_ptr = mxGetPr(plhs[1]);
  for(int i=0; i<cols; i++)
    P_ptr[i*2] = path[i].first, P_ptr[i*2+1] = path[i].second;
  
} 

#endif
