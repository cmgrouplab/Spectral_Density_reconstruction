//********************************************************
//** 2D construction of spectral function **
//********************************************************


//author: duyu chen
//version: 04/16/2016


//modified: 08/11/2022
//focus on vector-based stealthy system
//modifications including: (1) new data structure for saving k-vectors for the cosntrained region
//                         (2) use vector form chi_vector
//                         (3) only compute chi_vector in a small square region that just enclose the target region
//                         (4) do not consider shape function of individual pixels, only focus on the collective coordinates of the binary matrix 

//Questions worth investigating: (a) morphology of materials for different symmetry of vector chi 
//                               (b) realiability w.r.t. volume fraction and X 
//                               (c) meaningful defintiion of X for different shaped regions


using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#define MAXX 120 
#define MAXY 21474 //this cannot be very large... 
#define Nt MAXX/2 //sample length
#define MAXS (MAXX*MAXX) //number of lattic points...
#define SD 2
int NP; //number of black pixels 

//*******************************************
//cooling schedule

double alpha = 0.85;
//double beta = 0.85;

int TN = 150;
double T = 0.000001;

int Nevl = 10000;
int Nmax = 25000;

//*******************************************


double f1;
double f2;

int xind; 
int yind;
int indexn;
//int config[MAXX][MAXX];
//int best_config[MAXX][MAXX];

int* x;
int* y;

int* bestx;
int* besty;

int sitex[MAXS];
int sitey[MAXS];


int config[MAXX][MAXX]; //using the underlying configuration 
int inverseconfig[MAXX][MAXX];





int flag_iconfig; //start from random config or read-in

double p_acc_c; //the critcial trial move acceptance rate,
//when p_acc_tmp lower than this, start the DPNs selecting rule
double p_acc_tmp = 1.0;  
//0 means start surface opt from the beginning
//1 means start at p_acc_c
int flag_surface_opt = 0; //the flag for DPNs rule
int flag_setup_DPNs = 0;

int num_DPNs_counter[MAXS];

vector< vector<int> > DPNs_B; //the DPNs...
vector< vector<int> > DPNs_W; //the DPNs...

int num_DPNs_B[9]; //the number of black pixels for each DPNs
int num_DPNs_W[9]; //for white pixels

double prob_DPNs_B[9];//the probablity 
double prob_DPNs_W[9]; 
//to compute, generate a rand num in [0, 1], then find the first prob_DPNs that smaller than rand_num

double prob_M = 0.5; //this is for prob for the pixel with maximum number of DP neighbors
int neighbor_max_B;
int neighbor_max_W;

int sum_B = 0;
int sum_W = 0;
int moved_B_num; //the number of DPN of the moved pixels
int moved_W_num;

int moved_B_it; //the position in the array
int moved_W_it;

//the following are for spectral density 

int ct_region = 0; //number of k vectors that are constrained
int vector_region[MAXS][SD]; //store the k vectors that define the constrained region
                    //due to symmetry requirement of spectral density, it should possess inversion symmetry

#define K_bound 12 //this the largest linear size of the region to be confined, corresponding to hard edge length of a square
              //this square should bound the constrained region


double pi = 3.14159265358979;
int Nk = K_bound; //not larger than MAXX/2, for stealthy can be chosen smaller
double Kspace = 2 * pi / MAXX;

double Kbin = 0.036; //these two are for angular averaged
#define Kbin_num 150


double J_k_Re[2*K_bound+1][2*K_bound+1]; //2Nk+1, make sure consistency, only compute useful k's
double J_k_Im[2*K_bound+1][2*K_bound+1];
double dJ_k_Re[2*K_bound+1][2*K_bound+1];
double dJ_k_Im[2*K_bound+1][2*K_bound+1];

//the vector-based Chi_k, which is simply the square of the mod of collective coordinates 
//due to the symmetry constrants, only need to consider positive part of k_y axis, and the postive k_x half plane
//we make the matrix larger, so that it is consistent with J_K, which is the collective coordinates, roughly half of the matrix is empty
double Chi_k_vector[2*K_bound+1][2*K_bound+1];
double Chi_k_vectorT[2*K_bound+1][2*K_bound+1];
double Chi_k_vectorb[2*K_bound+1][2*K_bound+1];

double Chi_k_vector_plot[2*K_bound+1][2*K_bound+1]; //this is for plot, we fill all of the entries.

//these are for angular averaged chi_K
/*
double K_Pos[Kbin_num];
int K_Counter[Kbin_num];
double obj[Kbin_num];
double Chi_k[Kbin_num];
double Chi_k_T[Kbin_num];
double Chi_k_b[Kbin_num];
*/
//vector form

double average_intensity;

ofstream fout;



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void read_parameter()
{
  cout<<"Reading parameters for annealing reconstruction from standard input"<<endl;

  cout<<"Init config flag flag_iconfig = "; cin>>flag_iconfig; cout<<flag_iconfig<<endl;

  cout<<"Number of black pixels: NP = "; cin>>NP; cout<<NP<<endl;

  //now computing the volume fraction...
  f1 = (double)NP/(double)(MAXX*MAXX);
  f2 = 1.0 - f1;
  average_intensity = f1;
  //now for the cooling schedule
  cout<<"Starting temp T0 = "; cin>>T; cout<<T<<endl;
  cout<<"Decreasing ratio: alpha = "; cin>>alpha; cout<<alpha<<endl;
  cout<<"Number of decreasing T stages: TN = "; cin>>TN; cout<<TN<<endl;
  cout<<"Number of pxiel move per stage: Nevl = "; cin>>Nevl; cout<<Nevl<<endl;
  cout<<"Maximum number of accepted pixel moves per stage: Nmax = "; cin>>Nmax; cout<<Nmax<<endl;
  cout<<"Critical acceptance rate for surface opt: p_acc_c = "; cin>>p_acc_c; cout<<p_acc_c<<endl;
  //alloc the particle positions and the separation distance matrix
  //the separation distance matrix is not necessary for superlarge systems..
  x = new int[NP];
  y = new int[NP];
  bestx = new int[NP];
  besty = new int[NP];
  cout << "init successful!" << endl;

  
}


void get_region_circle()
{
	cout<<"computing k vectors for the constrained CIRCLE region..."<<endl;
	cout<<"K_bound = "<<K_bound<<endl;
	
	double R_c;
	cout<<"R_circle ="; cin>>R_c;
	if(R_c>=K_bound)
	{
		cout<<"R_c > K_bound! Abort!"<<endl;
		exit(1);
	}
	
	//now need find the k vectors that contribute to the target function
	
	ct_region = 0; //number of k vectors that are constrained

    //int vector_region[MAXS][SD];
	
	//first deal with the half positive plane
	for(int i = 1; i <= Nk; i++)
		for(int j = - Nk; j <= Nk; j++)
		{
			double temp_dist = sqrt(i*i+j*j); //distance to origin
			
			if(temp_dist<=R_c)
				{
					vector_region[ct_region][0] = i;
					vector_region[ct_region][1] = j; //this can be negative, but converted to positie by a translation
					
					ct_region++;
				}
		}
		
	//now deal with the positive y portion
	for(int j=1; j<=Nk; j++)
	{
		if(j<=R_c)
		{
			vector_region[ct_region][0] = 0;
			vector_region[ct_region][1] = j; //this can be negative, but converted to positie by a translation
					
			ct_region++;
		}
	}
	
	cout<<"ct_region = "<<ct_region<<endl;
	
	//print out the k vector to test
	fout.open("k_region.xls");
	
	for(int i=0; i<ct_region; i++)
		fout<<vector_region[i][0]<<'\t'<<vector_region[i][1]<<endl;
	fout.close();
	
}


//this is for radial averaged chi_K
/*
void get_obj()
{
	FILE* fp;

	if((fp = fopen("chi_k_obj.txt","r"))==NULL)
    {
      printf("Can not open objective file for spectral function! Abort!\n");
      exit(1);
    }

  //double fr = 0.5;//(double)NP/(double)N;

  //printf("here\n");

	for(int k = 0; k < Kbin_num; k++)
	{
		fscanf(fp, "%lf", &obj[k]);
	}


	fclose(fp);
}
*/



double GetInnerProduct(double Vector1[SD], double Vector2[SD])
{
	double sum = 0;

	for(int i = 0; i < SD; i++)
		sum += Vector1[i] * Vector2[i];
	return sum;
}

//this is for radial averaged chi_k
/*
void Get_Chi_k(double chi_k[Kbin_num])
{
	double Chi_k_temp;
	int t;
	double k_dis;
	double KPoint[SD];
	for(int i = 0; i < Kbin_num; i++)
	{
		chi_k[i] = 0.0;
	}
	for(int i = 1; i <= Nk; i++)
		for(int j = - Nk; j <= Nk; j++)
		{
			KPoint[0] = i * Kspace;
			KPoint[1] = j * Kspace;

			k_dis = 0.0;
			for(int d = 0; d < SD; d++)
				k_dis += KPoint[d] * KPoint[d];
			k_dis = sqrt(k_dis);
			t = floor(k_dis / Kbin);
			if(t < Kbin_num)
			{
				Chi_k_temp = (J_k_Re[i+Nk][j+Nk] * J_k_Re[i+Nk][j+Nk] + J_k_Im[i+Nk][j+Nk] * J_k_Im[i+Nk][j+Nk]) / (double)MAXS; //chi_k from collective coordinates
				
				Chi_k_temp = Chi_k_temp * 4.0 * sin(KPoint[0] / 2.0) * sin(KPoint[0] / 2.0) / KPoint[0] / KPoint[0];//shape function for square pixel
				
				if(j != 0)
					Chi_k_temp = Chi_k_temp * 4.0 * sin(KPoint[1] / 2.0) * sin(KPoint[1] / 2.0) / KPoint[1] / KPoint[1];			
				
				chi_k[t] += Chi_k_temp; 
				//i+Nk, j+Nk
			}
		}

	for(int i = 1; i <= Nk; i++)
	{
		KPoint[0] = 0.0;
		KPoint[1] = i * Kspace;
		k_dis = KPoint[1];
		t = floor(k_dis / Kbin);

		if(t < Kbin_num)
		{
			Chi_k_temp = (J_k_Re[Nk][i+Nk] * J_k_Re[Nk][i+Nk] + J_k_Im[Nk][i+Nk] * J_k_Im[Nk][i+Nk]) / (double)MAXS;
			Chi_k_temp = Chi_k_temp * 4.0 * sin(KPoint[1] / 2.0) * sin(KPoint[1] / 2.0) / KPoint[1] / KPoint[1];			
			chi_k[t] += Chi_k_temp;
		}
	
	}

	for(t = 0; t < Kbin_num; t++)
		if(K_Counter[t] != 0)
		{
			chi_k[t] = chi_k[t] / K_Counter[t];
		}
}
*/


void Get_Chi_k_vector(double chi_k_vector[2*K_bound+1][2*K_bound+1])
{
	double Chi_k_temp;
	
	//int t;
	//double k_dis;
	//double KPoint[SD];
	
	/*
	for(int i = 0; i < Kbin_num; i++)
	{
		chi_k[i] = 0.0;
	}
	*/
	
	for(int i = 1; i <= Nk; i++)
		for(int j = - Nk; j <= Nk; j++)
		{
			/*
			KPoint[0] = i * Kspace;
			KPoint[1] = j * Kspace;

			k_dis = 0.0;
			for(int d = 0; d < SD; d++)
				k_dis += KPoint[d] * KPoint[d];
			k_dis = sqrt(k_dis);
			t = floor(k_dis / Kbin);
			if(t < Kbin_num)
			{
				Chi_k_temp = (J_k_Re[i+Nk][j+Nk] * J_k_Re[i+Nk][j+Nk] + J_k_Im[i+Nk][j+Nk] * J_k_Im[i+Nk][j+Nk]) / (double)MAXS; //chi_k from collective coordinates
				
				//Chi_k_temp = Chi_k_temp * 4.0 * sin(KPoint[0] / 2.0) * sin(KPoint[0] / 2.0) / KPoint[0] / KPoint[0];//shape function for square pixel
				
				
				//if(j != 0)
				//	Chi_k_temp = Chi_k_temp * 4.0 * sin(KPoint[1] / 2.0) * sin(KPoint[1] / 2.0) / KPoint[1] / KPoint[1];			
				
				//the above should also for shape function
				
				chi_k[t] += Chi_k_temp; 
				//i+Nk, j+Nk
			}
			*/
			
			Chi_k_temp = (J_k_Re[i+Nk][j+Nk] * J_k_Re[i+Nk][j+Nk] + J_k_Im[i+Nk][j+Nk] * J_k_Im[i+Nk][j+Nk]) / (double)MAXS; //chi_k from collective coordinates
			
			chi_k_vector[i+Nk][j+Nk] = Chi_k_temp;
		}

	for(int i = 1; i <= Nk; i++)
	{
		/*
		KPoint[0] = 0.0;
		KPoint[1] = i * Kspace;
		k_dis = KPoint[1];
		t = floor(k_dis / Kbin);

		if(t < Kbin_num)
		{
			Chi_k_temp = (J_k_Re[Nk][i+Nk] * J_k_Re[Nk][i+Nk] + J_k_Im[Nk][i+Nk] * J_k_Im[Nk][i+Nk]) / (double)MAXS;
			Chi_k_temp = Chi_k_temp * 4.0 * sin(KPoint[1] / 2.0) * sin(KPoint[1] / 2.0) / KPoint[1] / KPoint[1];			
			chi_k[t] += Chi_k_temp;
		}
		*/
		
		Chi_k_temp = (J_k_Re[Nk][i+Nk] * J_k_Re[Nk][i+Nk] + J_k_Im[Nk][i+Nk] * J_k_Im[Nk][i+Nk]) / (double)MAXS;
		chi_k_vector[Nk][i+Nk] = Chi_k_temp;
	}

    /*
	for(t = 0; t < Kbin_num; t++)
		if(K_Counter[t] != 0)
		{
			chi_k[t] = chi_k[t] / K_Counter[t];
		}
		*/
}

void Get_J_k()
{
	
	double KPoint[SD];
	
	for(int i = 1; i <= Nk; i++)
		for(int j = - Nk; j <= Nk; j++)
		{
			KPoint[0] = i * Kspace;
			KPoint[1] = j * Kspace;
			double sum_cos = 0.0;
			double sum_sin = 0.0;
			double temp;
			double temp_center[SD];
	
			for(int m = 0; m < MAXX; m++)
				for(int n = 0; n < MAXX; n++)
			{
				temp_center[0] = (double)m;
				temp_center[1] = (double)n;
				temp = GetInnerProduct(temp_center, KPoint); 
				sum_cos += cos(temp) * (config[m][n] - average_intensity);
				sum_sin += sin(temp) * (config[m][n] - average_intensity);
			}

			J_k_Re[i+Nk][j+Nk] = sum_cos;
			J_k_Im[i+Nk][j+Nk] = sum_sin;
			
		}

	for(int i = 1; i <= Nk; i++)
	{
		KPoint[0] = 0.0;
		KPoint[1] = i * Kspace;
		double sum_cos = 0.0;
		double sum_sin = 0.0;
		double temp;
		double temp_center[SD];
	
		for(int m = 0; m < MAXX; m++)
			for(int n = 0; n < MAXX; n++)
			{
				temp_center[0] = (double)m;
				temp_center[1] = (double)n;
				temp = GetInnerProduct(temp_center, KPoint); 
				sum_cos += cos(temp) * (config[m][n] - average_intensity);
				sum_sin += sin(temp) * (config[m][n] - average_intensity);
			}

		J_k_Re[Nk][i+Nk] = sum_cos;
		J_k_Im[Nk][i+Nk] = sum_sin;

	}

	/*
	for(int i = 0; i <= Nk - 1; i++)
		for(int j = 0; j <= 2 * Nk; j++)
		{
			J_k_Re[i][j] = J_k_Re[2*Nk-i][2*Nk-j];
			J_k_Im[i][j] = J_k_Im[2*Nk-i][2*Nk-j];
		}

	for(int i = 0; i <= Nk - 1; i++)
	{
		J_k_Re[Nk][i] = J_k_Re[Nk][2*Nk-i];
		J_k_Im[Nk][i] = J_k_Im[Nk][2*Nk-i];
	}

	J_k_Re[Nk][Nk] = 0.0;
	J_k_Im[Nk][Nk] = 0.0;
	*/
}

int touch(int ind)
{
  int i = x[ind];
  int j = y[ind];

  if(config[i][j] == 1) 
	  return 1;
  else 
	  return 0;
}


void init_data()
{
	for(int i = 0; i < MAXX; i++)
		for(int j = 0; j < MAXX; j++)
		{
			sitex[MAXX*i+j] = i;
			sitey[MAXX*i+j] = j;
		}

	for(int i = 0; i <= 2 * Nk; i++)
		for(int j = 0; j <= 2 * Nk; j++)
		{
			J_k_Re[i][j] = 0.0;
			J_k_Im[i][j] = 0.0;
			dJ_k_Re[i][j] = 0.0;
			dJ_k_Im[i][j] = 0.0;
		}

   /*
	int t;
	double k_dis;
	double KPoint[2];
	for(int i = 0; i < Kbin_num; i++)
	{
		K_Pos[i] = 0.0;
		K_Counter[i] = 0;
	}
	for(int i = 1; i <= Nk; i++)
		for(int j = - Nk; j <= Nk; j++)
		{
			
			KPoint[0] = i * Kspace;
			KPoint[1] = j * Kspace;

			k_dis = 0.0;
			for(int d = 0; d < SD; d++)
				k_dis += KPoint[d] * KPoint[d];
			k_dis = sqrt(k_dis);
			t = floor(k_dis / Kbin);
			if(t < Kbin_num)
			{
				K_Pos[t] += k_dis;
				K_Counter[t] ++;
			}
		}

	for(int i = 1; i <= Nk; i++)
	{
		KPoint[0] = 0.0;
		KPoint[1] = i * Kspace;
		k_dis = KPoint[1];
		t = floor(k_dis / Kbin);

		if(t < Kbin_num)
		{			
			K_Pos[t] += k_dis;
			K_Counter[t] ++;
		}
	
	}

	for(t = 0; t < Kbin_num; t++)
		if(K_Counter[t] != 0)
		{
			K_Pos[t] = K_Pos[t] / K_Counter[t];
		}
	*/	
		
	Get_J_k();
}


void init_config()
{
	for(int i = 0; i < MAXX; i++)
		for(int j = 0; j < MAXX; j++)
      {
		config[i][j] = 0;
		inverseconfig[i][j] = -1;
      }
  //initalize the underlying configuration


	for(int i = 0; i < NP; i++)
    {
      do{
        x[i] = rand()% MAXX;
        y[i] = rand()% MAXX;
      }while(touch(i)==1);

      config[x[i]][y[i]] = 1;
	  inverseconfig[x[i]][y[i]] = i;
    }
  
 
	init_data();
}


void read_config()
{
	for(int i = 0; i < MAXX; i++)
		for(int j = 0; j < MAXX; j++)
		{
			config[i][j] = 0;
			inverseconfig[i][j] = -1;
		}


	FILE* fp;

	if((fp=fopen("Iconfig.txt","r"))==NULL)
    {
      printf("No Iconfig.txt is found! Abort!\n");
      exit(1);
    }

	int tempx;
	int tempy;

	for(int i = 0; i < NP; i++)
    {

      fscanf(fp, "%d", &tempx);
      fscanf(fp, "%d", &tempy);

      config[tempx][tempy] = 1;
	  x[i] = tempx;
	  y[i] = tempy;
	  inverseconfig[x[i]][y[i]] = i;
    }

	fclose(fp);

	init_data();
}



void setup_DPNs()
{
  flag_setup_DPNs = 1;
  //cout << "successfully set up DPN!" << endl;
  for(int i = 0; i < 9; i++)
    {
      //cout<<"here"<<endl;
      DPNs_B.push_back(vector<int>(1, -1));
      DPNs_W.push_back(vector<int>(1, -1));

      num_DPNs_B[i] = 0; //the number of black pixels for each DPNs
      num_DPNs_W[i] = 0; //for white pixels

      prob_DPNs_B[i] = 0.0; 
      prob_DPNs_W[i] = 0.0; 
    }

  //now loop over config to find DPNs for the first time
  int temp_num;
  int index1, index2;
  int p_index;
  
  for(int i = 0; i < MAXX; i++)
    for(int j = 0; j < MAXX; j++)
	{
	  temp_num = 0;

	  p_index = MAXX * i + j;

	 if(config[i][j] == 1)
	{
	      for(int m = -1; m <= 1; m++)
			for(int n = -1; n <= 1; n++)
		    {
		      index1 = i + m; 
		      if(index1<0) 
				  index1 = index1 + MAXX;
		      else 
				  if(index1>=MAXX) 
					  index1 = index1 - MAXX;

		      index2 = j + n; 
		      if(index2<0) 
				  index2 = index2 + MAXX;
		      else 
				  if(index2>=MAXX) 
					  index2 = index2 - MAXX;

		      if(config[index1][index2] != config[i][j])
				temp_num++;
		    }

	      DPNs_B[temp_num].push_back(p_index);
	      num_DPNs_B[temp_num]++;
	      num_DPNs_counter[p_index] = temp_num;
	      
	 }

	 else
	 {
	      for(int m = -1; m <= 1; m++)
			for(int n = -1; n <= 1; n++)
		    {
		      index1 = i + m; 
		      if(index1<0) index1 = index1 + MAXX;
		      else if(index1>=MAXX) index1 = index1 - MAXX;
		      
		      index2 = j + n; 
		      if(index2<0) 
				  index2 = index2 + MAXX;
		      else 
				  if(index2>=MAXX) 
					  index2 = index2 - MAXX;

		      if(config[index1][index2] != config[i][j])
				temp_num++;
		    }

	      DPNs_W[temp_num].push_back(p_index);
	      num_DPNs_W[temp_num]++;
	      num_DPNs_counter[p_index] = temp_num;
	 }
	}

  //now we find the probability for each pixel...

	for(int i = 8; i >= 0; i--)
    {
		if(num_DPNs_B[i] > 0)
		{
			neighbor_max_B = i;
			break;
		}
    }
	for(int i = 8; i >= 0; i--)
    {
		if(num_DPNs_W[i] > 0)
		{
			neighbor_max_W = i;
			break;
		}
    }


   prob_DPNs_B[neighbor_max_B] = prob_M;
   prob_DPNs_W[neighbor_max_W] = prob_M;

   

   double sum_PB =0; 
   double sum_PW = 0;
   
   for(int i = 0; i < neighbor_max_B; i++)
     sum_B += num_DPNs_B[i] * (i+1) * (i+1);

   for(int i = 0; i < neighbor_max_B; i++)
     {
       prob_DPNs_B[i] = (1 - prob_M) * num_DPNs_B[i] * (i+1) * (i+1) / sum_B;
       sum_PB += prob_DPNs_B[i];
     }


   for(int i=0; i<neighbor_max_W; i++)
     sum_W += num_DPNs_W[i] * (i+1) * (i+1);

   for(int i=0; i<neighbor_max_W; i++)
     {
       prob_DPNs_W[i] = (1 - prob_M) * num_DPNs_W[i] * (i+1) * (i+1) / sum_W;
       sum_PW += prob_DPNs_W[i];
     }
   
}


void change_config()
{
	if(flag_surface_opt == 0)
    {
		int s;
		do{
		int ind = rand()% NP;

		int mod = MAXX;
		int rx = rand()% mod;
		int ry = rand()% mod;

		int Tx = x[ind];
		int Ty = y[ind];

		x[ind] = x[ind] + rx;
		if(x[ind] >= MAXX) 
			x[ind] = x[ind] - MAXX;
		else 
			if(x[ind] < 0) 
				x[ind] = x[ind] + MAXX;
		y[ind] = y[ind] + ry;
		if(y[ind] >= MAXX) 
			y[ind] = y[ind] - MAXX;
		else 
			if(y[ind] < 0) 
				y[ind] = y[ind] + MAXX;

		 if(touch(ind) == 1)
		{
			x[ind] = Tx;
			y[ind] = Ty;

			s = 0;
		}
		 else     
		{
			indexn = ind;
			xind = Tx;
			yind = Ty;
			
			//change the underlying configuration
			config[x[indexn]][y[indexn]] = 1;
			config[xind][yind] = 0;
			inverseconfig[x[indexn]][y[indexn]] = indexn; // new posistion
			inverseconfig[xind][yind] = -1; //old posistion
			s = 1;
		}
		}while(s == 0);
	}

	else 
    {
      //this is the first time to setup DPNs 
		if(flag_setup_DPNs == 0)
		{	
			cout<<"set up DPNS"<<endl;

			setup_DPNs();

		}
	  //cout<<"max_W = "<<neighbor_max_W<<endl;
	  //cout<<"max_B = "<<neighbor_max_B<<endl;

	  //first select which group of pixels should be moved
	  double rand_B = (double)(rand()%MAXY)/(double)MAXY;
	  double rand_W = (double)(rand()%MAXY)/(double)MAXY;

	  //by default, we choose the pixel with the largest DPNs
	  moved_B_num = neighbor_max_B;
	  moved_W_num = neighbor_max_W;
	  double temp_PB = 0.0;
	  double temp_PW = 0.0;

	  for(int i = neighbor_max_B; i >= 0; i--)
	    {
		  temp_PB += prob_DPNs_B[i];
	      if(temp_PB > rand_B)
		{
		  moved_B_num = i;

		  break;
		}
	    }

	  for(int i = neighbor_max_W; i >= 0; i--)
	    {
		  temp_PW += prob_DPNs_W[i];
	      if(temp_PW > rand_W)
		{
		  moved_W_num = i;
		  
		  break;
		}
	    }

	  //now from the group, select a pixel to be moved
	  int size_B = num_DPNs_B[moved_B_num];
	  int size_W = num_DPNs_W[moved_W_num];

	  
	  moved_B_it = rand() % (size_B) + 1;
	  moved_W_it = rand() % (size_W) + 1;
	  
	  int temp_ind_B = DPNs_B[moved_B_num][moved_B_it];
	  int temp_ind_W = DPNs_W[moved_W_num][moved_W_it];

	 
	  

	  xind = sitex[temp_ind_B];
	  yind = sitey[temp_ind_B];
	  

	  indexn = inverseconfig[xind][yind];
	  x[indexn] = sitex[temp_ind_W];
	  y[indexn] = sitey[temp_ind_W];
	  inverseconfig[x[indexn]][y[indexn]] = indexn; 
	  config[x[indexn]][y[indexn]] = 1; // new posistion


	  config[xind][yind] = 0;
      inverseconfig[xind][yind] = -1; //old posistion
	 
	}
}

void update_DPNs()
{
  //cout<<"updating DPNS"<<endl;
	
	//int temp_index = 0;
	int index1, index2;
	DPNs_W[moved_W_num].erase(DPNs_W[moved_W_num].begin()+moved_W_it);//first update these when the order is still correct
	DPNs_B[moved_B_num].erase(DPNs_B[moved_B_num].begin()+moved_B_it);
	num_DPNs_B[moved_B_num]--;
	num_DPNs_W[moved_W_num]--;
	//first, check the black pixel
	int new_B_num = 0;

 
  
	if(config[x[indexn]][y[indexn]] != 1)
    {
      cout<<"the moved black pixel is not right! check again!"<<endl;
      exit(1);
    }

	//temp_index = DPNs_B[moved_B_num][moved_B_it];
  
	int temp_index_B;
	int temp_index_W;
  

	for(int i = -1; i <= 1; i++)
		for(int j = -1; j <= 1; j++)
	{
	  index1 = i + x[indexn]; 
	  if(index1<0) index1 = index1 + MAXX;
	  else if(index1>=MAXX) index1 = index1 - MAXX;
	  
	  index2 = j + y[indexn]; 
	  if(index2<0) index2 = index2 + MAXX;
	  else if(index2>=MAXX) index2 = index2 - MAXX;
	  
	  

	  temp_index_B = index1 * MAXX + index2;
	  if(config[index1][index2] != config[x[indexn]][y[indexn]])
	    {
			new_B_num++;
			if((index1 != xind) || (index2 != yind))
			{
			for(int v = num_DPNs_W[num_DPNs_counter[temp_index_B]]; v >= 1; v--)
			{
				if(temp_index_B == DPNs_W[num_DPNs_counter[temp_index_B]][v])
				{
					DPNs_W[num_DPNs_counter[temp_index_B]].erase(DPNs_W[num_DPNs_counter[temp_index_B]].begin()+v);
					break;
				}
			}
			
			DPNs_W[num_DPNs_counter[temp_index_B]+1].push_back(temp_index_B);

			num_DPNs_W[num_DPNs_counter[temp_index_B]]--;
			num_DPNs_W[num_DPNs_counter[temp_index_B]+1]++;
			num_DPNs_counter[temp_index_B]++;
			}
	    }
	  else
		  if(i != 0 || j != 0)
		{
			for(int v = num_DPNs_B[num_DPNs_counter[temp_index_B]]; v >= 1; v--)
			{
				if(temp_index_B == DPNs_B[num_DPNs_counter[temp_index_B]][v])
				{
					DPNs_B[num_DPNs_counter[temp_index_B]].erase(DPNs_B[num_DPNs_counter[temp_index_B]].begin()+v);
					break;
				}
			}
			
		  DPNs_B[num_DPNs_counter[temp_index_B]-1].push_back(temp_index_B);
		  num_DPNs_B[num_DPNs_counter[temp_index_B]]--;
		  num_DPNs_B[num_DPNs_counter[temp_index_B]-1]++;

		  num_DPNs_counter[temp_index_B]--;
		}
	}

  
  
	DPNs_B[new_B_num].push_back(x[indexn] * MAXX + y[indexn]);
    
    num_DPNs_B[new_B_num]++;
    
	num_DPNs_counter[x[indexn] * MAXX + y[indexn]] = new_B_num;

  //now the white pxiel

  //first, check the black pixel
	int new_W_num = 0;

	if(config[xind][yind] != 0)
    {
      cout<<"the moved white pixel is not right! check again!"<<endl;
      exit(1);
    }

	//temp_index = DPNs_W[moved_W_num][moved_W_it];
 
  

  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++)
	{
	  index1 = i + xind; 
	  if(index1<0) index1 = index1 + MAXX;
	  else if(index1>=MAXX) index1 = index1 - MAXX;
	  
	  index2 = j + yind; 
	  if(index2<0) index2 = index2 + MAXX;
	  else if(index2>=MAXX) index2 = index2 - MAXX;
	  
	  temp_index_W = index1 * MAXX + index2;

	  if(config[index1][index2] != config[xind][yind])
	  {
	    new_W_num++;
		if((index1 != x[indexn]) || (index2 != y[indexn]))
		{
		for(int v = num_DPNs_B[num_DPNs_counter[temp_index_W]]; v >= 1; v--)
			{
				if(temp_index_W == DPNs_B[num_DPNs_counter[temp_index_W]][v])
				{
					DPNs_B[num_DPNs_counter[temp_index_W]].erase(DPNs_B[num_DPNs_counter[temp_index_W]].begin()+v);
					break;
				}
			}
		DPNs_B[num_DPNs_counter[temp_index_W]+1].push_back(temp_index_W);
		num_DPNs_B[num_DPNs_counter[temp_index_W]]--;
		num_DPNs_B[num_DPNs_counter[temp_index_W]+1]++;

		num_DPNs_counter[temp_index_W]++;
		}
	  }
	  else
		 if(i != 0 || j != 0)
		{
			for(int v = num_DPNs_W[num_DPNs_counter[temp_index_W]]; v >= 1; v--)
			{
				if(temp_index_W == DPNs_W[num_DPNs_counter[temp_index_W]][v])
				{
					DPNs_W[num_DPNs_counter[temp_index_W]].erase(DPNs_W[num_DPNs_counter[temp_index_W]].begin()+v);
					break;
				}
			}
			DPNs_W[num_DPNs_counter[temp_index_W]-1].push_back(temp_index_W);

			num_DPNs_W[num_DPNs_counter[temp_index_W]]--;
			num_DPNs_W[num_DPNs_counter[temp_index_W]-1]++;
			num_DPNs_counter[temp_index_W]--;
		}
	}
  
  DPNs_W[new_W_num].push_back(xind * MAXX + yind);
 
 
  num_DPNs_W[new_W_num]++;
    
  
  num_DPNs_counter[xind * MAXX + yind] = new_W_num;
  //cout<<"finsihing updating"<<endl;
  
  //***************************************************************
  //now we update the probability for each pixel...
  for(int i = 0; i < 9; i++)
    {
       prob_DPNs_B[i] = 0.0; 
       prob_DPNs_W[i] = 0.0; 
    }

  for(int i = 8; i >= 0; i--)
    {
      if(num_DPNs_B[i] > 0)
	{
	  neighbor_max_B = i;
	  break;
	}
    }
   for(int i = 8; i >= 0; i--)
    {
      if(num_DPNs_W[i]>0)
	{
	  neighbor_max_W = i;
	  break;
	}
    }

   prob_DPNs_B[neighbor_max_B] = prob_M;
   prob_DPNs_W[neighbor_max_W] = prob_M;

   sum_B = 0;
   sum_W = 0;

   double sum_PB =0; 
   double sum_PW = 0;
   
   for(int i = 0; i < neighbor_max_B; i++)
     sum_B += num_DPNs_B[i] * (i+1) * (i+1);
   for(int i = 0; i < neighbor_max_B; i++)
     {
       prob_DPNs_B[i] = (1 - prob_M) * num_DPNs_B[i] * (i+1) * (i+1) / sum_B;
	   sum_PB += prob_DPNs_B[i];
     }


   for(int i=0; i<neighbor_max_W; i++)
     sum_W += num_DPNs_W[i] * (i+1) * (i+1);
   for(int i=0; i<neighbor_max_W; i++)
     {
       prob_DPNs_W[i] = (1-prob_M)* num_DPNs_W[i] * (i+1) * (i+1) / sum_W;
       sum_PW += prob_DPNs_W[i];
     }
   //cout << "temp_index: " << temp_index << endl;
}

void resume_config()
{
  
	config[x[indexn]][y[indexn]] = 0;
	config[xind][yind] = 1;

	inverseconfig[x[indexn]][y[indexn]] = -1;
	inverseconfig[xind][yind] = indexn;
  
	x[indexn] = xind;
	y[indexn] = yind;
	//retain the underlying configuration
}

//this is specifically for the stealthy system
double energy(double Chi_k_vector[2*K_bound+1][2*K_bound+1])
{

  double E = 0;

  //get_obj(); //this decides which structure we will generate
  
  for(int i = 0; i < ct_region; i++)
    {
      E = E + Chi_k_vector[vector_region[i][0]+Nk][vector_region[i][1]+Nk] * Chi_k_vector[vector_region[i][0]+Nk][vector_region[i][1]+Nk];
    }
  
  return E;
  
}


double d_energy(double Chi_k_vector[2*K_bound+1][2*K_bound+1], double Chi_k_vectorT[2*K_bound+1][2*K_bound+1])
{
  double d_E = 0;
  d_E = energy(Chi_k_vectorT) - energy(Chi_k_vector);
  return d_E;
}

double PE(double dE, double T) // the probability that should compare ...
{
  if(dE > 0) 
	  return exp(-dE/T);
  else 
	  return 1;
}

void free_mem()
{
	delete [] x;
	delete [] y;
	delete [] bestx;
	delete [] besty;
}

int main()
{
	srand(time(NULL));

	double energyt = 0;
	double energyb = 250000;

	//get_obj();

	get_region_circle();
	//exit(1);

	FILE * fp = fopen("Energy.txt","w");
	fclose(fp);


	

	read_parameter();
  
	if(flag_iconfig == 0)
		init_config();//initialize configuration. volume fraction preserved...
	else
		read_config();

  //init_config();

	Get_Chi_k_vector(Chi_k_vector);
	
	double energyi = energy(Chi_k_vector);
	
	
	if(energyi < energyb)
	{
	    energyb = energyi;
		for(int it = 0; it < NP; it++)
		{
			bestx[it] = x[it];
			besty[it] = y[it];
		}
		
		/*
		for(int n = 0; n < Kbin_num; n++)
		{
				Chi_k_b[n] = Chi_k[n];
		}     
		*/
		
		for(int i = -Nk; i <= Nk; i++)
			for(int j=-Nk; j<=Nk; j++)
			{	
				Chi_k_vectorb[i+Nk][j+Nk] = Chi_k_vector[i+Nk][j+Nk]; //update the full data set
			} 
	}
	
	
	cout << "initial energy: " << energyi << endl;
	
	//we simply print out this as matrix, without worring about the actual k vector values
	fp = fopen("TChi_k.txt","w");
	for(int i = -Nk; i <=Nk ; i++)
	{
		for(int j=-Nk; j<=Nk; j++)
			fprintf(fp, "%1.15f\t", Chi_k_vector[i+Nk][j+Nk]);
			
		fprintf(fp, "\n");
	}
	  
	fclose(fp);

  //printf("SN[0] = %d\t BN[0] = %d \n", SN[0], BN[0]);



  //simulated annealing procedure to evlove the system
  //**********************************************************************
  //*********************************************************************

	for(int q = 0; q < TN; q++)
    {
		T = alpha * T;
		if(p_acc_tmp < p_acc_c) //acceptance probablity computed last time
			flag_surface_opt = 1;
		int Nacc = 0;
		int Nmv = 0;
		
		for(int in = 0; in < Nevl; in++)
		{
	  
		
		change_config();
		Nmv++;
		
        double Center_old[SD];
		double Center_new[SD];
		Center_old[0] = xind;
		Center_old[1] = yind;
		Center_new[0] = x[indexn];
		Center_new[1] = y[indexn];
		for(int i = 1; i <= Nk; i++)
			for(int j = -Nk; j <= Nk; j++)
			{
				double KPoint[SD];
				KPoint[0] = i * Kspace;
				KPoint[1] = j * Kspace;
				dJ_k_Re[i+Nk][j+Nk] = cos(GetInnerProduct(Center_new, KPoint)) - cos(GetInnerProduct(Center_old, KPoint));
				dJ_k_Im[i+Nk][j+Nk] = sin(GetInnerProduct(Center_new, KPoint)) - sin(GetInnerProduct(Center_old, KPoint));

				J_k_Re[i+Nk][j+Nk] = J_k_Re[i+Nk][j+Nk] + dJ_k_Re[i+Nk][j+Nk];
				J_k_Im[i+Nk][j+Nk] = J_k_Im[i+Nk][j+Nk] + dJ_k_Im[i+Nk][j+Nk];
			}

		for(int i = 1; i <= Nk; i++)
		{
			double KPoint[SD];
			KPoint[0] = 0.0;
			KPoint[1] = i * Kspace;
			dJ_k_Re[Nk][i+Nk] = cos(GetInnerProduct(Center_new, KPoint)) - cos(GetInnerProduct(Center_old, KPoint));
			dJ_k_Im[Nk][i+Nk] = sin(GetInnerProduct(Center_new, KPoint)) - sin(GetInnerProduct(Center_old, KPoint));

			J_k_Re[Nk][i+Nk] = J_k_Re[Nk][i+Nk] + dJ_k_Re[Nk][i+Nk];
			J_k_Im[Nk][i+Nk] = J_k_Im[Nk][i+Nk] + dJ_k_Im[Nk][i+Nk];
		}
	  

		Get_Chi_k_vector(Chi_k_vectorT);

	  //Now we use Monte Carlo Steps

		double P = (double)(rand()%MAXY)/(double)MAXY;

		if( P > PE(d_energy(Chi_k_vector, Chi_k_vectorT), T))
	    {
			resume_config();
	      
			for(int i = 1; i <= Nk; i++)
				for(int j = -Nk; j <= Nk; j++)
			{
				J_k_Re[i+Nk][j+Nk] = J_k_Re[i+Nk][j+Nk] - dJ_k_Re[i+Nk][j+Nk];
				J_k_Im[i+Nk][j+Nk] = J_k_Im[i+Nk][j+Nk] - dJ_k_Im[i+Nk][j+Nk];
			}

			for(int i = 1; i <= Nk; i++)
			{	
			
				J_k_Re[Nk][i+Nk] = J_k_Re[Nk][i+Nk] - dJ_k_Re[Nk][i+Nk];
				J_k_Im[Nk][i+Nk] = J_k_Im[Nk][i+Nk] - dJ_k_Im[Nk][i+Nk];
			}
	      
		}
	      
		else 
	    {
			for(int i = -Nk; i <= Nk; i++)
			for(int j=-Nk; j<=Nk; j++)
			{	
				Chi_k_vector[i+Nk][j+Nk] = Chi_k_vectorT[i+Nk][j+Nk]; //update the full data set
			}
			if(flag_surface_opt == 1)
				update_DPNs();
			Nacc++; //acceptance rate
	    }
	  
	  
	  //compare and record the best energy and configuration...
		energyt = energy(Chi_k_vector);
	  
		if(energyt < energyb)
	    {
			energyb = energyt;
	      
	      
			for(int it = 0; it < NP; it++)
			{
				bestx[it] = x[it];
				besty[it] = y[it];
			}
		
 		/*
			for(int it = 0; it < Kbin_num; it++)
			{
				Chi_k_b[it] = Chi_k[it];
			}
		*/
			
		for(int i = -Nk; i <= Nk; i++)
			for(int j=-Nk; j<=Nk; j++)
			{	
				Chi_k_vectorb[i+Nk][j+Nk] = Chi_k_vector[i+Nk][j+Nk]; //update the full data set
			} 
	      
	    }
	    
	  //printf("%f   %d change has finished... \n", energyt, i+1 );
	  
		if(Nmv > Nmax)
			break;
	  

		}
		p_acc_tmp = (double)Nacc/(double)Nmv;
		printf("%d th change of temperature has finished...: T = %lf\n", q+1, T);
		cout<<"The acceptance rate: "<<p_acc_tmp<<endl;
		cout<<"The energy E = "<<energyb<<endl;
       

		fp = fopen("Energy.txt","a");
		fprintf(fp, "%1.12f \n", energyb);
		fclose(fp);

		fp = fopen("Mconfig.txt","w");
		for(int it=0; it<NP; it++)
		{
			fprintf(fp, "%d \t %d \n", x[it], y[it]);
		 }
		fclose(fp);

		fp = fopen("M_bconfig.xls","w");
		for(int it=0; it<NP; it++)
		{
			fprintf(fp, "%d \t %d \n", bestx[it], besty[it]);
		 }
		fclose(fp);

       
		fp = fopen("Chi_k_b.xls", "w");
		for(int i = -Nk; i <=Nk ; i++)
	  	{
		for(int j=-Nk; j<=Nk; j++)
			fprintf(fp, "%1.15f\t", Chi_k_vectorb[i+Nk][j+Nk]);
			
		fprintf(fp, "\n");
    	}
		fclose(fp);
	

		fp = fopen("Chi_k.xls", "w");
		for(int i = -Nk; i <=Nk ; i++)
		{
		for(int j=-Nk; j<=Nk; j++)
			fprintf(fp, "%1.15f\t", Chi_k_vector[i+Nk][j+Nk]);
			
		fprintf(fp, "\n");
		}
		fclose(fp);
  
		printf("*************************************************\n");
       
	}
  //********************************************************************
  //*******************************************************************
  //this is the end of simulated annealing...

	fp = fopen("Fconfig.txt","w");
	for(int it=0; it<NP; it++)
    {
	    fprintf(fp, "%d \t %d \n", x[it], y[it]);
    }
	fclose(fp);
  
 
	free_mem();
	return 0;

}



