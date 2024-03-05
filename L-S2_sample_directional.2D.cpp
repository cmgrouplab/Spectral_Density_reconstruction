//*******************************************************
// S2+L (re)construction using orthogonal sampling method
//*******************************************************

//updated May, 28th, 2008

//modified: 05/14/2012
//the prevous method sampling L was not complete, missing 
// important events...

//version 2.0:
//use the most efficient method for sampling S2

//modified: 03/26/19
//to seprately output correlation functions on 2 orthogonal directions

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>


#define MAXX 256 //this has to be an even number !!! since Nt = MAXX/2
#define MAXY 25000
#define Nt 120

int NP;
double f1;//volume fraction of black pixels
double f2;

//the following is to treat pix as particles
int pix_position[MAXX]; //the position of black pixels
int pix_counter; //the number of pix on each line/column/height...


int config[MAXX][MAXX];


int lineS2[MAXX][Nt];
int columeS2[MAXX][Nt];

int N2V[MAXX][Nt];
int N2H[MAXX][Nt];



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



void read_config()
{
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      config[i][j] = 0;


  FILE* fp;

  if((fp=fopen("Mconfig.txt","r"))==NULL)
    {
      printf("No Mconfig.txt is found! Abort!\n");
      exit(1);
    }

  int tempx;
  int tempy;

  fscanf(fp, "%d", &NP);
  int temp_NP = 0;

  for(int i=0; i<NP; i++)
    {

      fscanf(fp, "%d", &tempx);
      fscanf(fp, "%d", &tempy);

      if(tempx<MAXX && tempy<MAXX)
	{
	  config[tempx][tempy] = 1;

	  temp_NP++;
	}
    }


  NP = temp_NP;
  cout<<"NP = "<<NP<<endl;

  fclose(fp);
}






void sampleS2line(int index1)
{
  
  for(int r=0; r<Nt; r++)
    {
      lineS2[index1][r] = 0;
    }

  //serach the line for pixel positions
  pix_counter = 0;

  for(int i=0; i<MAXX; i++)
    {
      if(config[i][index1] == 1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	lineS2[index1][temp_dist]++;

	//if(temp_dist == 0) temp_np++;
      }

  //temp_npL += lineS2[index1][index2][0];
  

  //cout<<"s2Line_Np = "<<lineS2[index1][index2][0]<<endl;
}



void sampleS2colume(int index1)
{

  for(int r=0; r<Nt; r++)
    {
      columeS2[index1][r] = 0;
    }
  
  //serach the line for pixel positions
  pix_counter = 0;

  for(int i=0; i<MAXX; i++)
    {
      if(config[index1][i]==1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	columeS2[index1][temp_dist]++;

	//if(temp_dist == 0) temp_np++;
	
      }

  // temp_npC += columeS2[index1][index2][0];
}




void sample_horizontal(int lindex)
{
  for(int r=0; r<Nt; r++)
    N2H[lindex][r] = 0;

  int ener[MAXX];
  int flag_empty = 0;

  for(int i=0; i<MAXX; i++)
    {
      if(config[lindex][i]==0) ener[i]=-1;
      else
	{
	  int en = 0;
	  int neb1 = i - 1;
	  if(neb1<0) neb1 = neb1 + MAXX;
	  if(config[lindex][neb1]==1) en++;
	  int neb2 = i+1;
	  if(neb2>=MAXX) neb2 = neb2 - MAXX;
	  if(config[lindex][neb2]==1) en++;

	  ener[i] = en;
	  flag_empty ++;
	}
    }

  int position[MAXX];
  for(int i=0; i<MAXX; i++)
    {
      position[i] = -1;
    }

  int ctp = 0;
  for(int i=0; i<MAXX; i++)
    {
      if(ener[i]==1)
	{
	  position[ctp] = i;
	  ctp++;
	}
      else if(ener[i]==0) //this counts for the individula pixels...
	{
	  N2H[lindex][0]++;
	}
    }

  if(config[lindex][0]==1&&config[lindex][MAXX-1]==1)//when the cord is at the bd
    {
      if(ctp>2)
	{
	  for(int i=1; i<ctp-1; i=i+2)
	    {
	      int len = position[i+1]-position[i]+1;
	      for(int r=0; r<=len; r++)
		{
		  if(r<Nt) N2H[lindex][r] = N2H[lindex][r]+(len-r);
		}
	    }

	  int len = (position[0]+1)+(MAXX - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<Nt) N2H[lindex][r] = N2H[lindex][r]+(len-r);
	    }
	}
      else if(ctp ==2)
	{
	  int len = (position[0]+1)+(MAXX - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<Nt) N2H[lindex][r] = N2H[lindex][r]+(len-r);
	    }
	}
      else if(ctp == 0 & flag_empty != 0)
	{
	  int len = MAXX;
	  for(int r=0; r<=len; r++)
	    {
	      if(r<Nt) N2H[lindex][r] = N2H[lindex][r]+(len-r);
	    }
	}
      
    }
  else
    {
      for(int i=0; i<ctp; i=i+2)
	{
	  int len = position[i+1]-position[i]+1;
	  for(int r= 0; r<=len; r++)
	    {
	      if(r<Nt) N2H[lindex][r] = N2H[lindex][r]+(len-r);
	    }
	}
    }
  


}

void sample_vertical(int cindex)
{
 for(int r=0; r<Nt; r++)
    N2V[cindex][r] = 0;

  int ener[MAXX];
  int flag_empty = 0;
  for(int i=0; i<MAXX; i++)
    {
      if(config[i][cindex]==0) ener[i]=-1;
      else
	{
	  int en = 0;
	  int neb1 = i - 1;
	  if(neb1<0) neb1 = neb1 + MAXX;
	  if(config[neb1][cindex]==1) en++;
	  int neb2 = i+1;
	  if(neb2>=MAXX) neb2 = neb2 - MAXX;
	  if(config[neb2][cindex]==1) en++;

	  ener[i] = en;
	  flag_empty++;
	}
    }

  int position[MAXX];
  for(int i=0; i<MAXX; i++)
    {
      position[i] = -1;
    }

  int ctp = 0;
  for(int i=0; i<MAXX; i++)
    {
      if(ener[i]==1)
	{
	  position[ctp] = i;
	  ctp++;
	}
      else if(ener[i]==0) //this counts for the individula pixels...
	{
	  N2V[cindex][0]++;
	}
    }

  if(config[0][cindex]==1&&config[MAXX-1][cindex]==1)//when the cord is at the bd
    {
      if(ctp>2)
	{
	  for(int i=1; i<ctp-1; i=i+2)
	    {
	      int len = position[i+1]-position[i]+1;
	      for(int r=0; r<=len; r++)
		{
		  if(r<Nt) N2V[cindex][r] = N2V[cindex][r]+(len-r);
		}
	    }

	  int len = (position[0]+1) + (MAXX - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<Nt) N2V[cindex][r] = N2V[cindex][r]+(len-r);
	    }
	}
      else if(ctp ==2)
	{
	  int len = (position[0]+1)+(MAXX - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<Nt) N2V[cindex][r] = N2V[cindex][r]+(len-r);
	    }
	}
      else if(ctp == 0 & flag_empty != 0)
	{
	  int len = MAXX;
	  for(int r=0; r<=len; r++)
	    {
	      if(r<Nt) N2V[cindex][r] = N2V[cindex][r]+(len-r);
	    }
	}
      
    }
  else
    {
      for(int i=0; i<ctp; i=i+2)
	{
	  int len = position[i+1]-position[i]+1;
	  for(int r= 0; r<=len; r++)
	    {
	      if(r<Nt) N2V[cindex][r] = N2V[cindex][r]+(len-r);
	    }
	}
    }
  

 
}




main()
{
  double S2X[Nt];
  double S2Y[Nt];
  double ST2[Nt];
  
  double LX[Nt];
  double LY[Nt];
  double LT[Nt];

  int SS2X[Nt];
  int SS2Y[Nt];
  int SLX[Nt];
  int SLY[Nt];


  for(int i=0; i<Nt; i++)
    {
      S2X[i] = 0;
      S2Y[i] = 0;
      ST2[i] = 0;
      SS2X[i] = 0;
      SS2Y[i] = 0;

      LX[i] = 0;
      LY[i] = 0;
      LT[i] = 0;
      SLX[i] = 0;
      SLY[i] = 0;
    }

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<Nt; j++)
      {
    	lineS2[i][j] = 0;
    	columeS2[i][j] = 0;

    	N2H[i][j] = 0;
    	N2V[i][j] = 0;
      }




  
  
  read_config();
 
  //init_config();//initialize configuration. volume fraction preserved...

  // read_config();

  //**************************************************
  //now we sample S2 for the first time...

  cout<<"sampling S_2 now..."<<endl;

  for(int r=0; r<Nt; r++)
    {
      for(int i=0; i<MAXX; i++)
    	{
	      sampleS2line(i);
	      sampleS2colume(i);

	       SS2X[r] = SS2X[r] + lineS2[i][r];
	       SS2Y[r] = SS2Y[r] + columeS2[i][r];
         }

      S2X[r] = (double) SS2X[r]/(double)(MAXX*MAXX);
      S2Y[r] = (double) SS2Y[r]/(double)(MAXX*MAXX);

      printf("%d\t%f\t%f\n", r, S2X[r], S2Y[r]);

    }

  printf("**********************************************\n");

  FILE* file = fopen("TS2.txt", "w");
  fprintf(file, "************** S2X ************* \n");
  for(int r=0; r<Nt; r++)
    {
          fprintf(file, "%d\t%f\n", r, S2X[r]);
    }
  fprintf(file, "************** S2Y ************* \n");
  for(int r=0; r<Nt; r++)
    {
          fprintf(file, "%d\t%f\n", r, S2Y[r]);
    }
  fclose(file);


  //now we sample L for the first time...
  cout<<"sampling L now..."<<endl;


  for(int r=0; r<Nt; r++)
    {
      for(int i=0; i<MAXX; i++)
    	{
     	  sample_horizontal(i);
     	  sample_vertical(i);

    	  SLX[r] = SLX[r] + N2H[i][r] ;
    	  SLY[r] = SLY[r] + N2V[i][r];
    	}

      LX[r] = (double) SLX[r]/(double)(MAXX*MAXX);
      LY[r] = (double) SLY[r]/(double)(MAXX*MAXX);

      printf("%d\t%f\t%f\n", r, LX[r], LY[r]);

    }

   file = fopen("TL.txt", "w");
   fprintf(file, "************** LX ************* \n");
   for(int r=0; r<Nt; r++)
    {
          fprintf(file, "%d\t%f\n", r, LX[r]);
    }
   fprintf(file, "************** LY ************* \n");
   for(int r=0; r<Nt; r++)
    {
          fprintf(file, "%d\t%f\n", r, LY[r]);
    }
  fclose(file);

  printf("******************************************\n");
      

  //this is the end of the codes...

}
