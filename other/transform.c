//#####################################################################################################################
//#This c/c++ code serves two purposes:                                                                               #
//#1. Add gaps to sinograms from VersaPET scanner to match STIR requirement in both transverse and axial planes       #
//#2. Sort the sinogram data output from Michelogram to STIR projection data format                                   #
//#Functions of the code:                                                                                             #
//#1. The code takes the sinogram .scn from the scanner and convert it to Michelogram                                #
//#2. The code converts the Michelogram to STIR Projection data without adding gaps                                   #
//#3. The code adds gaps to Michelogram and then coverts the gap-added Michelogram to STIR projection data            #
//#Author: Alex Shouyi Wei, PhD student in Dept. Biomedical Engineering, Stony Brook University                       #
//#Feb. 14, 2016, Stony Brook, New York, USA                                                                          #
//#####################################################################################################################

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

using namespace std;

#define N_RINGS 16
#define MAX_D_RING 15
#define N_SLICE 256
#define N_DET 192
#define S_WIDTH 191
#define N_DET_GAP 216
#define S_WIDTH_GAP 215
#define N_RINGS_NEW 73
#define MAX_D_RING_NEW 72

//currently I don't know how exactly the 3D sinograms from VersaPET scanner are ordered, I defaulted the order of segments as 0, negative 1, positive 1, negative 2, positive 2,..., but it could be the opposite as 0, positive 1, negative 1,...
//following is a lookup table that corresponds the indexes of 3D sinograms from VersaPET scanner to those of Michelograms which is the default one
//int raw_sinogram_index[N_SLICE] = {1, 18, 35, 52, 69, 86, 103, 120, 137, 154, 171, 188, 205, 222, 239, 256, 2, 19, 36, 53, 70, 87, 104, 121, 138, 155, 172, 189, 206, 223, 240, 17, 34, 51, 68, 85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255, 3, 20, 37, 54, 71, 88, 105, 122, 139, 156, 173, 190, 207, 224, 33, 50, 67, 84, 101, 118, 135, 152, 169, 186, 203, 220, 237, 254, 4, 21, 38, 55, 72, 89, 106, 123, 140, 157, 174, 191, 208, 49, 66, 83, 100, 117, 134, 151, 168, 185, 202, 219, 236, 253, 5, 22, 39, 56, 73, 90, 107, 124, 141, 158, 175, 192, 65, 82, 99, 116, 133, 150, 167, 184, 201, 218, 235, 252, 6, 23, 40, 57, 74, 91, 108, 125, 142, 159, 176, 81, 98, 115, 132, 149, 166, 183, 200, 217, 234, 251, 7, 24, 41, 58, 75, 92, 109, 126, 143, 160, 97, 114, 131, 148, 165, 182, 199, 216, 233, 250, 8, 25, 42, 59, 76, 93, 110, 127, 144, 113, 130, 147, 164, 181, 198, 215, 232, 249, 9, 26, 43, 60, 77, 94, 111, 128, 129, 146, 163, 180, 197, 214, 231, 248, 10, 27, 44, 61, 78, 95, 112, 145, 162, 179, 196, 213, 230, 247, 11, 28, 45, 62, 79, 96, 161, 178, 195, 212, 229, 246, 12, 29, 46, 63, 80, 177, 194, 211, 228, 245, 13, 30, 47, 64, 193, 210, 227, 244, 14, 31, 48, 192, 209, 226, 15, 32, 208, 225, 16, 207}; 

//another possibility 
int raw_sinogram_index[N_SLICE] = {1, 18, 35, 52, 69, 86, 103, 120, 137, 154, 171, 188, 205, 222, 239, 256, 17, 34, 51, 68, 85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255, 2, 19, 36, 53, 70, 87, 104, 121, 138, 155, 172, 189, 206, 223, 240, 33, 50, 67, 84, 101, 118, 135, 152, 169, 186, 203, 220, 237, 254, 3, 20, 37, 54, 71, 88, 105, 122, 139, 156, 173, 190, 207, 224, 49, 66, 83, 100, 117, 134, 151, 168, 185, 202, 219, 236, 253, 4, 21, 38, 55, 72, 89, 106, 123, 140, 157, 174, 191, 208, 65, 82, 99, 116, 133, 150, 167, 184, 201, 218, 235, 252, 5, 22, 39, 56, 73, 90, 107, 124, 141, 158, 175, 192, 81, 98, 115, 132, 149, 166, 183, 200, 217, 234, 251, 6, 23, 40, 57, 74, 91, 108, 125, 142, 159, 176, 97, 114, 131, 148, 165, 182, 199, 216, 233, 250, 7, 24, 41, 58, 75, 92, 109, 126, 143, 160, 113, 130, 147, 164, 181, 198, 215, 232, 249,  8, 25, 42, 59, 76, 93, 110, 127, 144, 129, 146, 163, 180, 197, 214, 231, 248, 9, 26, 43, 60, 77, 94, 111, 128, 145, 162, 179, 196, 213, 230, 247, 10, 27, 44, 61, 78, 95, 112, 161, 178, 195, 212, 229, 246, 11, 28, 45, 62, 79, 96, 177, 194, 211, 228, 245, 12, 29, 46, 63, 80, 193, 210, 227, 244, 13, 30, 47, 64, 192, 209, 226, 14, 31, 48, 208, 225, 15, 32, 207, 16}; 

unsigned short ans;
float Mich_sfu[N_RINGS][N_RINGS][N_DET/2][S_WIDTH]={0};
float Mich_gapsfu[N_RINGS][N_RINGS][N_DET_GAP/2][S_WIDTH_GAP]={0};
//float Mich_tagap_sfu[N_RINGS_NEW][N_RINGS_NEW][N_DET_GAP/2][S_WIDTH_GAP]={0};

int main(int argc, char** argv){
    string filedir, inputfilename;
    string filename, Moutputfilename, Poutputfilename, Poriginaloutputfilename, Moriginaloutputfilename;
    int ring1, ring2;
    int id1=0;
    int id2=0;
    int phi_new=0;
    int u_new=0;
    int index1=0;
    int index2=0;
    int ind_r1 = 0;
    int ind_r2 = 0;
    FILE * file=fopen("crystalids_transform.txt","w");
    if(argc<2) {
    	cout<<" Right number of input argument please !! "<<endl ;
    	return 1;
    }
    
    inputfilename = argv[1];
    cout<<"input file name is"<<inputfilename<<endl;
    filename = argv[2];
    Moutputfilename = "Mich_" + filename + ".s";
    cout << "New Michelogram file name is = " << Moutputfilename << endl; 
    Poutputfilename = "Proj_"+ filename + ".s" ;
    cout << "New Projection file name is = " << Poutputfilename << endl ; //GE/STIR format sinogram
    Moriginaloutputfilename = "Michorigin_" + filename + ".s";     
    Poriginaloutputfilename = "Projorigin_" + filename + ".s";
    FILE *Mich_File;
    FILE *Proj_File;
   // FILE *newProj_File;
    FILE *newMich_File;
    FILE *newProj_File;
    Mich_File = fopen(Moutputfilename.c_str(), "wb");
    Proj_File = fopen(Poutputfilename.c_str(),"wb");
    newMich_File = fopen(Moriginaloutputfilename.c_str(), "wb");   
    newProj_File = fopen(Poriginaloutputfilename.c_str(), "wb");
    FILE *rawdata = fopen(inputfilename.c_str(), "rb");
    vector<float> buffer(N_SLICE*N_DET/2*S_WIDTH);
    fread(&buffer[0], sizeof(float), buffer.size(), rawdata);
    //here convert the 3D sinogram to Michelogram, and I can also add axial gaps to Michelogram
    for(int i=0; i<N_SLICE; i++){
      for(int j=0; j<N_DET/2; j++){
        for(int k=0; k<S_WIDTH; k++){
            index1 = (raw_sinogram_index[i]-1)/16;
            index2 = (raw_sinogram_index[i]-1)%16;  
            Mich_sfu[index1][index2][j][k] = (float)buffer[i*N_DET/2*S_WIDTH+j*S_WIDTH+k];
            
         }
       }
     }
 
    for(int h=0; h<N_RINGS; h++){
      for(int i=0; i<N_RINGS; i++){
	for(int j=0; j<N_DET/2; j++){
          for(int k=0; k<S_WIDTH; k++){
              if(Mich_sfu[h][i][j][k] != 0){
                for(int m=0; m<N_DET; m++){
                  for(int n=0; n<N_DET; n++){
                     if(((m+n+N_DET/2)%N_DET)/2 == j){
                          if (((m+n) < (3*N_DET/2)) && ((m+n) >= (N_DET/2))){
                              if(k==abs(m-n)-N_DET/2+S_WIDTH/2){// find a crystal pair
                              //change the m and n ID
                                     
		                         id1 = m+(int)((m+4)/8);
		                         id2 = n+(int)((n+4)/8); 
                                  
                                      phi_new = ((id1 + id2 + N_DET_GAP/2)%N_DET_GAP)/2;     
		                      if (((id1+id2) < (3*N_DET_GAP/2)) && ((id1+id2) >= (N_DET_GAP/2)))
					   u_new = abs(id1 - id2) -  N_DET_GAP/2 + S_WIDTH_GAP/2;
	     			      else u_new = -abs(id1 - id2) +  N_DET_GAP/2 + S_WIDTH_GAP/2;                                
                                      cout<<"h="<<h<<" "<<"i="<<i<<" "<<"j="<<j<<" "<<"k="<<k<<" "<<"m="<<m<<" "<<"n="<<n<<"id1="<<id1<<" "<<"id2="<<id2<<"phi="<<phi_new<<" "<<"u="<<u_new<<" "<<"counts="<<Mich_sfu[h][i][j][k]<<endl; 
                                      //fprintf (file, "phi = %d, u = %d, crystal1 = %d, crystal2 = %d\n", j, k, m, n);
                                      Mich_gapsfu[h][i][phi_new][u_new]=Mich_sfu[h][i][j][k];  //new Michelogram  
                              } 
                          }
                         else{
                                if(k==-abs(m-n)+N_DET/2+S_WIDTH/2){// find a crystal pair
                                //change the m and n ID
		                
		                         id1 = m+(int)((m+4)/8);
                           		 id2 = n+(int)((n+4)/8); 
                                       
                                      phi_new = ((id1 + id2 + N_DET_GAP/2)%N_DET_GAP)/2;
	      			     // if(phi_new>N_DET_GAP/4) phi_new = phi_new-N_DET_GAP/4; //offset
	     			     // else phi_new = phi_new+N_DET_GAP/4;   //offset           
		                      if (((id1+id2) < (3*N_DET_GAP/2)) && ((id1+id2) >= (N_DET_GAP/2)))
					   u_new = abs(id1 - id2) - N_DET_GAP/2 + S_WIDTH_GAP/2;
	     			      else u_new = -abs(id1 - id2) +  N_DET_GAP/2 + S_WIDTH_GAP/2; 
                                     
                                      cout<<"h="<<h<<" "<<"i="<<i<<" "<<"j="<<j<<" "<<"k="<<k<<" "<<"m="<<m<<" "<<"n="<<n<<"id1="<<id1<<" "<<"id2="<<id2<<"phi="<<phi_new<<" "<<"u="<<u_new<<" "<<"counts="<<Mich_sfu[h][i][j][k]<<endl; 
                                      //fprintf (file, "phi = %d, u = %d, crystal1 = %d, crystal2 = %d\n", j, k, m, n);
                                      Mich_gapsfu[h][i][phi_new][u_new]=Mich_sfu[h][i][j][k];  //new Michelogram
                                }  
                          } 
                          
 
                    }
              
                }
              }
             }
           }
         }

      }  
    }

   /* 
    for(int h=0; h<N_RINGS; h++){
      for(int i=0; i<N_RINGS; i++){
        ind_r1 = h*4+int(h/4)*3;
        ind_r2 = i*4+int(i/4)*3;   
	for(int j=0; j<N_DET_GAP/2; j++){
          for(int k=0; k<S_WIDTH_GAP; k++){
             Mich_tagap_sfu[ind_r1][ind_r2][j][k] = Mich_gapsfu[h][i][j][k];
          }
        }
       }
    }
*/
    ans = fwrite(Mich_gapsfu,4,((N_RINGS)*(N_RINGS)*N_DET_GAP/2*S_WIDTH_GAP),Mich_File);
    fclose(Mich_File);
   // ans = fwrite(Mich_tagap_sfu,4,((N_RINGS_NEW)*(N_RINGS_NEW)*N_DET_GAP/2*S_WIDTH_GAP),newMich_File);
    //ans = fwrite(Mich_gapsfu,4,((N_RINGS)*(N_RINGS)*N_DET_GAP/2*S_WIDTH_GAP),newMich_File);
  
    ans = fwrite(Mich_sfu,4,((N_RINGS )*(N_RINGS)*N_DET/2*S_WIDTH),newMich_File);
    fclose(newMich_File);
    int S_NUM; //number of segments

    for (int ii = 0 ; ii < 2*MAX_D_RING + 1 ; ii++)
    {
      if (ii <= MAX_D_RING) S_NUM = N_RINGS - MAX_D_RING + ii;
      else                 S_NUM = N_RINGS + MAX_D_RING - ii;
  
      for (int jj = 0 ; jj < N_DET_GAP/2 ; jj++)
	{ 
	  float Proj[S_NUM][S_WIDTH_GAP];
	  for (int kk = 0 ; kk < S_NUM ; kk++)
	    {
	      if (ii <= MAX_D_RING) ring1 = kk;
	      else                 ring2 = kk;

	      if (ii <= MAX_D_RING) ring2 = ring1 + MAX_D_RING - ii;
	      else                 ring1 = ring2 - MAX_D_RING + ii;

	      for (int ll = 0 ; ll < S_WIDTH_GAP ; ll++)
		Proj[kk][ll] = Mich_gapsfu[ring2][ring1][jj][ll];
	    }
	  ans = fwrite(Proj,4,(S_NUM*S_WIDTH_GAP),Proj_File); 
	}
    }
  
    fclose(Proj_File); //Projection data without gaps

/*
    for (int ii = 0 ; ii < 2*MAX_D_RING_NEW + 1 ; ii++)
    {
      if (ii <= MAX_D_RING_NEW) S_NUM = N_RINGS_NEW - MAX_D_RING_NEW + ii;
      else                 S_NUM = N_RINGS_NEW + MAX_D_RING_NEW - ii;
  
      for (int jj = 0 ; jj < N_DET_GAP/2 ; jj++)
	{ 
	  float newProj[S_NUM][S_WIDTH_GAP];
	  for (int kk = 0 ; kk < S_NUM ; kk++)
	    {
	      if (ii <= MAX_D_RING_NEW) ring1 = kk;
	      else                 ring2 = kk;

	      if (ii <= MAX_D_RING_NEW) ring2 = ring1 + MAX_D_RING_NEW - ii;
	      else                 ring1 = ring2 - MAX_D_RING_NEW + ii;

	      for (int ll = 0 ; ll < S_WIDTH_GAP ; ll++)
		newProj[kk][ll] = Mich_tagap_sfu[ring2][ring1][jj][ll];
	    }
	  ans = fwrite(newProj,4,(S_NUM*S_WIDTH_GAP),newProj_File); 
	}
    }
  
    fclose(newProj_File); //Projection data with gaps
*/ 
  for (int ii = 0 ; ii < 2*MAX_D_RING + 1 ; ii++)
    {
      if (ii <= MAX_D_RING) S_NUM = N_RINGS - MAX_D_RING + ii;
      else                 S_NUM = N_RINGS + MAX_D_RING - ii;
  
      for (int jj = 0 ; jj < N_DET/2 ; jj++)
	{ 
	  float Proj1[S_NUM][S_WIDTH];
	  for (int kk = 0 ; kk < S_NUM ; kk++)
	    {
	      if (ii <= MAX_D_RING) ring1 = kk;
	      else                 ring2 = kk;

	      if (ii <= MAX_D_RING) ring2 = ring1 + MAX_D_RING - ii;
	      else                 ring1 = ring2 - MAX_D_RING + ii;

	      for (int ll = 0 ; ll < S_WIDTH ; ll++)
		Proj1[kk][ll] = Mich_sfu[ring2][ring1][jj][ll];
	    }
	  ans = fwrite(Proj1,4,(S_NUM*S_WIDTH),newProj_File); 
	}
    }
  
    fclose(newProj_File); //Projection data without gaps


    return 0;

}
    
     








    
        
        
