//  *********************************************************************
//  * DISCLAIMER                                                        *
//  *                                                                   *
//  * Neither the authors of this software system, nor their employing  *
//  * institutes, nor the agencies providing financial support for this *
//  * work  make  any representation or  warranty, express or implied,  *
//  * regarding  this  software system or assume any liability for its  *
//  * use.                                                              *
//  *                                                                   *
//  * This  code  implementation is the  intellectual property  of the  *
//  * OpenGATE collaboration.                                           *
//  * By copying,  distributing  or modifying the Program (or any work  *
//  * based  on  the Program)  you indicate  your  acceptance of  this  *
//  * statement, and all its terms.                                     *
//  *********************************************************************
//
//######################################################################################
//# Program   : VersaPET_processing.c                                                  #
//# Modified by Alex Shouyi Wei for VersaPET geometry in Stony Brook University        #
//# Jan, 2016, NY, USA                                                                 #
//# Original code authors                                                              #
//# Authors   : Sadek A. Nehmeh, CR Schmidtlein                                        # 
//#                                                                                    #
//# Program   : Bin_GATE_v1.0.c  29-JUL-2010                                           # 
//#                                                                                    #   
//# Objective : To read the coincidences TTree from the .root file, and generates the  #
//#             corresponding Michelogram and Projection files.                        #
//#                                                                                    #
//# Input     : Monte Carlo data from GATE and egsPET                                  #
//#                                                                                    #
//# Output    : 1 Michelogram files according to various binning definitions           #
//#           : 2 Projection  files according to various binning definitions           #
//#                                                                                    #
//# HOW TO COMPILE:                                                                    #
//# 1) Compile using this command line in the terminal:                                #
//#       g++ (name of file) `root-config --cflags --libs`                             #
//#                                                                                    #
//# HOW TO RUN:                                                                        #
//# 1) After compiling the code type the following command line:                       #
//#       ./a.out 'directory name of root files' 'output file name'                    #
//#                                                                                    #
//######################################################################################
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TH2D.h" 
#include "TDirectory.h"
#include "TList.h"
#include "Rtypes.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TRandom.h"

using namespace std ;
/*
// DLS Parameters
#define  N_RINGS               18 // DLS: Number of rings
#define  N_SEG                 35 // DLS: Number of segments (2 N_RINGS - 1)
#define  N_DET                672 // DLS: Detectors per ring
#define  S_WIDTH              283 // DLS: Number of radial sinogram bins
#define  N_RSEC                56 // DLS: Number of rsector
#define  N_MODULE               6 // DLS: Number of modules (2x3)
#define  N_MOD_xy               2 // DLS: Number of tangential modules
#define  N_MOD_z                3 // DLS: Number of axial modules
#define  N_SUBMOD               1 // DLS: Number of submodules (1x1)
#define  N_SMOD_xy              1 // DLS: Number of tangential submodules
#define  N_SMOD_z               1 // DLS: Number of axial submodules
#define  N_CRYSTAL             36 // DLS: Number of crystals (6x6)
#define  N_CRY_xy               6 // DLS: Number of tangential crystals
#define  N_CRY_z                6 // DLS: Number of axial crystals
#define  MAX_D_RING            11 // DLS: Maximum ring difference
#define  N_PLANES             265 // DLS: Total number of sinograms
#define  FOV                  550 // DLS: Width of the FOV (mm)
#define  SLICES_PER_FOV        35 // DLS: Number of slices per FOV
#define  USE_OFFSET             1 // DLS: On/Off use of offset
#define  OFFSET               499 // DLS: Sets initial sinogram angle
*/
// VersaPET Parameters
#define N_RINGS                16 // VersaPET: Number of rings
#define N_SEG                  31 // VersaPET: Number of segments (2 * N_RINGS - 1)
#define N_DET                 216 // VersaPET: Detectors per ring
#define S_WIDTH               215 // VersaPET: Number of radial sinogram bins
#define N_RSEC                 24 // VersaPET: Number of rsector
#define N_MODULE                4 // VersaPET: Number of modules 1x4
#define N_MOD_xy                1 // VersaPET: Number of tangential modules 
#define N_MOD_z                 4 // VersaPET: Number of axial modules
#define N_CRYSTAL              32 // VersaPET: Number of crystals 8x4
#define N_CRY_xy                8 // VersaPET: Number of tangential crystals 
#define N_CRY_z                 4 // VersaPET: Number of axial crystals
#define MAX_D_RING             15 // VersaPET: Maximum ring difference
#define USE_OFFSET		0 // VersaPET: On/off use of offset
#define OFFSET                 99 // VersaPET: Sets initial sinogram angle

unsigned short ans,ans1;
unsigned short ans_SC, ans_RC, ans_SRC, ans_S, ans_R;
/*
  Current data constructions available for output.
  These are: Michelograms:  [ring1][ring2][phi][u]
  
  total counts, trues, scatter, and randoms can be independently collected
*/
Float_t    Mich_r1r2fu[N_RINGS][N_RINGS][N_DET/2][S_WIDTH]={0}; 

int main(int argc, char** argv)
{
  //---------------------------------------------------------------------------------
  // the first  argument (argv[1]) is the sub-directory of the input file
  // the second argument (argv[2]) is the name of the output file
  //
  //---------------------------------------------------------------------------------
  if(argc<2) {
    std::cout<<" Right number of input argument please !! "<<std::endl ;
    return 1;
  }
  
  Double_t PI = acos(-1.0);

  string filedir, inputfilename ;
  string filename, Moutputfilename, Poutputfilename ;
 
  filedir = argv[1] ; 
  inputfilename = filedir + "*.root" ;

  cout << "Input file name is " << inputfilename << endl;
  TChain *Coincidences = new TChain("Coincidences") ;
  Coincidences->Add(inputfilename.c_str()) ; 
  
  filename = argv[2] ; 
  Moutputfilename = "Mich_"+ filename + ".s" ;
  cout << "Michelogram file name is = " << Moutputfilename << endl ;
  Poutputfilename = "Proj_"+ filename + ".s" ;
  cout << "Projection file name is = " << Poutputfilename << endl ;

  FILE  *Mich_r1r2fuFile, *Proj_File;      
  Mich_r1r2fuFile = fopen(Moutputfilename.c_str(),"wb");
  Proj_File       = fopen(Poutputfilename.c_str(),"wb");

//#####################################################################
//#              Loop over the .root file in the directory "PATH"     #
//#####################################################################
  Int_t   Trues = 0, Scatters = 0, Randoms = 0;
  Int_t   nbytes = 0;
  
  //####################################################################
  //#             Declaration of leaves types - TTree Coincidences     #
  //####################################################################
  
  Float_t         axialPos, rotationAngle, sinogramS, sinogramTheta;
  Char_t          comptVolName1[40], comptVolName2[40];
  Int_t           compton1, compton2;
  Int_t           runID, sourceID1, sourceID2, eventID1, eventID2; 
  Int_t           layerID1, layerID2, crystalID1, crystalID2;
  Int_t           submoduleID1, submoduleID2, moduleID1, moduleID2, rsectorID1, rsectorID2;
  Int_t           comptonPhantom1, comptonPhantom2;
  Float_t         energy1, energy2;   
  Float_t         globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
  Float_t         sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2;
  Double_t        time1, time2;
  
  //######################################################################################
  //#                        Set branch addresses - TTree Coincidences                   #
  //######################################################################################
  
  Coincidences->SetBranchStatus("*",0);
  Coincidences->SetBranchAddress("axialPos",&axialPos);
  Coincidences->SetBranchAddress("comptVolName1",&comptVolName1);
  Coincidences->SetBranchAddress("comptVolName2",&comptVolName2);
  Coincidences->SetBranchAddress("comptonCrystal1",&compton1);
  Coincidences->SetBranchAddress("comptonCrystal2",&compton2);
  Coincidences->SetBranchAddress("crystalID1",&crystalID1);
  Coincidences->SetBranchAddress("crystalID2",&crystalID2);
  Coincidences->SetBranchAddress("comptonPhantom1",&comptonPhantom1);
  Coincidences->SetBranchAddress("comptonPhantom2",&comptonPhantom2);
  Coincidences->SetBranchAddress("energy1",&energy1);
  Coincidences->SetBranchAddress("energy2",&energy2);   
  Coincidences->SetBranchAddress("eventID1",&eventID1);
  Coincidences->SetBranchAddress("eventID2",&eventID2);
  Coincidences->SetBranchAddress("globalPosX1",&globalPosX1);
  Coincidences->SetBranchAddress("globalPosX2",&globalPosX2);
  Coincidences->SetBranchAddress("globalPosY1",&globalPosY1);
  Coincidences->SetBranchAddress("globalPosY2",&globalPosY2);
  Coincidences->SetBranchAddress("globalPosZ1",&globalPosZ1);
  Coincidences->SetBranchAddress("globalPosZ2",&globalPosZ2);
  Coincidences->SetBranchAddress("layerID1",&layerID1);
  Coincidences->SetBranchAddress("layerID2",&layerID2);
  Coincidences->SetBranchAddress("moduleID1",&moduleID1);
  Coincidences->SetBranchAddress("moduleID2",&moduleID2);
  Coincidences->SetBranchAddress("rotationAngle",&rotationAngle);
  Coincidences->SetBranchAddress("rsectorID1",&rsectorID1);
  Coincidences->SetBranchAddress("rsectorID2",&rsectorID2);
  Coincidences->SetBranchAddress("runID",&runID);
  Coincidences->SetBranchAddress("sinogramS",&sinogramS);
  Coincidences->SetBranchAddress("sinogramTheta",&sinogramTheta);
  Coincidences->SetBranchAddress("sourceID1",&sourceID1);
  Coincidences->SetBranchAddress("sourceID2",&sourceID2);
  Coincidences->SetBranchAddress("sourcePosX1",&sourcePosX1);
  Coincidences->SetBranchAddress("sourcePosX2",&sourcePosX2);
  Coincidences->SetBranchAddress("sourcePosY1",&sourcePosY1);
  Coincidences->SetBranchAddress("sourcePosY2",&sourcePosY2);
  Coincidences->SetBranchAddress("sourcePosZ1",&sourcePosZ1);
  Coincidences->SetBranchAddress("sourcePosZ2",&sourcePosZ2);
  Coincidences->SetBranchAddress("submoduleID1",&submoduleID1);
  Coincidences->SetBranchAddress("submoduleID2",&submoduleID2);
  Coincidences->SetBranchAddress("time1",&time1);
  Coincidences->SetBranchAddress("time2",&time2);
    
  Int_t nentries = (Int_t)(Coincidences->GetEntries());
  
  //#####################################################################
  //#             SINOGRAMS AND PROJECTION PLANES BINNING               #
  //#####################################################################
  Int_t    ring1, ring2, crystal1, crystal2;
  Int_t    phi, u;
  Int_t    Counts = 0;  
  int flip, swap, zi, c1, c2;

  printf("Total Number of Coincidence Events:= %d \n",nentries ); 
  
  for (Int_t i = 0 ; i < nentries ; i++)
    {      

      if ((i%250000)  == 0 && i!=0)  printf("... %d ",i);       
      if ((i%1000000) == 0 && i!=0)  printf("... %d\n",i);       

      nbytes += Coincidences->GetEntry(i);

      // Update the number of Trues and Randoms...
      //------------------------------------------
      
      if (eventID1 == eventID2)  
	{
	  if (comptonPhantom1 == 0 && comptonPhantom2 == 0) Trues++;
	  else Scatters++;
	}
      else Randoms++;       

      //-----------------------------------
      //  Identify the ring#...
      //-----------------------
     /* ring1 = (Int_t)(crystalID1/N_CRY_xy) 
	    + (Int_t)(submoduleID1/N_SMOD_xy)*N_CRY_z
	    + (Int_t)(moduleID1/N_MOD_xy)*N_SMOD_z*N_CRY_z;
      ring2 = (Int_t)(crystalID2/N_CRY_xy) 
	    + (Int_t)(submoduleID2/N_SMOD_xy)*N_CRY_z
	    + (Int_t)(moduleID2/N_MOD_xy)*N_SMOD_z*N_CRY_z;
      */
     // cout<<"moduleID1="<<moduleID1<<endl;
     // cout<<"moduleID2="<<moduleID2<<endl;
     // cout<<"crystalID1="<<crystalID1<<endl;
     // cout<<"crystalID2="<<crystalID2<<endl;
      ring1 = (Int_t)(crystalID1/N_CRY_xy)
            +(Int_t)(moduleID1/N_MOD_xy)*N_CRY_z;
      ring2 = (Int_t)(crystalID2/N_CRY_xy)
            +(Int_t)(moduleID2/N_MOD_xy)*N_CRY_z;
     // cout<<"ring1="<<ring1<<endl;
     // cout<<"ring2="<<ring2<<endl;
      if ( abs(ring1 - ring2) > MAX_D_RING )  {
         cout<<"ring1="<<ring1<<endl;
         cout<<"ring2="<<ring2<<endl;
         continue;  
      }
      //-----------------------
      //  Identify the crystal#...
      //-----------------------------------
     /* crystal1 = rsectorID1 * N_MOD_xy * N_SMOD_xy * N_CRY_xy
	       + (moduleID1%N_MOD_xy) * N_SMOD_xy * N_CRY_xy
	       + (submoduleID1%N_SMOD_xy) * N_CRY_xy
	       + (crystalID1%N_CRY_xy);
      crystal2 = rsectorID2 * N_MOD_xy * N_SMOD_xy * N_CRY_xy
	       + (moduleID2%N_MOD_xy) * N_SMOD_xy * N_CRY_xy
	       + (submoduleID2%N_SMOD_xy) * N_CRY_xy
	       + (crystalID2%N_CRY_xy);
      */
       crystal1 = rsectorID1 * N_MOD_xy * N_CRY_xy
                + (crystalID1%N_CRY_xy)+rsectorID1;
       crystal2 = rsectorID2 * N_MOD_xy * N_CRY_xy
                + (crystalID2%N_CRY_xy)+rsectorID2;        
     //  cout<<"crystal1="<<crystal1<<endl;
     //  cout<<"crystal2="<<crystal2<<endl;

      //-----------------------------------------------------
      //  Rotate the image correctly#...
      //--------------------------------
      if (USE_OFFSET == 1)
	{
	  crystal1 = crystal1 + OFFSET;
	  crystal2 = crystal2 + OFFSET;    
	  if (crystal1 >= N_DET)  crystal1 = crystal1 - N_DET;
	  if (crystal2 >= N_DET)  crystal2 = crystal2 - N_DET;
	}
      //--------------------------------
      //  Bin the crystal ring pairs into Michelograms
      //  u - radial sinogram component
      //  phi - azimuthal sinogram component
      //  ring pairs are sorted according to c1 < c2 else flip
      //  where c1 and c2 are crystals at phi(u = S_WIDTH/2)
      //--------------------------------
      phi = ((crystal1 + crystal2 + N_DET/2)%N_DET)/2;
             

      if (((crystal1 + crystal2) < (3*N_DET/2)) && ((crystal1 + crystal2) >= (N_DET/2)))
	u    =  abs(crystal1 - crystal2) -  N_DET/2 + S_WIDTH/2;
      else u = -abs(crystal1 - crystal2) +  N_DET/2 + S_WIDTH/2;         

      if ( u >= S_WIDTH || u < 0 ) continue;

      if (u%2 == 0) 
	{
	  zi = (N_DET/2 - (crystal1 - crystal2) - 1)/2;
	  if (zi >=  N_DET/4) zi = zi - N_DET/2 + 1;
	  if (zi <= -N_DET/4) zi = zi + N_DET/2 - 1;
	}
      else          
	{
	  zi = (N_DET/2 - (crystal1 - crystal2))/2;
	  if (zi >=  N_DET/4) zi = zi - N_DET/2;
	  if (zi <= -N_DET/4) zi = zi + N_DET/2;
	}

      c1 = crystal1 + zi;
      c2 = crystal2 - zi;
      if (c1 >= N_DET) c1 = c1 - N_DET;
      if (c1 < 0)      c1 = c1 + N_DET;
      if (c2 >= N_DET) c2 = c2 - N_DET;
      if (c2 < 0)      c2 = c2 + N_DET;
      
      if (c1 < c2) flip = 0;
      else         flip = 1;

      if (flip) 
	{
	  swap  = ring1;
	  ring1 = ring2;
	  ring2 = swap;
	}

     // Update the different arrays...
     //-------------------------------
     //***ALL EVENTS
     Mich_r1r2fu[ring2][ring1][phi][u] += 1.;

     if (eventID1 == eventID2)
       {
	 //***true+scatter
	 if (comptonPhantom1 == 0 && comptonPhantom1 == 0)
	   {
	     //***true
	   }
	 else //azimuthal sinogram component
	   {
	     //***scatter
	   }
       }
     else
       {
	 //***random
       }
     Counts++;
	                 
    }
  printf("\n");

  // Write the data to disk, and then close Michelogram file...
  //----------------------------------------------------------------

  ans = fwrite(Mich_r1r2fu,4,(N_RINGS*N_RINGS*N_DET/2*S_WIDTH),Mich_r1r2fuFile);
  fclose(Mich_r1r2fuFile);
  
  // Generate projection files...
  //-----------------------------
  // From Segment number -MAX_D_RING to +MAX_D_RING
  // After phi
  // After z
  // After r
  Int_t S_NUM;
  for (int i = 0 ; i < 2*MAX_D_RING + 1 ; i++)
    {
      if (i <= MAX_D_RING) S_NUM = N_RINGS - MAX_D_RING + i;
      else                 S_NUM = N_RINGS + MAX_D_RING - i;
  
      for (int j = 0 ; j < N_DET/2 ; j++)
	{ 
	  float Proj[S_NUM][S_WIDTH];
	  for (int k = 0 ; k < S_NUM ; k++)
	    {
	      if (i <= MAX_D_RING) ring1 = k;
	      else                 ring2 = k;

	      if (i <= MAX_D_RING) ring2 = ring1 + MAX_D_RING - i;
	      else                 ring1 = ring2 - MAX_D_RING + i;

	      for (int l = 0 ; l < S_WIDTH ; l++)
		Proj[k][l] = Mich_r1r2fu[ring2][ring1][j][l];
	    }
	  ans = fwrite(Proj,4,(S_NUM*S_WIDTH),Proj_File); 
	}
    }
  
  fclose(Proj_File);
  printf("Total Number of True Coincidence Events:= %d \n",Trues ); 
  //-----------------------------
  
  return(0);
}

