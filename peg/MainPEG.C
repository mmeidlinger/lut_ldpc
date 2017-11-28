/************************************************************************/
/*                                                                      */
/*        Free software: Progressive edge-growth (PEG) algorithm        */
/*        Created by Xiaoyu Hu                                          */
/*                   Evangelos Eletheriou                               */
/*                   Dieter Arnold                                      */
/*        IBM Research, Zurich Research Lab., Switzerland               */
/*                                                                      */
/*        The C++ sources files have been compiled using xlC compiler   */
/*        at IBM RS/6000 running AIX. For other compilers and platforms,*/
/*        minor changes might be needed.                                */
/*                                                                      */
/*        Bug reporting to: xhu@zurich.ibm.com                          */
/**********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "BigGirth.h"
#include "Random.h"
#include "CyclesOfGraph.h"

#define EPS  1e-6

using namespace std;

int main(int argc, char * argv[]){
  int i, j, m, N, M;
  int sglConcent=1;  // default to non-strictly concentrated parity-check distribution
  int targetGirth=100000; // default to greedy PEG version 
  char codeName[100], degFileName[100];
  int *degSeq, *deg;
  double *degFrac;
  BigGirth *bigGirth;
  CyclesOfGraph *cog;

  int numArgs=(argc-1)/2;
  if (argc<9) {
  USE:
    cout<<"*******************************************************************************************"<<endl;
    cout<<" Usage Reminder: MainPEG -numM M -numN N -codeName CodeName -degFileName DegFileName " <<endl;
    cout<<"         option:         -sglConcent SglConcent                                     " <<endl; 
    cout<<"                         sglConcent==0 ----- strictly concentrated parity-check      " <<endl;
    cout<<"                                       degree distribution (including regular graphs)" <<endl;
    cout<<"                         sglConcent==1 ----- Best-effort concentrated (DEFAULT)      " <<endl;
    cout<<"         option:         -tgtGirth TgtGirth                                          " <<endl; 
    cout<<"                  TgtGirth==4, 6 ...; if very large, then greedy PEG (DEFAULT)       " <<endl;
    cout<<"                  IF sglConcent==0, TgtGirth is recommended to be set relatively small" <<endl;
    cout<<"                                                                                       " <<endl;
    cout<<" Remarks: File CodeName stores the generated PEG Tanner graph. The first line contains"<<endl;
    cout<<"          the block length, N. The second line defines the number of parity-checks, M."<<endl;
    cout<<"          The third line defines the number of columns of the compressed parity-check "<<endl;
    cout<<"          matrix. The following M lines are then the compressed parity-check matrix.  "<<endl;
    cout<<"          Each of the M rows contains the indices (1 ... N) of 1's in the compressed  "<<endl;
    cout<<"          row of parity-check matrix. If not all column entries are used, the column  "<<endl;
    cout<<"          is filled up with 0's.                                                      "<<endl;
    cout<<"                                                                                      "<<endl;
    cout<<"          File DegFileName is the input file to specify the degree distribution (node "<<endl;
    cout<<"          perspective). The first line contains the number of various degrees. The second"<<endl;
    cout<<"          defines the row vector of degree sequence in the increasing order. The vector"<<endl;
    cout<<"          of fractions of the corresponding degree is defined in the last line.         "<<endl;
    cout<<"                                                                                       "<<endl;
    cout<<"          A log file called 'leftHandGirth.dat' will also be generated and stored in the"<<endl;
    cout<<"          current directory, which gives the girth of the left-hand subgraph of j, where"<<endl;
    cout<<"          1<=j<=N. The left-hand subgraph of j is defined as all the edges emanating from"<<endl;
    cout<<"          bit nodes {1 ... j} and their associated nodes.                                "<<endl; 
    cout<<"                                                                                         "<<endl;
    cout<<"          The last point is, when strictly concentrated parity-check degree distribution"<<endl;
    cout<<"          is invoked, i.e. sglConcent==0, the girth might be weaken to some extent as    "<<endl;
    cout<<"          compared to the generic PEG algorithm.                                         "<<endl;
    cout<<"**********************************************************************************************"<<endl;
    exit(-1);
  }else {
    for(i=0;i<numArgs;i++){
      if (strcmp(argv[2*i+1], "-numM")==0) {
	M=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-numN")==0) {
	N=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-codeName")==0) {
	strcpy(codeName, argv[2*i+2]); 
      } else if(strcmp(argv[2*i+1], "-degFileName")==0) {
	strcpy(degFileName, argv[2*i+2]); 
      } else if(strcmp(argv[2*i+1], "-sglConcent")==0) {
	sglConcent=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-tgtGirth")==0) {
	targetGirth=atoi(argv[2*i+2]);
      } else{
    goto USE;
      }
    }
    if(M>N) {
      cout<<"Warning: M must be samller than N"<<endl;
      exit(-1);
    }
  }

  degSeq=new int[N];

  ifstream infn(degFileName);
  if (!infn) {cout << "\nCannot open file " << degFileName << endl; exit(-1); } 
  infn >>m;
  deg=new int[m];
  degFrac=new double[m];
  for(i=0;i<m;i++) infn>>deg[i];
  for(i=0;i<m;i++) infn>>degFrac[i];
  infn.close();  
  double dtmp=0.0;
  for(i=0;i<m;i++) dtmp+=degFrac[i];
  cout.setf(ios::fixed, ios::floatfield);
  if(fabs(dtmp-1.0)>EPS) {
    cout.setf(ios::fixed, ios::floatfield);
    cout <<"\n Invalid degree distribution (node perspective): sum != 1.0 but "<<setprecision(10)<<dtmp<<endl; exit(-1); 
  } 
  for(i=1;i<m;i++) degFrac[i]+=degFrac[i-1];
  for(i=0;i<N;i++) {
    dtmp=(double)i/N;
    for(j=m-1;j>=0;j--) {
      if(dtmp>degFrac[j]) break;
    }
    if(dtmp<degFrac[0]) degSeq[i]=deg[0];
    else degSeq[i]=deg[j+1];
  }

  bigGirth=new BigGirth(M, N, degSeq, codeName, sglConcent, targetGirth);

  (*bigGirth).writeToFile_Hcompressed();
  //(*bigGirth).writeToFile_Hmatrix()        //  different output format
  //(*bigGirth).writeToFile();               //  different output format: including generator matrix (compressed)
  
  //computing local girth distribution  
  if(N<10000) {
    cout<<" Now computing the local girth on the global Tanner graph setting. "<<endl;
    cout<<"     might take a bit long time. Please wait ...                   "<<endl;
    (*bigGirth).loadH();
    cog=new CyclesOfGraph(M, N, (*bigGirth).H);
    (*cog).getCyclesTable();
    (*cog).printCyclesTable();
    delete cog;
    cog=NULL;
  }

  delete [] degSeq;  degSeq=NULL;
  delete [] deg; deg=NULL;
  delete [] degFrac; degFrac=NULL;
  delete bigGirth;
}




