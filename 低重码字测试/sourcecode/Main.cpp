/************************************************************************/
/*                                                                      */
/*        Free software: An approximate algorithm for computing MinDist */
/*                                                                      */
/*        Created by: Xiaoyu Hu                                         */
/*        IBM Research, Zurich Research Lab., Switzerland               */
/*                    and Marc P. C. Fossorier                          */
/*        Dept. Elect. Engineering, Univ. Hawaii at Manoa, USA          */
/*                                                                      */
/*        The C++ sources files have been compiled using xlC compiler   */
/*        at IBM RS/6000 running AIX. For other compilers and platforms,*/
/*        minor changes might be needed.                                */
/*                                                                      */
/*        Bug reporting to: xhu@zurich.ibm.com                          */
/************************************************************************/
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include "BlockCodeGenerator.h"
#include "BlockCodeGenerator.C"
#include "FactorGraph.h"
#include "FactorGraph.C"
#include "BlaumSpectra.h"
#include "BlaumSpectra.C"
#include "SoftDemodulator.h"
#include "SoftDemodulator.C"
#include "AWGN.h"
#include "AWGN.C"
#include "ChannelTransition.h"
#include "ChannelTransition.C"
#include "FunctionNode.h"
#include "FunctionNode.C"
#include "Node.h"
#include "Node.C"
#include "ParityCheck.h"
#include "ParityCheck.C"
#include "ParityCheckNode.h"
#include "ParityCheckNode.C"
#include "VariableNode.h"
#include "VariableNode.C"
int main(int argc, char * argv[]) {
  int i, j, k, l, s, absSyndrome, numIters;
  int *codeWord;
  int threshWeight=30, speedupFactor=1;
  double *errorImpulse, *(*probabilities), *codeWordLLR;

  char codeName[20], outNameCW[20], outNameSpectra[20];
  BlockCodeGenerator *blockCodeGenerator;
  BlaumSpectra       *blaum;
  FactorGraph        *factorGraph;
  SoftDemodulator    *softDemod;

  const int   maxIterations=80;
  const int   minIterations=1;

  int opt=1;

  int numArgs=(argc-1)/2;

  if (argc<7) { // 
  USAGE:
    cout<<"----------------------------------------------------------------------------------"<<endl;
    cout<<"-                              Usage Reminder:                                   -"<<endl;
    cout<<"-    Main -code codeName  -outFileCW outNameCW -thresh estLowWeight              -"<<endl;
    cout<<"-      Option: -switch opt  -speedup inum                                        -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"-              opt==0 --- Berrou's impulse error followed by                     -"<<endl;
    cout<<"-                         single-position bit-reversing followed by              -"<<endl;
    cout<<"-                         two-position bit reversing                             -"<<endl;
    cout<<"-              opt==1 --- single- and two-position bit reversing (DEFAULT)       -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"-              inum==1 --- normal speed at each iteration (DEFAULT)              -"<<endl;
    cout<<"-                    2 --- faster, but low-weight codewords may be missed        -"<<endl;
    cout<<"-                    3 --- fastest, but low-weight codewords may be missed       -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"- This free software is to compute approximately the minimum distance and the    -"<<endl;
    cout<<"- corresponding multiplicity of iteratively decodeable linear codes, for instance"<<endl;
    cout<<"- LDPC codes. Recall that the general MinDist problem of linear codes is NP-hard,"<<endl;
    cout<<"- the software thus produces an estimate and an upper bound of the minimum distance"<<endl;
    cout<<"- based on true (low-weight) codewords found by a fine-tuned local search. The   -"<<endl;
    cout<<"- algorithm proved to be powerful in the sense that it found 2 codewords of weight"<<endl;
    cout<<"- 20 for MacKay (3, 6)-regular (504, 252) code, 1 codeword of weight 34 for MacKay"<<endl;
    cout<<"- (3, 6)-regular (1008, 504) code, 66 codewords of weight 40 for the Margulis (p=11)"<<endl;
    cout<<"- code, 2184 codewords of weight 14 for the Ramanujan-Margulis (13, 5) code, and -"<<endl;
    cout<<"- 204 codewords of weight 24 for the Ramanujan-Margulis (17, 5) code.            -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"- The main progamm Main reads in an LDPC code from 'codeName'. The format is as  -"<<endl;
    cout<<"- follows: The first line contains the number block length, N. The second line   -"<<endl;
    cout<<"- defines the number of parity-checks, M. The third line defines the number of   -"<<endl;
    cout<<"- columns of the compressed parity-check. The following M lines are then the     -"<<endl;
    cout<<"- compressed parity-check matrix. Each of the M rows contains the indices (1 ... N)"<<endl;
    cout<<"- of 1's in the compressed row of parity-check matrix. If not all column entries -"<<endl;
    cout<<"- are used, the column is filled up with 0's.                                    -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"- One needs to supply the programm with an estimate of the minimum distance of the"<<endl;
    cout<<"- code by specifying -thresh estLowWeight. The specific value of estLowWeight    -"<<endl;
    cout<<"- does not affect the ability of the algorithm to capture low-weight codewords,  -"<<endl;
    cout<<"- but it prescribes the largest weight of which codewords (found) are writen into-"<<endl;
    cout<<"- outNameCW --- for a codewith a minimum distance of 20, no one is interested in -"<<endl;
    cout<<"- the codeword of weight over 50. Therefore, estLowWeight shall be slightly larger"<<endl;
    cout<<"- the true minimum distance. If estLowWeight is smaller than the lowest weight   -"<<endl;
    cout<<"- found by the programm, no codewords are written into the output codeword file, -"<<endl;
    cout<<"- but the program will tell you the lowest weight ever found. Note that the output"<<endl;
    cout<<"- codeword file 'outNameCW' is updated from time to time, one can terminate the  -"<<endl;
    cout<<"- program at any time when he think the lowest weight is satisfactory.           -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"- A log file called 'outNameCW.spect' records the lower-part spectrum of the code "<<endl;
    cout<<"- based on the codewords found so far. The posotion denotes the distance, and    -"<<endl;
    cout<<"- the number the multiplicity.                                                   -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"-                                                                                -"<<endl;
    cout<<"-                             version 1.1, by X. Y. Hu and M. Fossorier, 09/2003 -"<<endl;
    cout<<"-                                                         version 1.0, 3/08/2002 -"<<endl;
    cout<<"----------------------------------------------------------------------------------"<<endl;
    exit(-1);
  } else {
    for(i=0;i<numArgs;i++){
      if (strcmp(argv[2*i+1], "-code")==0) {
	strcpy(codeName, argv[2*i+2]);
      } else if (strcmp(argv[2*i+1], "-outFileCW")==0) {
	strcpy(outNameCW, argv[2*i+2]);
	strcpy(outNameSpectra, outNameCW);
	strcat(outNameSpectra, ".spect");
      } else if (strcmp(argv[2*i+1], "-thresh")==0) {
	threshWeight=atoi(argv[2*i+2]); 
      } else if (strcmp(argv[2*i+1], "-switch")==0) {
	opt=atoi(argv[2*i+2]);  
      } else if (strcmp(argv[2*i+1], "-speedup")==0) {
	speedupFactor=atoi(argv[2*i+2]);  
      } else{
	cout<<"Wrong usage!"<<endl;
	goto USAGE;
      }
    }
  }

  if(speedupFactor<=0) speedupFactor=1;
  else if(speedupFactor>3) speedupFactor=3;
  blockCodeGenerator=new BlockCodeGenerator(codeName); 
  factorGraph       =new FactorGraph(blockCodeGenerator);
  blaum             =new BlaumSpectra(blockCodeGenerator, threshWeight);
  (*blaum).setSpeedupFact(speedupFactor);
  softDemod = new SoftDemodulator(2.0, 2); // 2 means the number of alphabet {1, -1}

  errorImpulse=new double [(*blockCodeGenerator).N];
  for(i=0;i<(*blockCodeGenerator).N;i++) errorImpulse[i]=-1.0;
  probabilities = new double* [(*blockCodeGenerator).N];
  for (i = 0; i < (*blockCodeGenerator).N; i++) {
    probabilities[i] = new double [2];
  }  
  codeWord=new int[(*blockCodeGenerator).N];
  codeWordLLR=new double[(*blockCodeGenerator).N];

  if(opt!=0) goto SINGLE_BIT_REVERSE;

 ERROR_IMPULSE:
  cout<<"---------------------------------------------------"<<endl<<endl;
  cout<<"         Berrou's Impulse Error                    "<<endl;
  cout<<"---------------------------------------------------"<<endl;

  for(i=0;i<(*blockCodeGenerator).N;i++) {
    for(j=1;j<((*blockCodeGenerator).N-(*blockCodeGenerator).K+3);j++) {
      errorImpulse[i]=j-1.0; 
      (*softDemod).getProbabilities(errorImpulse, (*blockCodeGenerator).N, probabilities);
      (*factorGraph).initialize(probabilities); 
      numIters=0;
      do {
	(*factorGraph).LDPCupdate();
	(*factorGraph).getCodeWord(codeWord, codeWordLLR);
	absSyndrome = (*factorGraph).getAbsSyndrome(codeWord);
	numIters++;
      } while(((absSyndrome != 0) && (numIters <= maxIterations)) || numIters < minIterations );
      s=0;
      for(k=0;k<(*blockCodeGenerator).N;k++) {
	if(codeWord[k]!=0) {s=1;break;}
      }
      if(s!=0) {
	break;
      }
    }
    if(s==0) {cout<<"Error-impulse method invalid, divert to one bit-reverse method"<<endl;goto SINGLE_BIT_REVERSE;}

    (*blaum).initTmpDmin();
    (*factorGraph).initialize(probabilities); 
    numIters=0;      
    do {
      (*factorGraph).LDPCupdate();
      (*factorGraph).getCodeWord(codeWord, codeWordLLR);
      if((numIters%speedupFactor)==0) (*blaum).erasureDecoding(codeWord, codeWordLLR);
      absSyndrome = (*factorGraph).getAbsSyndrome(codeWord);
      numIters++;      
    } while(((absSyndrome != 0) && (numIters <= maxIterations)) || numIters < minIterations );
    cout<<i<<"-th bit of impulse error="<<j<<", local low-weight found="<<(*blaum).tmpMinWeight<<", min dis. so far="<<(*blaum).minWeight<<endl;
    (*blaum).writeToFile(outNameCW, outNameSpectra);
    errorImpulse[i]=-1.0; // restoring the previous tested bit
  }

 SINGLE_BIT_REVERSE:
  cout<<"---------------------------------------------------"<<endl<<endl;
  cout<<"         Single-position bit-reverse               "<<endl;
  cout<<"---------------------------------------------------"<<endl;

  for(i=0;i<(*blockCodeGenerator).N;i++) errorImpulse[i]=-1.0;
  for(i=0;i<(*blockCodeGenerator).N;i++) {
    errorImpulse[i]=1.0;
    (*softDemod).getProbabilities(errorImpulse, (*blockCodeGenerator).N, probabilities);
    (*factorGraph).initialize(probabilities); 
    (*blaum).initTmpDmin();
    numIters=0;      
    do {
      (*factorGraph).LDPCupdateExBitK(i); // i-th bit node not being updated in the SPA decoder
      (*factorGraph).getCodeWord(codeWord, codeWordLLR);
      if((numIters%speedupFactor)==0) (*blaum).erasureDecoding(codeWord, codeWordLLR);
      absSyndrome = (*factorGraph).getAbsSyndrome(codeWord);
      numIters++;
    } while(((absSyndrome != 0) && (numIters <= maxIterations)) || numIters < minIterations );      
    cout<<"Low-weight found by "<<i<<"-th bit reversing := "<<(*blaum).tmpMinWeight<<"; minimum low-weight found so far :="<<(*blaum).minWeight<<endl;
    (*blaum).writeToFile(outNameCW, outNameSpectra);
    errorImpulse[i]=-1.0; // restoring the previous tested bit
  }

 DOUBLE_BIT_REVERSE:  
  cout<<"---------------------------------------------------"<<endl<<endl;
  cout<<"           Bit Reverse (double positions)          "<<endl;
  cout<<"---------------------------------------------------"<<endl;
  for(i=0;i<(*blockCodeGenerator).N;i++) errorImpulse[i]=-1.0;  // initialize
  for(i=0;i<(*blockCodeGenerator).N-1;i++) {    
    errorImpulse[i]=1.0;
    for(s=i+1;s<(*blockCodeGenerator).N;s++) {
      errorImpulse[s]=1.0;
      (*softDemod).getProbabilities(errorImpulse, (*blockCodeGenerator).N, probabilities);
      (*factorGraph).initialize(probabilities); 
      (*blaum).initTmpDmin();
      numIters=0;      
      do {
	(*factorGraph).LDPCupdateExBitK(i, s); // (i, s)-th bit node not being updated in the SPA decoder
	(*factorGraph).getCodeWord(codeWord, codeWordLLR);
	if((numIters%speedupFactor)==0) (*blaum).erasureDecoding(codeWord, codeWordLLR);
	absSyndrome = (*factorGraph).getAbsSyndrome(codeWord);
	numIters++;
      } while(((absSyndrome != 0) && (numIters <= maxIterations)) || numIters < minIterations );
      errorImpulse[s]=-1.0;
      cout<<"Low-weight found by double bit-reverse at ( "<<i<<", "<<s<<" ) :="<<(*blaum).tmpMinWeight<< "; minimum low-weight found so far="<<(*blaum).minWeight<<endl;
      (*blaum).writeToFile(outNameCW, outNameSpectra);
    }
    errorImpulse[i]=-1.0; // restoring the previous tested bit
  }
return 0;
 END:
  delete [] errorImpulse;

  for(i=0;i<(*blockCodeGenerator).N; i++)
    delete [] probabilities[i];
  delete [] probabilities;
  delete [] codeWord;
  delete [] codeWordLLR;
  
  delete factorGraph;
  delete blaum;
  delete blockCodeGenerator;
  delete softDemod;
return 0;
}



