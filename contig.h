#pragma once

#ifndef CONTIG_H_INCLUDED 
#define CONTIG_H_INCLUDED 


#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>


#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using namespace BamTools;
using namespace std;

typedef struct Contig {
	char* contig;
	long int contigLength;
	long int contigId;	
	char* name;		
	bool repeativeContig;
	bool shortContig;

	Contig() {
		contig = NULL;
		contigLength = 0;
		contigId = 0;
		name = NULL;
		repeativeContig = false;
		shortContig = false;
	}
}Contig;

typedef struct ContigSetHead {
	Contig* contigSet;
	long int contigCount;
	long int allContigLength;
	bool* repeatContigIndex;
	Contig* minContigSet;
	long int minContigCount;
	long int minAllContigLength;
	int contigLengthIgnore;
	int minAlignmentScore;
	int minOverlapLength;
	int overlapContigCount;
	int minAlignmentRevised;
	bool* visited;
	int* flag;
}ContigSetHead;




ContigSetHead* GetContigSet(char* contigSetFile, int contigLengthThreshold, int insertSize);
void outputContigSet(ContigSetHead* contigSetHead, int contigLengthThreshold);
//void DetectRepeativeContigInSet(ContigSetHead* contigSetHead, char* bamFileName, float ratio);




#endif