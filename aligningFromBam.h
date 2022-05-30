#ifndef aligningFromBam_H_INCLUDED 
#define aligningFromBam_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "contig.h"

using namespace std;
using namespace BamTools;


typedef struct ReadMapPosition {
	long int leftShortCoverage;
	long int rightShortCoverage;
	long int* leftShortIndex;
	long int* rightShortIndex;

	long int leftLongCoverage;
	long int rightLongCoverage;
	long int* leftLongIndex;
	long int* rightLongIndex;
	ReadMapPosition(){
		leftShortCoverage = 0;
		rightShortCoverage = 0;
		leftShortIndex = NULL;
		rightShortIndex = NULL;
		leftLongCoverage = 0;
		rightLongCoverage = 0;
		leftLongIndex = NULL;
		rightLongIndex = NULL;
	}

}ReadMapPosition;

typedef struct AligningResult {  
	long int leftPosition;  
	long int rightPosition;	
	long int leftContigIndex;
	long int rightContigIndex;
	bool leftOrientation;	
	bool rightOrientation;
	long int leftLength;
	long int rightLength;
	double gapDiatanceForSR;

}AligningResult;

typedef struct AligningResultHead {
	AligningResult* aligningResult;
	int allocateAligningResultCount;
	long int aligningResultCount;
}AligningResultHead;

typedef struct LongResult {
	long int readStartPosition;
	long int readEndPosition;
	long int contigStartPosition;
	long int contigEndPosition;

	int orientation;
	int orientation1;
	int contigIndex;
	int contigIndex1;
	int overlapLength;
	int quality;
	double gapdistace;
}LongResult;

typedef struct LongResultHead {
	LongResult* longResult;
	int allocateLongResultCount;
	int longResultCount;
	int aligningShortContigResultCount;
}LongResultHead;

//----------------------------------------
typedef struct SimpleLongResult {
	int leftIndex;
	int rightIndex;
	int simpleAlignCount;
	int count1;
	int count2;
	int count3;
	int count4;
	double gapSimpleLong;
	int overlapSimpleLong;
	int* oriLongCount;

	bool leftOritation;
	bool rightOritation;
	long int leftPosition;
	long int rightPosition;
}SimpleLongResult;

typedef struct SimpleShortResult {
	int leftContigIndex;
	int rightContigIndex;
	int simpleShortAlignCount;
	int count1;
	int count2;
	int count3;
	int count4;
	int flag;
	double gapSimpleShort;
	int overlapSimpleShort;
	int* oriShortCount;
}SimpleShortResult;


typedef struct SimpleResultHead {
	SimpleShortResult* simpleShortResult;
	SimpleLongResult* simpleLongResult;
	long int simpleShortLineCount;
	long int simpleLongLineCount;
}SimpleResult;


ReadMapPosition* GetReadInformation(ContigSetHead* contigSetHead, char* leftShortBamFile, char* rightShortBamFile, char* longBamFile, char* originalShortResult, char* originalLongResult, char* originalShortContig, int insertSize, int readType);
bool GetAlignLongResultOneLine(LongResultHead* longResultHead, BamAlignment alignment, ContigSetHead* contigSetHead, long int index);
void outPutLongAlignResult(LongResultHead* longResultHead, ContigSetHead* contigSetHead, FILE* fp, FILE* fp1, long int readIndex, long int readLength);

AligningResultHead* GetAligningResultHead(ContigSetHead* contigSetHead, char* file, int readType, char* originalShortResult, int insertSize, int maxSize, char* line);
int cmp(const void* arg1, const void* arg2);
void OptimaizeShortAlignResult(AligningResultHead* aligningResultHead);
double CaculateShortDistance(AligningResultHead* aligningResultHead, SimpleResultHead* simpleResultHead, long int& count, long int leftIndex, long int rightIndex, long int startIndex, int* type, int max);
void OptimizeLongRead(ContigSetHead* contigSetHead, char* originalfile, char* line, int maxSize, FILE* fp1);

void AddEdgeInSimpleGraph(SimpleResultHead* simpleResultHead, long int simpleLongAllCount, int leftIndex, bool leftOrientation, int rightIndex, int rightOrientation, int distance, int minOverlapLength);
void OptimaizeLongAlignResult(ContigSetHead* contigSetHead, char* originalLongResult, SimpleResultHead* simpleResultHead);
int cmpl(const void* arg1, const void* arg2);
SimpleResultHead* deleteDuplicateInSimpleLong(SimpleResultHead* simpleResultHead);
void mixOrientation(SimpleResultHead* simpleResultHead, int start, int end);

SimpleLongResult* GetLineIndex1(ContigSetHead* contigSetHead, bool* lineIndex, char* file, char* line, int maxSize, long int& lineCount, long int& simpleGraphNodeCount);
void AddEdgeInSimpleGraph1(SimpleLongResult* simpleLongResult, long int simpleLongAllCount, int leftIndex, bool leftOrientation, int rightIndex, int rightOrientation, int distance, int minOverlapLength);

SimpleResultHead* OptimaizeLongAlign(ContigSetHead* contigSetHead, char* optimizeLongResult);


#endif