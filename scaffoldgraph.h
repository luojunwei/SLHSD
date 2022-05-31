#pragma once
#ifndef SCAFFOLDGRAPH_H_INCLUDED 
#define SCAFFOLDGRAPH_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "contig.h"
#include "aligningFromBam.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using namespace BamTools;
using namespace std;

//typedef struct EdgeInfo {
//	bool leftOritation;
//	bool rightOritation;
//	long int leftPosition;
//	long int rightPosition;
//
//}EdgeInfo;


typedef struct ScaffoldGraphEdge {
	int index;
	int nodeIndex;
	int allReadCount;
	int shortReadCount;
	int longReadCount;
	double weight;
	int orientation;
	double gapDistance;
	int overlapLength;
	ScaffoldGraphEdge* next;
	bool visited;
	int flag;

	bool* ori;
	bool indexOri;
	bool nodeIndexOri;
	ScaffoldGraphEdge() {
		index = -1;
		nodeIndex = -1;
		allReadCount = 0;
		shortReadCount = 0;
		longReadCount = 0;
		weight = 0;
		gapDistance = 0;
		overlapLength = 0;
		next = NULL;
		visited = false;
		//flag = 0;
		ori = NULL;
	}
}ScaffoldGraphEdge;

typedef struct ScaffoldGraph {
	ScaffoldGraphEdge* outEdge;
	ScaffoldGraphEdge* inEdge;
	long int edgeCount;
}ScaffoldGraph;

typedef struct ScaffoldGraphHead {
	ScaffoldGraph* scaffoldGraph;
	long int scaffoldGraphNodeCount;
}ScaffoldGraphHead;


typedef struct LocalScaffoldSet {
	int leftContigIndex;
	int righrContigIndex;
	int localScaffoldFlag;
	double localGapDistance;
	int localOverlap;

	bool leftOritation;
	bool rightOritation;

	long int leftPosition;
	long int rightPosition;

}LocalScaffoldSet;

typedef struct LocalScaffoldSetHead {
	long int localScaffoldNum;
	LocalScaffoldSet* localScaffoldSet;
}LocalScaffoldSetHead;


ScaffoldGraphHead* GetScaffoldGraphHeadFromAlignResult(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, AligningResultHead* aligningResultHead, SimpleResultHead* simpleResultHead, ReadMapPosition* readMapPosition, int readType, int insertSize, LocalScaffoldSetHead* localScaffoldSetHead);
void InsertEdge(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, int index, long int leftNodeIndex, long int rightNodeIndex, long int gapDistance, int overlapLength, int shortCount, int longCount, ReadMapPosition* readMapPosition, double w);

double GetEdgeWeight(ContigSetHead* contigSetHead, long int leftNodeIndex, long int rightNodeIndex, int shortOriIndex, int longOriIndex, int shortCount, int longCount, ReadMapPosition* readMapPosition, long int gapDistance, int readType, int insertSize);


int cmplocal(const void* arg1, const void* arg2);
LocalScaffoldSetHead* mixShortAndLong(ContigSetHead* contigSetHead, SimpleLongResult* simpleLongResult, AligningResult* aligningResult, long int shortLen, long int longLen);
ScaffoldGraphEdge* OptimizeEdge(ScaffoldGraphEdge* edge, ContigSetHead* contigSetHead, int nodeIndex);
ScaffoldGraphEdge* MergeEdges(ScaffoldGraphEdge* edge, int contigLength, int contigLength0);
void DeleteMergeEdge(ScaffoldGraphEdge* edge);
ScaffoldGraphEdge* MergeMultipleEdges(ScaffoldGraphEdge* edge, int contigLength, int contigLength0);
ScaffoldGraphEdge* OptimizeEdgeInScaffoldGraph(ScaffoldGraphEdge* edge, ContigSetHead* contigSetHead, int nodeIndex);
void OptimizeScaffoldGraph1(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, int readType);
void OptimizeScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, int readType);

ScaffoldGraphEdge* DeleteEdgeNode(ScaffoldGraphEdge* edge, int nodeIndex, bool orientation);
void DeleteMinEdgeWeight(ScaffoldGraphHead* scaffoldGraphHead, int min);
int CountAverageReadCount(ScaffoldGraphHead* scaffoldGraphHead);
int DeleteSpecialScaffoldEdge(ScaffoldGraph* scaffoldGraph, long int index, long int index1);
void RemoveCycleInScaffold(ScaffoldGraphHead* scaffoldGraphHead);
void DeleteMinEdgeWeight1(ScaffoldGraphHead* scaffoldGraphHead, int minReadCount);
void OutputScaffoldGraph(ScaffoldGraphHead* scaffoldGraphHead);
void OutputScaffold(ScaffoldGraph* scaffoldGraph, int nodeCount);

void InsertOutOrInEdge(ScaffoldGraphHead* scaffoldGraphHead, int readIndex, int leftNodeIndex, bool leftOrientation, int rightNodeIndex, bool rightOrientation, int gapDistance, int overlapLength);
void BuildGraphWithSmallIS(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, AligningResultHead* aligningResultHead, SimpleResultHead* simpleResultHead, ReadMapPosition* readMapPosition, int readType, int insertSize, LocalScaffoldSetHead* localScaffoldSetHead, char* file);

ScaffoldGraphHead* BuildGraph(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, AligningResultHead* aligningResultHead, SimpleResultHead* simpleResultHead, ReadMapPosition* readMapPosition, int readType, int insertSize, LocalScaffoldSetHead* localScaffoldSetHead, char* file);
#endif