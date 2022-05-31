#pragma once
#ifndef SCAFFOLDING_H_INCLUDED 
#define SCAFFOLDING_H_INCLUDED 
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
#include "scaffoldgraph.h"
#include "lp/lp_lib.h"

using namespace std;

typedef struct ContigSequence {
    int index;
    int orientation;
    int gapDistance;
    ContigSequence* next;
}ContigSequence;

typedef struct ScaffoldSet {
    int length;
    ContigSequence* contigSequence;
    int* contigIndex;
    int contigNum;
    int sequenceCount;
    ScaffoldSet* next;
}ScaffoldSet;

typedef struct ScaffoldSetHead {
    ScaffoldSet* scaffoldSet;
}ScaffoldSetHead;



long int DetermineOrientationOfContigs(ScaffoldGraph* scaffoldGraph, long int contigCount, bool* contigOrientation);
long int* IterativeDetermineOrderOfContigs(ContigSetHead* contigSetHead, ScaffoldGraph* scaffoldGraph, long int contigCount, bool* contigOrientation, long int* contigPosition, long int allContigLength);

ScaffoldSetHead* GetScaffoldSet(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, SimpleResultHead* simpleResultHead, char* line, int maxSize, char* file, int insertSize);

void SortNode(int* sortNode, int* length, long int left, long int right);
ScaffoldGraphEdge* GetOptimizeNodeIndex(ScaffoldGraphHead* scaffoldGraphHead, int nodeIndex, bool orientation, ContigSequence* tempContigSequence, bool right, bool* printIndex);
int SearchScaffoldEdge(int index, ContigSequence* contigSequence);
void OptimizeScaffoldSetCongtigSequence(ScaffoldSet* scaffoldSet, int contigNum);
int GetContigSequenceNum(ContigSequence* tempContigSequence);

char* ReverseComplement(char* temp);
void OutputScaffoldSet(ScaffoldSet* scaffoldSet, ContigSetHead* contigSetHead, char* result);

bool SearchInsertContigSequence(ScaffoldSetHead* tempInsertSequenceHead, int startIndex, int endIndex, int* contigIndex, int* distance, bool* orientation, int* overlapLength, int count);
void InsertContigToSequence(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, bool* printIndex, char* file, char* line, int maxSize);
void InsertShortContigToSequence(ContigSequence* tempContigSequence, ContigSequence* preContigSequence, char* file, char* line, int maxSize, bool* printIndex, long int lineCount);



#endif