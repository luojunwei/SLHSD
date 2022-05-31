#include<cstdlib>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include <stdlib.h>
#include<cstring>
#include <getopt.h>
#include <sstream>

#include <sys/types.h>   
#include <dirent.h>
#include <sys/stat.h>

#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "scaffolding.h"

using namespace std;

int main(int argc, char** argv) {
	int maxSize = 50000;

	char* line = new char[maxSize];
	char* resultOutPutDirectory = new char[50];

	strcpy(resultOutPutDirectory, "./SLHSD-OUTPut-Directory/");

	long int readCount = 2000;
	int minAlignmentRevised = 150;
	int minAlignmentScore = 30;
	int contigLengthIgnore = 0;
	int minOverlapLength = 100;
	int overlapContigCount = 2;
	int readType = 1;
	int insertSize = 0;
	int readLength = 0;


	//contig
	char* contigSetFile = NULL;
	char* contigFile = NULL;

	//shortBamFile
	char* originalShortReadResult = NULL;
	char* rightShortBamFile = NULL;
	char* leftShortBamFile = NULL;
	char aligning = false;

	//longBamFile
	char* longBamFile = NULL;
	char* longFileName = NULL;
	long int readLengthCutOff = 1000;

	//scaffold
	int minPairedReadCount = 1;
	double ratio = 1;

	int ch = 0;
	while ((ch = getopt(argc, argv, "c:b:r:l:p:d:a:")) != -1) {
		switch (ch) {
		case 'c':contigSetFile = (char*)(optarg); break;
		case 'b':longBamFile = (char*)(optarg); break;
		case 'r':rightShortBamFile = (char*)(optarg); break;
		case 'l':leftShortBamFile = (char*)(optarg); break;
		case 'p':readType = atoi(optarg); break;		//1-pair-end reads; 2-mate-pair reads
		case 'd':insertSize = atoi(optarg); break;
		case 'a':readLength = atoi(optarg); break;
		default: break;
		}
	}

	if (opendir(resultOutPutDirectory) == NULL) {
		mkdir(resultOutPutDirectory, 0777);
	}

	int fileNameLen = strlen(resultOutPutDirectory);

	ContigSetHead* contigSetHead = GetContigSet(contigSetFile, contigLengthIgnore, insertSize);
	contigSetHead->contigLengthIgnore = contigLengthIgnore;	
	contigSetHead->minAlignmentScore = minAlignmentScore;	
	contigSetHead->minOverlapLength = minOverlapLength;		
	contigSetHead->overlapContigCount = overlapContigCount;	
	contigSetHead->minAlignmentRevised = minAlignmentRevised;
	long int contigCount = contigSetHead->contigCount;


	char* file1 = NULL;
	if (file1 == NULL) {
		file1 = (char*)malloc(sizeof(char) * 50);
		strcpy(file1, resultOutPutDirectory);
		strcat(file1, "/originalShortResult.fa");
	}

	char* file2 = NULL;
	if (file2 == NULL) {
		file2 = (char*)malloc(sizeof(char) * fileNameLen + 50);
		strcpy(file2, resultOutPutDirectory);
		strcat(file2, "/originalPathInLongRead.fa");
	}

	char* file4 = NULL;
	if (file4 == NULL) {
		file4 = (char*)malloc(sizeof(char) * fileNameLen + 50);
		strcpy(file4, resultOutPutDirectory);
		strcat(file4, "/originalPathInShortContig.fa");
	}

	char* file5 = (char*)malloc(sizeof(char) * fileNameLen + 50);
	strcpy(file5, resultOutPutDirectory);
	strcat(file5, "/optimaizePathInLongRead.fa");
	

	FILE* fp5;
	if ((fp5 = fopen(file5, "w")) == NULL) {
		printf("%s, does not exist!", file5);
		exit(0);
	}

	char* file6 = (char*)malloc(sizeof(char) * fileNameLen + 50);
	strcpy(file6, resultOutPutDirectory);
	strcat(file6, "/optimaizePathInShortContig.fa");


	FILE* fp6;
	if ((fp6 = fopen(file6, "w")) == NULL) {
		printf("%s, does not exist!", file6);
		exit(0);
	}
	
	ReadMapPosition* readMapPosition = GetReadInformation(contigSetHead, leftShortBamFile, rightShortBamFile, longBamFile, file1, file2, file4, insertSize, readType);
	
	OptimizeLongRead(contigSetHead, file2, line, maxSize, fp5);

	OptimizeLongRead(contigSetHead, file4, line, maxSize, fp6);
	fflush(fp5);
	fflush(fp6);

	fclose(fp5);
	fclose(fp6);
	
	char* file3 = (char*)malloc(sizeof(char) * 50);
	strcpy(file3, resultOutPutDirectory);
	strcat(file3, "/shortReadAlignResult.fa");

	AligningResultHead* aligningResultHead = GetAligningResultHead(contigSetHead, file3, readType, file1, insertSize, maxSize, line);

	SimpleResultHead* simpleResultHead = OptimaizeLongAlign(contigSetHead, file5);
	

	ScaffoldGraphHead* scaffoldGraphHead = (ScaffoldGraphHead*)malloc(sizeof(ScaffoldGraphHead));
	scaffoldGraphHead->scaffoldGraph = (ScaffoldGraph*)malloc(sizeof(ScaffoldGraph) * contigCount);
	scaffoldGraphHead->scaffoldGraphNodeCount = contigCount;
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		scaffoldGraphHead->scaffoldGraph[i].outEdge = NULL;
		scaffoldGraphHead->scaffoldGraph[i].inEdge = NULL;
	}
	LocalScaffoldSetHead* localScaffoldSetHead = NULL;

	switch (readType)
	{
	case 1:
		BuildGraphWithSmallIS(scaffoldGraphHead, contigSetHead, aligningResultHead, simpleResultHead, readMapPosition, readType, insertSize, localScaffoldSetHead, file5);
		OptimizeScaffoldGraph1(scaffoldGraphHead, contigSetHead, readType);

		break;
	case 2:
		GetScaffoldGraphHeadFromAlignResult(scaffoldGraphHead, contigSetHead, aligningResultHead, simpleResultHead, readMapPosition, readType, insertSize, localScaffoldSetHead);
		OptimizeScaffoldGraph(scaffoldGraphHead, contigSetHead, readType);

		break;
	default:
		cout << "please enter the read type in '1'(paired-end) or '2'(mate-pair)";
		break;
	}


	ScaffoldSetHead* scaffoldSetHead = GetScaffoldSet(scaffoldGraphHead, contigSetHead, simpleResultHead, line, maxSize, file6, insertSize);

	char* scaffoldTagFileName = (char*)malloc(sizeof(char) * fileNameLen + 50);
	strcpy(scaffoldTagFileName, resultOutPutDirectory);
	strcat(scaffoldTagFileName, "/scaffold");
	OutputScaffoldSet(scaffoldSetHead->scaffoldSet, contigSetHead, scaffoldTagFileName);

}


