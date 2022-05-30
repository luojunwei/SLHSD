#ifndef BUILDSCAFFOLDGRAPH_CPP_INCLUDED 
#define BUILDSCAFFOLDGRAPH_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <exception>
#include<malloc.h>

#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"

using namespace std;



int cmplocal(const void* arg1, const void* arg2) {
	LocalScaffoldSet* a1 = (LocalScaffoldSet*)arg1;
	LocalScaffoldSet* a2 = (LocalScaffoldSet*)arg2;

	if (a1->leftContigIndex != a2->leftContigIndex) {
		return a1->leftContigIndex - a2->leftContigIndex;
	}
	else {
		return a1->righrContigIndex - a2->righrContigIndex;
	}
}


LocalScaffoldSetHead* mixShortAndLong(ContigSetHead* contigSetHead,  SimpleLongResult* simpleLongResult, AligningResult* aligningResult, long int shortLen, long int longLen) {
	LocalScaffoldSetHead* localScaffoldSetHead = (LocalScaffoldSetHead*)malloc(sizeof(LocalScaffoldSetHead));
	localScaffoldSetHead->localScaffoldNum = shortLen + longLen;
	cout << shortLen << "," << longLen << endl;

	localScaffoldSetHead->localScaffoldSet = (LocalScaffoldSet*)malloc(sizeof(LocalScaffoldSet) * localScaffoldSetHead->localScaffoldNum);
	for (long int i = 0; i < localScaffoldSetHead->localScaffoldNum; i++) {
		localScaffoldSetHead->localScaffoldSet[i].leftContigIndex = -1;
		localScaffoldSetHead->localScaffoldSet[i].righrContigIndex = -1;
		localScaffoldSetHead->localScaffoldSet[i].localScaffoldFlag = 0;
		localScaffoldSetHead->localScaffoldSet[i].localGapDistance = 0;
		localScaffoldSetHead->localScaffoldSet[i].localOverlap = 0;
		localScaffoldSetHead->localScaffoldSet[i].leftOritation = false;
		localScaffoldSetHead->localScaffoldSet[i].rightOritation = false;
		localScaffoldSetHead->localScaffoldSet[i].leftPosition = 0;
		localScaffoldSetHead->localScaffoldSet[i].rightPosition = 0;
	}

	long int i=0;
	long int j=0;
	long int nCount = 0;
	for (i = 0; i < shortLen; i++) {
		localScaffoldSetHead->localScaffoldSet[nCount].leftContigIndex = aligningResult[i].leftContigIndex;
		localScaffoldSetHead->localScaffoldSet[nCount].righrContigIndex = aligningResult[i].rightContigIndex;
		localScaffoldSetHead->localScaffoldSet[nCount].localScaffoldFlag = 1;
		localScaffoldSetHead->localScaffoldSet[nCount].leftOritation = aligningResult[i].leftOrientation;
		localScaffoldSetHead->localScaffoldSet[nCount].rightOritation = aligningResult[i].rightOrientation;
		localScaffoldSetHead->localScaffoldSet[nCount].localGapDistance = aligningResult[i].gapDiatanceForSR;
		localScaffoldSetHead->localScaffoldSet[nCount].localOverlap = 0;
		localScaffoldSetHead->localScaffoldSet[nCount].leftPosition = aligningResult[i].leftPosition;
		localScaffoldSetHead->localScaffoldSet[nCount].rightPosition = aligningResult[i].rightPosition;

		nCount++;
	}

	for (j = 0; j < longLen; j++) {
		localScaffoldSetHead->localScaffoldSet[nCount].leftContigIndex = simpleLongResult[j].leftIndex;
		localScaffoldSetHead->localScaffoldSet[nCount].righrContigIndex = simpleLongResult[j].rightIndex;
		localScaffoldSetHead->localScaffoldSet[nCount].localScaffoldFlag = 2;
		localScaffoldSetHead->localScaffoldSet[nCount].leftOritation = simpleLongResult[j].leftOritation;
		localScaffoldSetHead->localScaffoldSet[nCount].rightOritation = simpleLongResult[j].rightOritation;

		localScaffoldSetHead->localScaffoldSet[nCount].localGapDistance = simpleLongResult[j].gapSimpleLong;
		localScaffoldSetHead->localScaffoldSet[nCount].localOverlap = simpleLongResult[j].overlapSimpleLong;

		localScaffoldSetHead->localScaffoldSet[nCount].leftPosition = simpleLongResult[j].leftPosition;
		localScaffoldSetHead->localScaffoldSet[nCount].rightPosition = simpleLongResult[j].rightPosition;
		nCount++;
	}

	localScaffoldSetHead->localScaffoldNum = nCount;

	qsort(localScaffoldSetHead->localScaffoldSet, localScaffoldSetHead->localScaffoldNum, sizeof(localScaffoldSetHead->localScaffoldSet[0]), cmplocal);

	return localScaffoldSetHead;

}


ScaffoldGraphHead* GetScaffoldGraphHeadFromAlignResult(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, AligningResultHead* aligningResultHead, SimpleResultHead* simpleResultHead, ReadMapPosition* readMapPosition, int readType, int insertSize) {
	
	LocalScaffoldSetHead* localScaffoldSetHead = NULL;
	long int shortLen = 0;
	long int longLen = 0;
	if (insertSize > 1000) {
		shortLen = aligningResultHead->aligningResultCount;
		longLen = simpleResultHead->simpleLongLineCount;
		localScaffoldSetHead = mixShortAndLong(contigSetHead, simpleResultHead->simpleLongResult, aligningResultHead->aligningResult, shortLen, longLen);
	}
	else {
		shortLen = 0;
		longLen = simpleResultHead->simpleLongLineCount;
		localScaffoldSetHead = mixShortAndLong(contigSetHead, simpleResultHead->simpleLongResult, aligningResultHead->aligningResult, shortLen, longLen);
	}
	long int i = 0;

	long int leftContigIndex = -1;
	long int rightContigIndex = -1;
	long int start = 0;
	long int end = -1;
	long int m = 0;

	long int shortCount = 0;
	long int longCount = 0;
	double gapDistance = 0;
	long int overlap = 0;

	int* so = (int*)malloc(sizeof(int) * 4);
	int* lo = (int*)malloc(sizeof(int) * 4);
	double* sd = (double*)malloc(sizeof(double) * 4);
	double* ld = (double*)malloc(sizeof(double) * 4);
	int* lp = (int*)malloc(sizeof(int) * 4);


	for (int u = 0; u < 4; u++) {
		so[u] = 0;
		lo[u] = 0;
		sd[u] = 0;
		ld[u] = 0;
		lp[u] = 0;
	}

	int max = 0;
	int maxIndex = -1;
	int second = 0;
	int secondIndex = -1;

	int max0 = 0;
	int maxIndex0 = -1;
	int second0 = 0;
	int secondIndex0 = -1;

	bool leftOrientation;
	bool rightOrientation;

	int index = 0;
	double weight = 0;


	for (long int i = 0; i < localScaffoldSetHead->localScaffoldNum + 1; i++) {
		if (localScaffoldSetHead->localScaffoldSet[i].leftContigIndex == leftContigIndex && localScaffoldSetHead->localScaffoldSet[i].righrContigIndex == rightContigIndex) {
			end++;
		}
		else {
			if (end - start + 1 > 0) {
				for (m = start; m <= end; m++) {
					gapDistance = gapDistance + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
					if (overlap < localScaffoldSetHead->localScaffoldSet[m].localOverlap) {
						overlap = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
					} 

					if (localScaffoldSetHead->localScaffoldSet[m].localScaffoldFlag == 1) {
						shortCount++;
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 0 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 0) {
							so[0]++;
							sd[0] = sd[0] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 0 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 1) {
							so[1]++;
							sd[1] = sd[1] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 1 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 0) {
							so[2]++;
							sd[2] = sd[2] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 1 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 1) {
							so[3]++;
							sd[3] = sd[3] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
						}
					}
					if (localScaffoldSetHead->localScaffoldSet[m].localScaffoldFlag == 2 && readType == 2) {
						longCount++;
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 0 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 1) {
							if (localScaffoldSetHead->localScaffoldSet[m].leftPosition <= localScaffoldSetHead->localScaffoldSet[m].rightPosition) {
								lo[3]++;
								ld[3] = ld[3] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[3]) {
									lp[3] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}
							else {
								lo[0]++;
								ld[0] = ld[0] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[0]) {
									lp[0] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}

						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 0 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 0) {
							lo[2]++;
							ld[2] = ld[2] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
							if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[2]) {
								lp[2] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
							}
						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 1 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 1) {
							lo[1]++;
							ld[1] = ld[1] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
							if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[1]) {
								lp[1] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
							}

						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 1 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 0) {
							if (localScaffoldSetHead->localScaffoldSet[m].leftPosition <= localScaffoldSetHead->localScaffoldSet[m].rightPosition) {
								lo[0]++;
								ld[0] = ld[0] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[0]) {
									lp[0] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}
							else {
								lo[3]++;
								ld[3] = ld[3] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[3]) {
									lp[3] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}
						}
					}
					else if (localScaffoldSetHead->localScaffoldSet[m].localScaffoldFlag == 2 && readType == 1) {
						longCount++;
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 0 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 1) {
							if (localScaffoldSetHead->localScaffoldSet[m].leftPosition <= localScaffoldSetHead->localScaffoldSet[m].rightPosition) {
								lo[0]++;
								ld[0] = ld[0] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[0]) {
									lp[0] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}
							else {
								lo[3]++;
								ld[3] = ld[3] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[3]) {
									lp[3] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}

						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 0 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 0) {
							lo[1]++;
							ld[1] = ld[1] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
							if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[1]) {
								lp[1] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
							}


						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 1 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 1) {
							lo[2]++;
							ld[2] = ld[2] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
							if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[2]) {
								lp[2] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
							}

						}
						if (localScaffoldSetHead->localScaffoldSet[m].leftOritation == 1 && localScaffoldSetHead->localScaffoldSet[m].rightOritation == 0) {
							if (localScaffoldSetHead->localScaffoldSet[m].leftPosition <= localScaffoldSetHead->localScaffoldSet[m].rightPosition) {
								lo[3]++;
								ld[3] = ld[3] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[3]) {
									lp[3] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}
							else {
								lo[0]++;
								ld[0] = ld[0] + localScaffoldSetHead->localScaffoldSet[m].localGapDistance;
								if (localScaffoldSetHead->localScaffoldSet[m].localOverlap > lp[0]) {
									lp[0] = localScaffoldSetHead->localScaffoldSet[m].localOverlap;
								}
							}

						}
					}
				}
				for (int u = 0; u < 4; u++) {
					if (so[u] != 0) {
						sd[u] = sd[u] / (double)(so[u]);
					}
					if (lo[u] != 0) {
						ld[u] = ld[u] / (double)(lo[u]);
					}
				}
				gapDistance = gapDistance / (double)(end - start + 1);
				
				max = 0;
				maxIndex = -1;
				second = 0;
				secondIndex = 0;
				max0 = 0;
				maxIndex0 = -1;
				second0 = 0;
				secondIndex0 = 0;
				for (int p = 0; p < 4; p++) {
					if (so[p] > max) {
						max = so[p];
						maxIndex = p;
					}
				}
				for (int p = 0; p < 4; p++) {
					if (so[p] > second && p != maxIndex) {
						second = so[p];
						secondIndex = p;
					}
				}
				for (int p = 0; p < 4; p++) {
					if (lo[p] > max0) {
						max0 = lo[p];
						maxIndex0 = p;
					}
				}
				for (int p = 0; p < 4; p++) {
					if (lo[p] > second0 && p != maxIndex0) {
						second0 = lo[p];
						secondIndex0 = p;
					}
				}
				weight = 0;
								
				if (shortCount >= 10 && longCount >= 1) {

					if (second == 0 && second0 == 0) {	
						if (maxIndex == maxIndex0) {
							index = maxIndex0;
							gapDistance = (sd[maxIndex] + ld[maxIndex0]) / 2;
							weight = GetEdgeWeight(contigSetHead, leftContigIndex, rightContigIndex, maxIndex, maxIndex0, shortCount, longCount, readMapPosition, gapDistance, readType, insertSize);
						}
						else {
							if (max0 >= 1 && lp[maxIndex0] >= 1000) {	
								index = maxIndex0;
								shortCount = 0;
								gapDistance = ld[index];
								weight = GetEdgeWeight(contigSetHead, leftContigIndex, rightContigIndex, -1, index, shortCount, longCount, readMapPosition, gapDistance, readType, insertSize);
							}
							else {
								index = maxIndex;
								longCount = 0;
								gapDistance = sd[index];
								weight = GetEdgeWeight(contigSetHead, leftContigIndex, rightContigIndex, index, -1, shortCount, longCount, readMapPosition, gapDistance, readType, insertSize);
							}
						}
						
					}
					else {
						if (max0 == second0) {
							if (maxIndex == maxIndex0) {
								index = maxIndex0;
								longCount = longCount - second0;
							}
							if (maxIndex == secondIndex0) {
								index = secondIndex0;
								longCount = longCount - max0;
							}
						}
						else {
							if (maxIndex != maxIndex0 && maxIndex == secondIndex0 && max + second0 >= max0) {
								index = secondIndex0;
								longCount = longCount - max0;
							}
							else {
								index = maxIndex0;
								longCount = longCount - second0;
							}
						}
						gapDistance = (sd[index] + ld[index]) / 2;
						weight = GetEdgeWeight(contigSetHead, leftContigIndex, rightContigIndex, maxIndex, index, shortCount, longCount, readMapPosition, gapDistance, readType, insertSize);
					}
					InsertEdge(scaffoldGraphHead, contigSetHead ,index, leftContigIndex, rightContigIndex, gapDistance, overlap, shortCount, longCount, readMapPosition, weight);
				}

				if (shortCount >= 30 && longCount == 0) {
					if (max == second) {
						bool a = abs(sd[maxIndex]) < abs(sd[secondIndex]) ? true : false;
						if (a == true) {
							index = maxIndex;
							shortCount = shortCount - second;
						}
						else {
							index = secondIndex;
							shortCount = shortCount - max;
						}
					}
					else {
						index = maxIndex;
						shortCount = shortCount - second;
					}
					gapDistance = sd[index];
					weight = GetEdgeWeight(contigSetHead, leftContigIndex, rightContigIndex, index, -1, shortCount, longCount, readMapPosition, sd[index], readType, insertSize);
					InsertEdge(scaffoldGraphHead, contigSetHead, index, leftContigIndex, rightContigIndex, gapDistance, overlap, shortCount, longCount, readMapPosition, weight);

				}

				if (shortCount == 0 && longCount >= 3) {
					if (max0 == second0) {
						bool b = lp[maxIndex0] > lp[secondIndex0] ? true : false;
						if (b == true) {
							index = maxIndex0;
							longCount = longCount - second0;
						}
						else {
							index = secondIndex0;
							longCount = longCount - max0;
						}
					}
					else {
						index = maxIndex0;
						longCount = longCount - second0;
					}
					gapDistance = ld[index];
					weight = GetEdgeWeight(contigSetHead, leftContigIndex, rightContigIndex, -1, index, shortCount, longCount, readMapPosition, ld[index], readType, insertSize);
					InsertEdge(scaffoldGraphHead, contigSetHead, index, leftContigIndex, rightContigIndex, gapDistance, overlap, shortCount, longCount, readMapPosition, weight);
				}
				
			}
			start = i;
			end = i;
			leftContigIndex = localScaffoldSetHead->localScaffoldSet[i].leftContigIndex;
			rightContigIndex = localScaffoldSetHead->localScaffoldSet[i].righrContigIndex;
			shortCount = 0;
			longCount = 0;
			for (int u = 0; u < 4; u++) {
				so[u] = 0;
				lo[u] = 0;
				sd[u] = 0;
				ld[u] = 0;
				lp[u] = 0;
			}
			gapDistance = 0;
			weight = 0;
		}
	}
	


	return scaffoldGraphHead;
}


void InsertEdge(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, int index, long int leftNodeIndex, long int rightNodeIndex, long int gapDistance, int overlapLength, int shortCount, int longCount, ReadMapPosition* readMapPosition, double w) {

	bool leftOrientation;
	bool rightOrientation;

	ScaffoldGraphEdge* leftNode = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
	leftNode->index = leftNodeIndex;
	leftNode->nodeIndex = rightNodeIndex;
	leftNode->gapDistance = gapDistance;
	leftNode->overlapLength = overlapLength;
	leftNode->allReadCount = shortCount + longCount;
	leftNode->shortReadCount = shortCount;
	leftNode->longReadCount = longCount;
	leftNode->weight = w; 
	leftNode->next = NULL;

	ScaffoldGraphEdge* rightNode = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
	rightNode->index = rightNodeIndex;
	rightNode->nodeIndex = leftNodeIndex;
	rightNode->gapDistance = gapDistance;
	rightNode->overlapLength = overlapLength;
	rightNode->allReadCount = shortCount + longCount;
	rightNode->shortReadCount = shortCount;
	rightNode->longReadCount = longCount;
	rightNode->weight = w; 
	rightNode->next = NULL;


	if (index == 0 || index == 3) {
		leftNode->orientation = 0;
		rightNode->orientation = 0;
		if (index == 0) {
			leftOrientation = true;
			rightOrientation = false;
		}
		if (index == 3) {
			leftOrientation = false;
			rightOrientation = true;
		}
	}
	else {
		leftNode->orientation = 1;
		rightNode->orientation = 1;
		if (index == 1) {
			leftOrientation = true;
			rightOrientation = true;
		}
		if (index == 2) {
			leftOrientation = false;
			rightOrientation = false;
		}
	}

	ScaffoldGraphEdge* temp = NULL;
	ScaffoldGraphEdge* temp1 = NULL;
	if (leftOrientation == true) {
		if (scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge != NULL) {
			leftNode->next = scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge;
		}
		scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge = leftNode;
	}
	else {
		if (scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge != NULL) {
			leftNode->next = scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge;
		}
		scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge = leftNode;
	}
	if (rightOrientation == true) {
		if (scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge != NULL) {
			rightNode->next = scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge;
		}
		scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge = rightNode;
	}
	else {
		if (scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge != NULL) {
			rightNode->next = scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge;
		}
		scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge = rightNode;
	}

	

}


double GetEdgeWeight(ContigSetHead* contigSetHead, long int leftNodeIndex, long int rightNodeIndex, int shortOriIndex, int longOriIndex, int shortCount, int longCount, ReadMapPosition* readMapPosition, long int gapDistance, int readType, int insertSize) {
	int tempLength = 0;

	double shortCov = 0;
	double longCov = 0;
	long int* leftCount = new long int[4];
	long int* rightCount = new long int[4];
	for (long int i = 0; i < 4; i++) {	
		leftCount[i] = 0;
		rightCount[i] = 0;
	}

	if (contigSetHead->contigSet[leftNodeIndex].contigLength >= insertSize * 1.5) {
		tempLength = insertSize * 1.5 - abs(gapDistance);
		for (int i = 0; i < tempLength; i++) {
			leftCount[0] = leftCount[0] + readMapPosition[leftNodeIndex].leftShortIndex[i];
			leftCount[2] = leftCount[2] + readMapPosition[leftNodeIndex].leftLongIndex[i];
		}
		for (int j = contigSetHead->contigSet[leftNodeIndex].contigLength; j > contigSetHead->contigSet[leftNodeIndex].contigLength - tempLength; j--) {
			leftCount[1] = leftCount[1] + readMapPosition[leftNodeIndex].rightShortIndex[j];
			leftCount[3] = leftCount[3] + readMapPosition[leftNodeIndex].rightLongIndex[j];
		}
	}
	else {
		tempLength = contigSetHead->contigSet[leftNodeIndex].contigLength * 0.7;
		for (int i = 0; i < tempLength; i++) {
			leftCount[0] = leftCount[0] + readMapPosition[leftNodeIndex].leftShortIndex[i];
			leftCount[2] = leftCount[2] + readMapPosition[leftNodeIndex].leftLongIndex[i];
		}
		for (int j = contigSetHead->contigSet[leftNodeIndex].contigLength; j > contigSetHead->contigSet[leftNodeIndex].contigLength - tempLength; j--) {
			leftCount[1] = leftCount[1] + readMapPosition[leftNodeIndex].rightShortIndex[j];
			leftCount[3] = leftCount[3] + readMapPosition[leftNodeIndex].rightLongIndex[j];
		}
	}
	if (contigSetHead->contigSet[rightNodeIndex].contigLength >= insertSize * 1.5) {
		tempLength = insertSize * 1.5 - abs(gapDistance);
		for (int i = 0; i < tempLength; i++) {
			rightCount[0] = rightCount[0] + readMapPosition[rightNodeIndex].leftShortIndex[i];
			rightCount[2] = rightCount[2] + readMapPosition[rightNodeIndex].leftLongIndex[i];
		}
		for (int j = contigSetHead->contigSet[rightNodeIndex].contigLength; j > contigSetHead->contigSet[rightNodeIndex].contigLength - tempLength; j--) {
			rightCount[1] = rightCount[1] + readMapPosition[rightNodeIndex].rightShortIndex[j];
			rightCount[3] = rightCount[3] + readMapPosition[rightNodeIndex].rightLongIndex[j];
		}
	}
	else {
		tempLength = contigSetHead->contigSet[rightNodeIndex].contigLength * 0.7;
		for (int i = 0; i < tempLength; i++) {
			rightCount[0] = rightCount[0] + readMapPosition[rightNodeIndex].leftShortIndex[i];
			rightCount[2] = rightCount[2] + readMapPosition[rightNodeIndex].leftLongIndex[i];
		}
		for (int j = contigSetHead->contigSet[rightNodeIndex].contigLength; j > contigSetHead->contigSet[rightNodeIndex].contigLength - tempLength; j--) {
			rightCount[1] = rightCount[1] + readMapPosition[rightNodeIndex].rightShortIndex[j];
			rightCount[3] = rightCount[3] + readMapPosition[rightNodeIndex].rightLongIndex[j];
		}
	}


	for (int i = 0; i < contigSetHead->contigSet[leftNodeIndex].contigLength; i++) {
		leftCount[2] = leftCount[2] + readMapPosition[leftNodeIndex].leftLongIndex[i]; 
		leftCount[3] = leftCount[3] + readMapPosition[leftNodeIndex].rightLongIndex[i];
	}
	for (int i = 0; i < contigSetHead->contigSet[rightNodeIndex].contigLength; i++) {
		rightCount[2] = rightCount[2] + readMapPosition[rightNodeIndex].leftLongIndex[i];
		rightCount[3] = rightCount[3] + readMapPosition[rightNodeIndex].rightLongIndex[i];
	}


	if (shortOriIndex != -1 && readType == 1) {
		if (shortOriIndex == 2) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[1]) * ((double)shortCount / (double)rightCount[0]));
		}
		if (shortOriIndex == 3) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[1]) * ((double)shortCount / (double)rightCount[1]));
		}
		if (shortOriIndex == 1) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[0]) * ((double)shortCount / (double)rightCount[1]));
		}
		if (shortOriIndex == 0) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[0]) * ((double)shortCount / (double)rightCount[0]));
		}
	}
	else if (shortOriIndex != -1 && readType == 2) {
		if (shortOriIndex == 1) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[1]) * ((double)shortCount / (double)rightCount[0]));
		}
		if (shortOriIndex == 0) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[1]) * ((double)shortCount / (double)rightCount[1]));
		}
		if (shortOriIndex == 2) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[0]) * ((double)shortCount / (double)rightCount[1]));
		}
		if (shortOriIndex == 3) {
			shortCov = sqrt(((double)shortCount / (double)leftCount[0]) * ((double)shortCount / (double)rightCount[0]));
		}
	}
	else {
		shortCov = 0;
	}



	if (longOriIndex != -1 && readType == 1) {
		if (longOriIndex == 0) {
			longCov = sqrt(((double)longCount / (double)leftCount[3]) * ((double)longCount / (double)rightCount[2]));
		}
		if (longOriIndex == 1) {
			longCov = sqrt(((double)longCount / (double)leftCount[3]) * ((double)longCount / (double)rightCount[3]));
		}
		if (longOriIndex == 2) {
			longCov = sqrt(((double)longCount / (double)leftCount[2]) * ((double)longCount / (double)rightCount[2]));
		}
		if (longOriIndex == 3) {
			longCov = sqrt(((double)longCount / (double)leftCount[2]) * ((double)longCount / (double)rightCount[3]));
		}
	}
	else if (longOriIndex != -1 && readType == 2) {
		if (longOriIndex == 0) {
			longCov = sqrt(((double)longCount / (double)leftCount[2]) * ((double)longCount / (double)rightCount[3]));
		}
		if (longOriIndex == 1) {
			longCov = sqrt(((double)longCount / (double)leftCount[2]) * ((double)longCount / (double)rightCount[2]));
		}
		if (longOriIndex == 2) {
			longCov = sqrt(((double)longCount / (double)leftCount[3]) * ((double)longCount / (double)rightCount[3]));
		}
		if (longOriIndex == 3) {
			longCov = sqrt(((double)longCount / (double)leftCount[3]) * ((double)longCount / (double)rightCount[2]));
		}
	}
	else {
		longCov = 0;
	}

	delete[] leftCount;
	delete[] rightCount;


	if (shortCov == 0) {
		return longCov;
	}
	else if (longCov == 0) {
		return shortCov;
	}
	else {
		return ((shortCov + longCov) / 2);
	}
}


void DeleteMergeEdge(ScaffoldGraphEdge* edge) {
	ScaffoldGraphEdge* tempEdge = edge;
	while (edge != NULL) {
		tempEdge = edge->next;
		edge->next = NULL;
		free(edge);
		edge = tempEdge;
	}
}

ScaffoldGraphEdge* MergeMultipleEdges(ScaffoldGraphEdge* edge, int contigLength, int contigLength0) {

	int forwardCount1 = 0;
	double forwardGapDistance1 = 0;
	int forwardOverlapLength = 0;
	int forwardShortCount = 0;
	int forwardLongCount = 0;
	double forwardWeight = 0;

	int reverseCount1 = 0;
	double reverseGapDistance1 = 0;
	int reverseOverlapLength = 0;
	int reverseShortCount = 0;
	int reverseLongCount = 0;
	double reverseWeight = 0;

	ScaffoldGraphEdge* result = NULL;
	ScaffoldGraphEdge* tempEdge = edge;

	bool orientation = tempEdge->orientation;
	bool previous;

	while (tempEdge != NULL) {
		if (tempEdge->orientation == 1) {
			forwardCount1++;
			forwardGapDistance1 = forwardGapDistance1 + tempEdge->gapDistance;
			forwardShortCount = forwardShortCount + tempEdge->shortReadCount;
			forwardLongCount = forwardLongCount + tempEdge->longReadCount;
			forwardWeight = forwardWeight + tempEdge->weight;
			if (tempEdge->overlapLength > forwardOverlapLength) {
				forwardOverlapLength = tempEdge->overlapLength;
			}
		}
		else {
			reverseCount1++;
			reverseGapDistance1 = reverseGapDistance1 + tempEdge->gapDistance;
			reverseShortCount = reverseShortCount + tempEdge->shortReadCount;
			reverseLongCount = reverseLongCount + tempEdge->longReadCount;
			reverseWeight = reverseWeight + tempEdge->weight;
			if (tempEdge->overlapLength > reverseOverlapLength) {
				reverseOverlapLength = tempEdge->overlapLength;
			}
		}
		tempEdge = tempEdge->next;
	}
	if (forwardCount1 > 0) {
		forwardGapDistance1 = forwardGapDistance1 / forwardCount1;

	}
	if (reverseCount1 > 0) {
		reverseGapDistance1 = reverseGapDistance1 / reverseCount1;
	}


	if (forwardCount1 > 0) {
		ScaffoldGraphEdge* result1 = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
		result1->index = edge->index;
		result1->nodeIndex = edge->nodeIndex;
		result1->gapDistance = forwardGapDistance1;
		result1->allReadCount = forwardShortCount + forwardLongCount;
		result1->shortReadCount = forwardShortCount;
		result1->longReadCount = forwardLongCount;
		result1->orientation = true;
		result1->weight = forwardWeight;
		result1->next = NULL;
		result = result1;
	}

	if (reverseCount1 > 0) {
		ScaffoldGraphEdge* result1 = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
		result1->index = edge->index;
		result1->nodeIndex = edge->nodeIndex;
		result1->gapDistance = reverseGapDistance1;
		result1->allReadCount = reverseShortCount + reverseLongCount;
		result1->shortReadCount = reverseShortCount;
		result1->longReadCount = reverseLongCount;
		result1->orientation = false;
		result1->weight = reverseWeight;
		result1->next = NULL;
		result = result1;
	}
	DeleteMergeEdge(edge);
	return result;
}


ScaffoldGraphEdge* MergeEdges(ScaffoldGraphEdge* edge, int contigLength, int contigLength0) {

	int forwardCount1 = 0;
	double forwardGapDistance1 = 0;
	int forwardOverlapLength = 0;

	int reverseCount1 = 0;
	double reverseGapDistance1 = 0;
	int reverseOverlapLength = 0;

	ScaffoldGraphEdge* result = NULL;
	ScaffoldGraphEdge* tempEdge = edge;

	bool orientation = tempEdge->orientation;
	bool previous;

	while (tempEdge != NULL) {
		if (tempEdge->orientation == 1) {
			forwardCount1++;
			forwardGapDistance1 = forwardGapDistance1 + tempEdge->gapDistance;
			if (tempEdge->overlapLength > forwardOverlapLength) {
				forwardOverlapLength = tempEdge->overlapLength;
			}
		}
		else {
			reverseCount1++;
			reverseGapDistance1 = reverseGapDistance1 + tempEdge->gapDistance;
			if (tempEdge->overlapLength > reverseOverlapLength) {
				reverseOverlapLength = tempEdge->overlapLength;
			}
		}
		tempEdge = tempEdge->next;
	}
	if (forwardCount1 > 0) {
		forwardGapDistance1 = forwardGapDistance1 / forwardCount1;

	}
	if (reverseCount1 > 0) {
		reverseGapDistance1 = reverseGapDistance1 / reverseCount1;
	}


	if (forwardCount1 > 0) {
		ScaffoldGraphEdge* result1 = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
		result1->index = edge->index;
		result1->nodeIndex = edge->nodeIndex;
		result1->gapDistance = forwardGapDistance1;
		result1->allReadCount = forwardCount1;
		result1->shortReadCount = 0;
		result1->longReadCount = forwardCount1;
		result1->orientation = true;
		result1->weight = forwardOverlapLength;
		result1->next = NULL;
		forwardCount1 = 0;
		result = result1;
	}

	if (reverseCount1 > 0) {
		ScaffoldGraphEdge* result1 = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
		result1->index = edge->index;
		result1->nodeIndex = edge->nodeIndex;
		result1->gapDistance = reverseGapDistance1;
		result1->allReadCount = reverseCount1;
		result1->shortReadCount = 0;
		result1->longReadCount = reverseCount1;
		result1->orientation = false;
		result1->weight = reverseOverlapLength;
		result1->next = NULL;
		reverseCount1 = 0;
		result = result1;
	}
	DeleteMergeEdge(edge);
	return result;
}



ScaffoldGraphEdge* OptimizeEdgeInScaffoldGraph(ScaffoldGraphEdge* edge, ContigSetHead* contigSetHead, int nodeIndex) {
	ScaffoldGraphEdge* previousEdge = NULL;
	ScaffoldGraphEdge* previousEdge1 = NULL;
	ScaffoldGraphEdge* result = NULL;
	ScaffoldGraphEdge* tempStart = edge;

	while (tempStart != NULL) {
		ScaffoldGraphEdge* tempEdge = tempStart;
		tempStart = NULL;
		previousEdge1 = NULL;
		ScaffoldGraphEdge* tempEdge1 = tempEdge->next;
		previousEdge = tempEdge;
		bool a = false;
		while (tempEdge1 != NULL) {
			if (tempEdge1->nodeIndex == tempEdge->nodeIndex) {
				previousEdge->next = tempEdge1;
				previousEdge = tempEdge1;
			}
			else if (a != true) {
				tempStart = tempEdge1;
				previousEdge1 = tempStart;
				a = true;
			}
			else {
				previousEdge1->next = tempEdge1;
				previousEdge1 = tempEdge1;
			}
			tempEdge1 = tempEdge1->next;
		}
		if (previousEdge1 != NULL) {
			previousEdge1->next = NULL;
		}
		if (previousEdge != NULL) {
			previousEdge->next = NULL;
		}


		tempEdge = MergeMultipleEdges(tempEdge, contigSetHead->contigSet[nodeIndex].contigLength, contigSetHead->contigSet[tempEdge->nodeIndex].contigLength);

		if (tempEdge != NULL) {
			if (result != NULL) {
				ScaffoldGraphEdge* tempEdge2 = tempEdge->next;
				while (tempEdge2 != NULL) {
					if (tempEdge2->next == NULL) {
						break;
					}
				}
				if (tempEdge2 != NULL) {
					tempEdge2->next = result;
				}
				else {
					tempEdge->next = result;
				}
			}
			result = tempEdge;
		}

	}
	return result;
}



ScaffoldGraphEdge* OptimizeEdge(ScaffoldGraphEdge* edge, ContigSetHead* contigSetHead, int nodeIndex) {
	ScaffoldGraphEdge* previousEdge = NULL;
	ScaffoldGraphEdge* previousEdge1 = NULL;
	ScaffoldGraphEdge* result = NULL;
	ScaffoldGraphEdge* tempStart = edge;

	while (tempStart != NULL) {
		ScaffoldGraphEdge* tempEdge = tempStart;
		tempStart = NULL;
		previousEdge1 = NULL;
		ScaffoldGraphEdge* tempEdge1 = tempEdge->next;
		previousEdge = tempEdge;
		bool a = false;
		while (tempEdge1 != NULL) {
			if (tempEdge1->nodeIndex == tempEdge->nodeIndex) {
				previousEdge->next = tempEdge1;
				previousEdge = tempEdge1;
			}
			else if (a != true) {
				tempStart = tempEdge1;
				previousEdge1 = tempStart;
				a = true;
			}
			else {
				previousEdge1->next = tempEdge1;
				previousEdge1 = tempEdge1;
			}
			tempEdge1 = tempEdge1->next;
		}
		if (previousEdge1 != NULL) {
			previousEdge1->next = NULL;
		}
		if (previousEdge != NULL) {
			previousEdge->next = NULL;
		}


		tempEdge = MergeEdges(tempEdge, contigSetHead->contigSet[nodeIndex].contigLength, contigSetHead->contigSet[tempEdge->nodeIndex].contigLength);

		if (tempEdge != NULL) {
			if (result != NULL) {
				ScaffoldGraphEdge* tempEdge2 = tempEdge->next;
				while (tempEdge2 != NULL) {
					if (tempEdge2->next == NULL) {
						break;
					}
				}
				if (tempEdge2 != NULL) {
					tempEdge2->next = result;
				}
				else {
					tempEdge->next = result;
				}
			}
			result = tempEdge;
		}

	}
	return result;
}



int CountAverageReadCount(ScaffoldGraphHead* scaffoldGraphHead) {
	ScaffoldGraphEdge* tempEdge = NULL;
	ScaffoldGraphEdge* tempEdge1 = NULL;
	int averageShortCount = 0;
	int averageLongCount = 0;
	int averageCount = 0;
	int edgeCount = 0;

	for (int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {

		if (scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL) {
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			while (tempEdge != NULL) {
				averageShortCount = tempEdge->shortReadCount + averageShortCount;
				averageLongCount = tempEdge->longReadCount + averageLongCount;
				averageCount = tempEdge->allReadCount + averageCount;
				edgeCount++;
				tempEdge = tempEdge->next;
			}
		}
		if (scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL) {
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			while (tempEdge != NULL) {
				averageShortCount = tempEdge->shortReadCount + averageShortCount;
				averageLongCount = tempEdge->longReadCount + averageLongCount;
				averageCount = tempEdge->allReadCount + averageCount;
				edgeCount++;
				tempEdge = tempEdge->next;
			}
		}
	}
	return averageCount / edgeCount;
}


void OptimizeScaffoldGraph1(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead) {

	for (int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		if (scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL) {
			scaffoldGraphHead->scaffoldGraph[i].outEdge = OptimizeEdgeInScaffoldGraph(scaffoldGraphHead->scaffoldGraph[i].outEdge, contigSetHead, i);
		}

		if (scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL) {
			scaffoldGraphHead->scaffoldGraph[i].inEdge = OptimizeEdgeInScaffoldGraph(scaffoldGraphHead->scaffoldGraph[i].inEdge, contigSetHead, i);
		}
	}
	RemoveCycleInScaffold(scaffoldGraphHead);


	ScaffoldGraphEdge* tempEdge = NULL;
	int averageCount = 0;
	int edgeCount = 0;

	//DeleteMinEdgeWeight(scaffoldGraphHead, min);
	//OutputScaffoldGraph(scaffoldGraphHead);
	//cout << "------------------------------" << endl;
}


 void RemoveCycleInScaffold(ScaffoldGraphHead* scaffoldGraphHead){
 	ScaffoldGraphEdge* tempEdge=NULL;
 	ScaffoldGraphEdge* tempEdge1=NULL;
 	ScaffoldGraphEdge* tempEdgeNext=NULL;
 	for(int i=0;i<scaffoldGraphHead->scaffoldGraphNodeCount;i++){
 		tempEdge=scaffoldGraphHead->scaffoldGraph[i].outEdge;
 		while(tempEdge!=NULL){
 			tempEdgeNext=tempEdge->next;
 			tempEdge1=tempEdge->next;
 			bool del=false;
 			while(tempEdge1!=NULL){
 				del=false;
 				if(tempEdge->nodeIndex==tempEdge1->nodeIndex){
 					int nodeIndex=tempEdge1->nodeIndex;
 					bool orientation=tempEdge->orientation;
 					bool orientation1=tempEdge1->orientation;
 					if(tempEdge->weight > tempEdge1->weight){
 						scaffoldGraphHead->scaffoldGraph[i].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge,tempEdge1->nodeIndex,tempEdge1->orientation);
 						if(orientation1==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation1);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation1);
 						}
 					}
 					else if(tempEdge->weight < tempEdge1->weight){
 						scaffoldGraphHead->scaffoldGraph[i].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge,tempEdge->nodeIndex,tempEdge->orientation);
 						if(orientation==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation);
 						}
 					}
 					else{
 						scaffoldGraphHead->scaffoldGraph[i].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge,tempEdge1->nodeIndex,tempEdge1->orientation);
 						if(orientation1==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation1);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation1);
 						}
 						scaffoldGraphHead->scaffoldGraph[i].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge,tempEdge->nodeIndex,tempEdge->orientation);
 						if(orientation==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation);
 						}
 					}
 					del=true;
 					break;
 				}
 				tempEdge1=tempEdge1->next;
 			}
 			if(del!=false){
 				tempEdge=scaffoldGraphHead->scaffoldGraph[i].outEdge;
 			}
 			else{
 				tempEdge=tempEdge->next;
 			}
 		}
 		tempEdge=scaffoldGraphHead->scaffoldGraph[i].inEdge;
 		while(tempEdge!=NULL){
 			tempEdgeNext=tempEdge->next;
 			tempEdge1=tempEdge->next;
 			bool del=false;
 			while(tempEdge1!=NULL){
 				del=false;
 				if(tempEdge->nodeIndex==tempEdge1->nodeIndex){
 					int nodeIndex=tempEdge1->nodeIndex;
 					bool orientation=tempEdge->orientation;
 					bool orientation1=tempEdge1->orientation;
 					if(tempEdge->weight > tempEdge1->weight){
 						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge,tempEdge1->nodeIndex,tempEdge1->orientation);
 						if(orientation1==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation1);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation1);
 						}
 					}
 					else if(tempEdge->weight < tempEdge1->weight){
 						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge,tempEdge->nodeIndex,tempEdge->orientation);
 						if(orientation==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation);
 						}
 					}
 					else{
 						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge,tempEdge1->nodeIndex,tempEdge1->orientation);
 						if(orientation==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation1);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation1);
 						}
 						scaffoldGraphHead->scaffoldGraph[i].inEdge= DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge,tempEdge->nodeIndex,tempEdge->orientation);
 						if(orientation==true){
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge,i,orientation);
 						}else{
 							scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge,i,orientation);
 						}
 					}
 					del=true;
 					break;
 				}
 				tempEdge1=tempEdge1->next;
 			}
 			if(del!=false){
 				tempEdge=scaffoldGraphHead->scaffoldGraph[i].inEdge;
 			}
 			else{
 				tempEdge=tempEdge->next;
 			}
 		}
 	}
 }


ScaffoldGraphEdge* DeleteEdgeNode(ScaffoldGraphEdge* edge, int nodeIndex, bool orientation) {
	ScaffoldGraphEdge* firstEdge = edge;
	ScaffoldGraphEdge* previousEdge = NULL;

	while (edge != NULL) {
		if (edge->nodeIndex == nodeIndex && edge->orientation == orientation) {
			if (previousEdge == NULL) {
				firstEdge = edge->next;
			}
			else {
				previousEdge->next = edge->next;
			}
			edge->next = NULL;
			free(edge);
			break;
		}
		previousEdge = edge;
		edge = edge->next;
	}
	return firstEdge;
}


void DeleteMinEdgeWeight(ScaffoldGraphHead* scaffoldGraphHead, int min) {
	ScaffoldGraphEdge* tempedge = NULL;
	ScaffoldGraphEdge* tempEdge1 = NULL;

	for (int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		ScaffoldGraphEdge* temp = NULL;
		if (scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL) {
			tempedge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			while (tempedge != NULL) {
				int nodeIndex = tempedge->nodeIndex;
				bool orientation = tempedge->orientation;

				if (tempedge->weight == 0) {
					tempEdge1 = tempedge->next;
					if (orientation == true) {
						temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge;
						if (temp != NULL) {
							if (temp->weight == 0) {
								scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempedge->nodeIndex, orientation);
								scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
							}
						}
						else {
							scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempedge->nodeIndex, orientation);
						}
					}
					else{
						temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge;
						if (temp != NULL) {
							if (temp->weight == 0) {
								scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempedge->nodeIndex, orientation);
								scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
							}
						}
						else {
							scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, tempedge->nodeIndex, orientation);
						}
					}
					tempedge = tempEdge1;
					continue;
				}

				tempedge = tempedge->next;
			}
		}

		if (scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL) {
			tempedge = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			while (tempedge != NULL) {
				int nodeIndex = tempedge->nodeIndex;
				bool orientation = tempedge->orientation;

				if (tempedge->weight == 0) {
					tempEdge1 = tempedge->next;
					if (orientation == true) {
						temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge;
						if (temp != NULL) {
							if (temp->weight == 0) {
								scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempedge->nodeIndex, orientation);
								scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
							}
						}
						else {
							scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempedge->nodeIndex, orientation);
						}
					}
					else {
						temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge;
						if (temp != NULL) {
							if (temp->weight == 0) {
								scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempedge->nodeIndex, orientation);
								scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
							}
						}
						else {
							scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, tempedge->nodeIndex, orientation);
						}
					}
					tempedge = tempEdge1;
					continue;
				}

				tempedge = tempedge->next;
			}
		}
	}

}



int DeleteSpecialScaffoldEdge(ScaffoldGraph* scaffoldGraph, long int index, long int index1) {
	scaffoldGraph[index].outEdge = DeleteEdgeNode(scaffoldGraph[index].outEdge, index1, 0);
	scaffoldGraph[index].outEdge = DeleteEdgeNode(scaffoldGraph[index].outEdge, index1, 1);
	scaffoldGraph[index].inEdge = DeleteEdgeNode(scaffoldGraph[index].inEdge, index1, 0);
	scaffoldGraph[index].inEdge = DeleteEdgeNode(scaffoldGraph[index].inEdge, index1, 1);
}

void BuildGraphWithSmallIS(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, AligningResultHead* aligningResultHead, SimpleResultHead* simpleResultHead, ReadMapPosition* readMapPosition, int readType, int insertSize, char* file) {
	cout << "readType=1" << endl;
	long int lineCount = 0;
	FILE* fp;
	if ((fp = fopen(file, "r")) == NULL) {
		printf("%s, does not exist!", file);
		exit(0);
	}
	char* p;
	const char* split = ",";
	int maxSize = 50000;
	char* line = new char[maxSize];

	cout << "gggg" << endl;
	int readIndex = 0;
	lineCount = 0;
	while ((fgets(line, maxSize, fp)) != NULL) {
		p = strtok(line, split);
		int count = atoi(p);
		if (count <= 1) {	
			continue;
		}
		p = strtok(NULL, split);
		int startContigIndex = atoi(p);
		p = strtok(NULL, split);
		int startPosition = atoi(p);
		p = strtok(NULL, split);
		int distance = atoi(p);
		p = strtok(NULL, split);
		bool startOrientation = atoi(p);
		p = strtok(NULL, split);
		int startOverlapLength = atoi(p);

		int a = 2;
		while (a <= count) {
			p = strtok(NULL, split);
			int endContigIndex = atoi(p);
			p = strtok(NULL, split);
			int endPosition = atoi(p);
			p = strtok(NULL, split);
			int distance1 = atoi(p);
			p = strtok(NULL, split);
			bool endOrientation = atoi(p);
			p = strtok(NULL, split);
			int endOverlapLength = atoi(p);

			if (startContigIndex != -1) {
				if (endContigIndex == -1) {
					distance = distance + distance1 + contigSetHead->contigSet[endContigIndex].contigLength;
				}
				else {
					int min = 0;
					if (startOverlapLength > endOverlapLength) {
						min = endOverlapLength;
					}
					else {
						min = startOverlapLength;
					}
					InsertOutOrInEdge(scaffoldGraphHead, readIndex, startContigIndex, startOrientation, endContigIndex, endOrientation, distance, min);

					distance = distance1;
					startContigIndex = endContigIndex;
					startOrientation = endOrientation;
					startOverlapLength = endOverlapLength;
					startPosition = endPosition;
				}
			}
			else {
				distance = distance1;
				startContigIndex = endContigIndex;
				startOrientation = endOrientation;
				startOverlapLength = endOverlapLength;
				startPosition = endPosition;
			}
			a++;
		}
		readIndex++;
	}
	fclose(fp);


}

void InsertOutOrInEdge(ScaffoldGraphHead* scaffoldGraphHead, int readIndex, int leftNodeIndex, bool leftOrientation, int rightNodeIndex, bool rightOrientation, int gapDistance, int overlapLength) {

	ScaffoldGraphEdge* leftNode = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
	leftNode->index = leftNodeIndex;
	leftNode->nodeIndex = rightNodeIndex;
	leftNode->gapDistance = gapDistance;
	leftNode->overlapLength = overlapLength;
	leftNode->allReadCount = 0;
	leftNode->shortReadCount = 0;
	leftNode->longReadCount = overlapLength;
	leftNode->weight = 0;
	leftNode->next = NULL;

	ScaffoldGraphEdge* rightNode = (ScaffoldGraphEdge*)malloc(sizeof(ScaffoldGraphEdge));
	rightNode->index = rightNodeIndex;
	rightNode->nodeIndex = leftNodeIndex;
	rightNode->gapDistance = gapDistance;
	rightNode->overlapLength = overlapLength;
	rightNode->allReadCount = 0;
	rightNode->shortReadCount = 0;
	rightNode->longReadCount = overlapLength;
	rightNode->weight = 0;
	rightNode->next = NULL;


	if (leftOrientation != rightOrientation) {
		leftNode->orientation = 0;
		rightNode->orientation = 0;
	}
	else {
		leftNode->orientation = 1;
		rightNode->orientation = 1;
	}

	if (leftOrientation == true) {
		if (scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge != NULL) {
			leftNode->next = scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge;
		}
		scaffoldGraphHead->scaffoldGraph[leftNodeIndex].outEdge = leftNode;
	}
	else {
		if (scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge != NULL) {
			leftNode->next = scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge;
		}
		scaffoldGraphHead->scaffoldGraph[leftNodeIndex].inEdge = leftNode;
	}

	if (rightOrientation == true) {
		if (scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge != NULL) {
			rightNode->next = scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge;
		}
		scaffoldGraphHead->scaffoldGraph[rightNodeIndex].inEdge = rightNode;
	}
	else {
		if (scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge != NULL) {
			rightNode->next = scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge;
		}
		scaffoldGraphHead->scaffoldGraph[rightNodeIndex].outEdge = rightNode;
	}

}

void OptimizeGraph(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead) {

	for (int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		if (scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL) {
			scaffoldGraphHead->scaffoldGraph[i].outEdge = OptimizeEdge(scaffoldGraphHead->scaffoldGraph[i].outEdge, contigSetHead, i);
		}

		if (scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL) {
			scaffoldGraphHead->scaffoldGraph[i].inEdge = OptimizeEdge(scaffoldGraphHead->scaffoldGraph[i].inEdge, contigSetHead, i);
		}
	}

	RemoveCycleInScaffold(scaffoldGraphHead);

}




#endif