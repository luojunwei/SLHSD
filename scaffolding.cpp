#ifndef SCAFFOLDING_CPP_INCLUDED 
#define SCAFFOLDING_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <typeinfo>

#include "scaffoldgraph.h"
#include "scaffolding.h"
#include "lp/lp_lib.h"



long int DetermineOrientationOfContigs(ScaffoldGraph* scaffoldGraph, long int contigCount, bool* contigOrientation) {

	long int i = 0;
	long int j = 0;
	long int p = 1;
	long int c = 1000;

	bool* contigIndex = new bool[contigCount];

	long int edgeNumber = 0;
	long int edgeNumberout = 0;
	long int edgeNumberin = 0;
	long int constraintNumber = 0;
	ScaffoldGraphEdge* tempEdge = NULL;
	ScaffoldGraphEdge* tempEdge1 = NULL;

	for (int i = 0; i < contigCount; i++) {
		tempEdge = scaffoldGraph[i].outEdge;
		while (tempEdge != NULL) {
			tempEdge->visited = false;
			tempEdge = tempEdge->next;
		}
		tempEdge = scaffoldGraph[i].inEdge;
		while (tempEdge != NULL) {
			tempEdge->visited = false;
			tempEdge = tempEdge->next;
		}
	}
	for (i = 0; i < contigCount; i++) {
		tempEdge = scaffoldGraph[i].outEdge;
		while (tempEdge != NULL) {
			edgeNumber++;
			edgeNumberout++;
			tempEdge = tempEdge->next;
		}
		tempEdge = scaffoldGraph[i].inEdge;
		while (tempEdge != NULL) {
			edgeNumber++;
			edgeNumberin++;
			tempEdge = tempEdge->next;
		}
		contigIndex[i] = false;
	}

	edgeNumber = (edgeNumber / 2);

	constraintNumber = 0;

	lprec* lp;
	int Ncol, * colno = NULL, ret = 0;
	REAL* row = NULL;

	Ncol = edgeNumber + contigCount;
	lp = make_lp(0, Ncol);
	if (lp == NULL) {
		printf("couldn't construct a new linear programming model");
		exit(0);
	}
	double* weight = new double[edgeNumber];
	long int* edgeLeftNode = new long int[edgeNumber];
	long int* edgeRightNode = new long int[edgeNumber];

	colno = (int*)malloc(Ncol * sizeof(*colno));
	row = (REAL*)malloc(Ncol * sizeof(*row));
	if ((colno == NULL) || (row == NULL)) {
		printf("couldn't new colno and row");
		exit(0);
	}
	set_add_rowmode(lp, TRUE);
	int ttt = 0;

	for (i = 0; i < contigCount; i++) {
		for (int q = 0; q < 2; q++) {
			if (q == 0) {
				tempEdge = scaffoldGraph[i].outEdge;
			}
			else {
				tempEdge = scaffoldGraph[i].inEdge;
			}
			while (tempEdge != NULL) {
				if (tempEdge->visited == false) {
					if (tempEdge->orientation == false) {
						j = 0;
						colno[j] = i + 1;
						row[j++] = 1;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = 1;

						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c + 1)) {
							printf("couldn't add_constraintex");
							exit(0);
						}

						j = 0;
						colno[j] = i + 1;
						row[j++] = -1;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = -1;

						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c - 1)) {
							printf("couldn't add_constraintex");
							exit(0);
						}

						ttt++;
						constraintNumber = constraintNumber + 2;

					}
					else {

						j = 0;
						colno[j] = i + 1;
						row[j++] = 1;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = -1;

						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c)) {
							printf("couldn't add_constraintex");
							exit(0);
						}

						j = 0;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = 1;
						colno[j] = i + 1;
						row[j++] = -1;

						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c)) {
							printf("couldn't add_constraintex");
							exit(0);
						}

						constraintNumber = constraintNumber + 2;
						ttt++;
					}


					j = 0;
					colno[j] = contigCount + p;
					row[j++] = 1;
					add_constraintex(lp, j, row, colno, LE, 1);
					j = 0;
					colno[j] = contigCount + p;
					row[j++] = 1;
					add_constraintex(lp, j, row, colno, GE, 0);
					constraintNumber = constraintNumber + 2;

					edgeLeftNode[p - 1] = i;
					edgeRightNode[p - 1] = tempEdge->nodeIndex;
					p++;
					weight[p - 2] = tempEdge->weight;

					tempEdge->visited = true;
					tempEdge1 = scaffoldGraph[tempEdge->nodeIndex].outEdge;
					while (tempEdge1 != NULL) {
						if (tempEdge1->nodeIndex == i) {
							tempEdge1->visited = true;
							break;
						}
						tempEdge1 = tempEdge1->next;
					}
					tempEdge1 = scaffoldGraph[tempEdge->nodeIndex].inEdge;
					while (tempEdge1 != NULL) {
						if (tempEdge1->nodeIndex == i) {
							tempEdge1->visited = true;
							break;
						}
						tempEdge1 = tempEdge1->next;
					}

				}
				tempEdge = tempEdge->next;
			}
		}
	}

	for (i = 0; i < contigCount; i++) {
		set_binary(lp, i + 1, TRUE);
	}

	p--;

	j = 0;
	for (i = 0; i < p; i++) {
		colno[j] = contigCount + i + 1;
		row[j] = weight[j];
		j++;
	}

	if (!set_obj_fnex(lp, j, row, colno)) {
		printf("couldn't set_obj_fnex");
		exit(0);
	}
	set_add_rowmode(lp, FALSE);


	set_timeout(lp, 1200);

	set_maxim(lp);
	set_scaling(lp, 128);


	ret = solve(lp);	//ÕâÀï

	if (!(ret == 0 || ret == 1)) {
		cout << "ee:" << ret << endl;
	}
	REAL* pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
	get_primal_solution(lp, pv);

	double temp = 1;
	long int result = 0;
	for (i = contigCount + constraintNumber + 1; i < constraintNumber + contigCount + p + 1; i++) {
		if (pv[i] == 0) {
			DeleteSpecialScaffoldEdge(scaffoldGraph, edgeLeftNode[i - contigCount - constraintNumber - 1], edgeRightNode[i - contigCount - constraintNumber - 1]);
			DeleteSpecialScaffoldEdge(scaffoldGraph, edgeRightNode[i - contigCount - constraintNumber - 1], edgeLeftNode[i - contigCount - constraintNumber - 1]);
			result++;

		}

	}

	for (i = constraintNumber + 1; i < contigCount + constraintNumber + 1; i++) {
		contigOrientation[i - constraintNumber - 1] = pv[i];
	}

	delete[] contigIndex;
	delete[] weight;
	delete[] edgeLeftNode;
	delete[] edgeRightNode;
	delete[] pv;

	delete_lp(lp);

	return result;

}

long int* IterativeDetermineOrderOfContigs(ContigSetHead* contigSetHead, ScaffoldGraph* scaffoldGraph, long int contigCount, bool* contigOrientation, long int* contigPosition, long int allContigLength) {

	long int i = 0;
	long int j = 0;
	long int p = 1;
	long int c = 2 * allContigLength;

	bool* contigVisited = new bool[contigCount];

	for (i = 0; i < contigCount; i++) {
		contigVisited[i] = false;
	}

	long int edgeNumber = 0;
	long int edgeNumberOut = 0;
	long int edgeNumberIn = 0;
	ScaffoldGraphEdge* tempEdge = NULL;
	ScaffoldGraphEdge* tempEdge1 = NULL;

	for (int i = 0; i < contigCount; i++) {
		tempEdge = scaffoldGraph[i].outEdge;
		while (tempEdge != NULL) {
			tempEdge->visited = false;
			tempEdge = tempEdge->next;
		}
		tempEdge = scaffoldGraph[i].inEdge;
		while (tempEdge != NULL) {
			tempEdge->visited = false;
			tempEdge = tempEdge->next;
		}
	}

	for (i = 0; i < contigCount; i++) {
		tempEdge = scaffoldGraph[i].outEdge;
		while (tempEdge != NULL) {
			edgeNumber++;
			edgeNumberOut++;
			tempEdge = tempEdge->next;
		}
		tempEdge = scaffoldGraph[i].inEdge;
		while (tempEdge != NULL) {
			edgeNumber++;
			edgeNumberIn++;
			tempEdge = tempEdge->next;
		}
	}

	edgeNumber = (edgeNumber / 2);

	lprec* lp;
	int Ncol, * colno = NULL, ret = 0;
	REAL* row = NULL;

	Ncol = contigCount + edgeNumber;
	lp = make_lp(0, Ncol);
	if (lp == NULL) {
		printf("couldn't construct a new linear programming model");
		exit(0);
	}
	double* weight = new double[edgeNumber];
	long int* edgeLeftNode = new long int[edgeNumber];
	long int* edgeRightNode = new long int[edgeNumber];
	long int* gapDistance = new long int[edgeNumber];


	colno = (int*)malloc(Ncol * sizeof(*colno));
	row = (REAL*)malloc(Ncol * sizeof(*row));
	if ((colno == NULL) || (row == NULL)) {
		printf("couldn't new colno and row");
		exit(0);
	}

	set_add_rowmode(lp, TRUE);

	long int constraintNumber = 0;

	for (i = 0; i < contigCount; i++) {
		for (int q = 0; q < 2; q++) {
			if (q == 0) {
				tempEdge = scaffoldGraph[i].outEdge;
			}
			else {
				tempEdge = scaffoldGraph[i].inEdge;
			}
			long int start = 0;
			while (tempEdge != NULL) {
				if (tempEdge->visited == false) {

					if ((contigOrientation[i] == 1 && q == 0) || (contigOrientation[i] == 0 && q == 1)) {
						j = 0;
						colno[j] = i + 1;
						row[j++] = -1;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = 1;


						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c + contigSetHead->contigSet[i].contigLength + tempEdge->gapDistance)) {
							printf("couldn't add_constraintex");
							exit(0);
						}


						j = 0;
						colno[j] = i + 1;
						row[j++] = 1;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = -1;
						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c - contigSetHead->contigSet[i].contigLength - tempEdge->gapDistance)) {
							printf("couldn't add_constraintex");
							exit(0);
						}

						edgeLeftNode[p - 1] = i;
						edgeRightNode[p - 1] = tempEdge->nodeIndex;

						contigVisited[i] = true;
						contigVisited[tempEdge->nodeIndex] = true;
					}
					if ((contigOrientation[i] == 0 && q == 0) || (contigOrientation[i] == 1 && q == 1)) {

						j = 0;
						colno[j] = i + 1;
						row[j++] = 1;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = -1;


						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c + contigSetHead->contigSet[tempEdge->nodeIndex].contigLength + tempEdge->gapDistance)) {
							printf("couldn't add_constraintex");
							exit(0);
						}



						j = 0;
						colno[j] = i + 1;
						row[j++] = -1;
						colno[j] = tempEdge->nodeIndex + 1;
						row[j++] = 1;

						colno[j] = contigCount + p;
						row[j++] = c;
						if (!add_constraintex(lp, j, row, colno, LE, c - contigSetHead->contigSet[tempEdge->nodeIndex].contigLength - tempEdge->gapDistance)) {
							printf("couldn't add_constraintex");
							exit(0);
						}
						edgeLeftNode[p - 1] = tempEdge->nodeIndex;
						edgeRightNode[p - 1] = i;

						contigVisited[i] = true;
						contigVisited[tempEdge->nodeIndex] = true;
					}

					j = 0;
					colno[j] = contigCount + p;
					row[j++] = 1;
					add_constraintex(lp, j, row, colno, LE, 1);
					j = 0;
					colno[j] = contigCount + p;
					row[j++] = 1;
					add_constraintex(lp, j, row, colno, GE, 0);

					weight[p - 1] = tempEdge->weight;
					gapDistance[p - 1] = tempEdge->gapDistance;
					p++;
					constraintNumber = constraintNumber + 4;
					tempEdge->visited = true;
					tempEdge1 = scaffoldGraph[tempEdge->nodeIndex].outEdge;
					while (tempEdge1 != NULL) {
						if (tempEdge1->nodeIndex == i) {
							tempEdge1->visited = true;
							break;
						}
						tempEdge1 = tempEdge1->next;
					}
					tempEdge1 = scaffoldGraph[tempEdge->nodeIndex].inEdge;
					while (tempEdge1 != NULL) {
						if (tempEdge1->nodeIndex == i) {
							tempEdge1->visited = true;
							break;
						}
						tempEdge1 = tempEdge1->next;
					}

				}
				tempEdge = tempEdge->next;
			}

		}
	}

	for (i = 0; i < contigCount; i++) {
		set_int(lp, i + 1, TRUE);
	}

	p--;



	j = 0;
	for (i = 0; i < p; i++) {
		colno[j] = contigCount + i + 1;
		row[j] = weight[j];
		j++;
	}

	if (!set_obj_fnex(lp, j, row, colno)) {
		printf("couldn't set_obj_fnex");
		exit(0);
	}
	set_add_rowmode(lp, FALSE);


	set_maxim(lp);
	set_timeout(lp, 1200);

	ret = solve(lp);

	if (!(ret == 0 || ret == 1)) {

		set_break_at_first(lp, true);
		ret = solve(lp);
		if (!(ret == 0 || ret == 1)) {
			return NULL;
		}
	}

	REAL* pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
	get_primal_solution(lp, pv);



	for (i = contigCount + constraintNumber + 1; i < contigCount + constraintNumber + p + 1; i++) {
		if (pv[i] < 1) {
			long int d = pv[1 + constraintNumber + edgeRightNode[i - contigCount - constraintNumber - 1]]
				- pv[1 + constraintNumber + edgeLeftNode[i - contigCount - constraintNumber - 1]]
				- contigSetHead->contigSet[edgeLeftNode[i - contigCount - constraintNumber - 1]].contigLength;
			double varD = double(labs(d - gapDistance[i - contigCount - constraintNumber - 1])) / 500;
			if (varD > 3) {
				DeleteSpecialScaffoldEdge(scaffoldGraph, edgeLeftNode[i - contigCount - constraintNumber - 1], edgeRightNode[i - contigCount - constraintNumber - 1]);
				DeleteSpecialScaffoldEdge(scaffoldGraph, edgeRightNode[i - contigCount - constraintNumber - 1], edgeLeftNode[i - contigCount - constraintNumber - 1]);
				continue;
			}
		}

	}


	for (i = 0; i < contigCount; i++) {

		contigPosition[i] = pv[1 + constraintNumber + i];
	}

	long int trueNumber = 0;
	long int realTrueNumber = 0;


	delete[] contigVisited;
	delete[] weight;
	delete[] edgeLeftNode;
	delete[] edgeRightNode;
	delete[] gapDistance;
	delete[] pv;

	delete_lp(lp);

}


void SortNode(int* sortNode, int* length, long int left, long int right) {
	if (left >= right) {
		return;
	}

	long int i = left;
	long int j = right;
	long int key = length[left];
	long int index = sortNode[left];

	while (i < j) {
		while (i < j && (key < length[j])) {
			j--;
		}

		if (i < j) {
			length[i] = length[j];
			sortNode[i] = sortNode[j];
			i++;
		}

		while (i < j && (key > length[i])) {
			i++;
		}

		if (i < j) {
			length[j] = length[i];
			sortNode[j] = sortNode[i];
			j--;
		}
	}

	length[i] = key;
	sortNode[i] = index;

	SortNode(sortNode, length, left, i - 1);
	SortNode(sortNode, length, i + 1, right);
}

ScaffoldGraphEdge* GetOptimizeNodeIndex(ScaffoldGraphHead* scaffoldGraphHead, int nodeIndex, bool orientation, ContigSequence* tempContigSequence, bool right, bool* printIndex) {
	ScaffoldGraphEdge* temp = NULL;
	ScaffoldGraphEdge* tempOriginal = NULL;
	int edgeNum = 0;
	if ((orientation == 1 && right == 1) || (orientation == 0 && right == 0)) {
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge;
		while (temp != NULL) {
			edgeNum++;
			temp = temp->next;
		}
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge;
	}
	else {
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge;
		while (temp != NULL) {
			edgeNum++;
			temp = temp->next;
		}
		temp = scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge;

	}
	return temp;
}


int SearchScaffoldEdge(int index, ContigSequence* contigSequence) {
	int num = 0;
	while (contigSequence != NULL) {
		if (contigSequence->index == index) {
			num++;
		}
		contigSequence = contigSequence->next;
	}

	return num;
}


int GetContigSequenceNum(ContigSequence* tempContigSequence) {
	int contigSequenceNum = 0;
	while (tempContigSequence != NULL) {
		contigSequenceNum++;
		tempContigSequence = tempContigSequence->next;
	}
	return contigSequenceNum;
}


void OptimizeScaffoldSetCongtigSequence(ScaffoldSet* scaffoldSet, int contigNum) {
	ScaffoldSet* tempScaffoldSet = scaffoldSet;
	ContigSequence* tempContigSequence = NULL;

	tempScaffoldSet = scaffoldSet;
	while (tempScaffoldSet != NULL) {
		int contigSequenceNum = GetContigSequenceNum(tempScaffoldSet->contigSequence);
		tempScaffoldSet->contigNum = contigSequenceNum;
		if (contigSequenceNum > 0) {
			tempScaffoldSet->contigIndex = (int*)malloc(sizeof(int) * contigSequenceNum);
		}
		else {
			tempScaffoldSet->contigIndex = NULL;
		}

		tempScaffoldSet = tempScaffoldSet->next;
	}
	tempScaffoldSet = scaffoldSet;
	while (tempScaffoldSet != NULL) {
		tempContigSequence = tempScaffoldSet->contigSequence;
		int p = 0;
		while (tempContigSequence != NULL) {
			tempScaffoldSet->contigIndex[p] = tempContigSequence->index;
			p++;
			tempContigSequence = tempContigSequence->next;
		}
		tempScaffoldSet = tempScaffoldSet->next;
	}

}


bool SearchInsertContigSequence(ScaffoldSetHead* tempInsertSequenceHead, int startIndex, int endIndex, int* contigIndex, int* distance, bool* orientation, int* overlapLength, int count) {
	ScaffoldSet* tempScaffoldSet = tempInsertSequenceHead->scaffoldSet;
	while (tempScaffoldSet != NULL) {
		if (tempScaffoldSet->contigSequence == NULL || tempScaffoldSet->contigNum != labs(endIndex - startIndex) + 1) {
			tempScaffoldSet = tempScaffoldSet->next;
			continue;
		}

		ContigSequence* tempContigSequence = tempScaffoldSet->contigSequence;
		if (startIndex < endIndex) {
			bool token = false;
			for (int i = startIndex; i <= endIndex; i++) {
				if (tempContigSequence->index != contigIndex[i] || tempContigSequence->orientation != orientation[i]) {
					token = true;
					break;
				}
				tempContigSequence = tempContigSequence->next;

			}
			if (token == true) {
				return false;
			}
			tempScaffoldSet->sequenceCount++;
			return true;
		}

		if (endIndex < startIndex) {
			bool token = false;
			for (int i = startIndex; i >= endIndex; i--) {
				if (tempContigSequence->index != contigIndex[i] || tempContigSequence->orientation == orientation[i]) {
					token = true;
					break;
				}
				tempContigSequence = tempContigSequence->next;
			}
			if (token == true) {
				return false;
			}
			tempScaffoldSet->sequenceCount++;
			return true;
		}

		tempScaffoldSet = tempScaffoldSet->next;
	}
	return false;

}

void InsertShortContigToSequence(ContigSequence* tempContigSequence, ContigSequence* preContigSequence, char* file, char* line, int maxSize, bool* printIndex, long int lineCount) {
	FILE* fp;
	if ((fp = fopen(file, "r")) == NULL) {
		printf("%s, does not exist!", file);
		exit(0);
	}
	long int index = -1;

	const char* split = ",";
	char* p;

	int maxCount = 100;
	int contigIndex[maxCount];
	int distance[maxCount];
	bool orientation[maxCount];
	int overlapLength[maxCount];
	ScaffoldSetHead* tempInsertSequenceHead = (ScaffoldSetHead*)malloc(sizeof(ScaffoldSetHead));
	tempInsertSequenceHead->scaffoldSet = NULL;

	while ((fgets(line, maxSize, fp)) != NULL) {
		index++;
		p = strtok(line, split);
		int count = atoi(p);
		if (count < 3) {
			continue;
		}

		if (count > maxCount) {
			continue;
		}
		for (int i = 0; i < maxCount; i++) {
			contigIndex[i] = -1;
			distance[i] = -1;
			orientation[i] = false;
			overlapLength[i] = -1;
		}
		int a = 1;
		while (a <= count) {
			p = strtok(NULL, split);
			contigIndex[a - 1] = atoi(p);
			p = strtok(NULL, split);
			distance[a - 1] = atoi(p);
			p = strtok(NULL, split);
			orientation[a - 1] = atoi(p);
			p = strtok(NULL, split);
			overlapLength[a - 1] = atoi(p);
			a++;
		}

		int endIndex = -1;
		int startIndex = -1;
		for (int i = 0; i < count; i++) {
			if (contigIndex[i] == tempContigSequence->index) {
				endIndex = i;
			}
			if (contigIndex[i] == preContigSequence->index) {
				startIndex = i;
			}

		}
		if (startIndex == -1 || endIndex == -1 || labs(startIndex - endIndex) == 1) {
			continue;
		}

		bool token = SearchInsertContigSequence(tempInsertSequenceHead, startIndex, endIndex, contigIndex, distance, orientation, overlapLength, count);

		if (token == true) {
			continue;
		}

		if (orientation[startIndex] == preContigSequence->orientation && orientation[endIndex] == tempContigSequence->orientation && endIndex > startIndex) {

			ScaffoldSet* tempScaffoldSet = (ScaffoldSet*)malloc(sizeof(ScaffoldSet));
			tempScaffoldSet->contigSequence = NULL;
			tempScaffoldSet->sequenceCount = 1;
			tempScaffoldSet->contigNum = endIndex - startIndex + 1;
			tempScaffoldSet->next = NULL;
			ContigSequence* tempContigSequence0 = tempScaffoldSet->contigSequence;
			for (int i = startIndex; i <= endIndex; i++) {
				ContigSequence* tempContigSequence1 = (ContigSequence*)malloc(sizeof(ContigSequence));
				tempContigSequence1->index = contigIndex[i];
				tempContigSequence1->orientation = orientation[i];
				tempContigSequence1->gapDistance = distance[i];
				tempContigSequence1->next = NULL;
				if (tempContigSequence0 != NULL) {
					tempContigSequence0->next = tempContigSequence1;
				}
				else {
					tempScaffoldSet->contigSequence = tempContigSequence1;
				}
				tempContigSequence0 = tempContigSequence1;
			}
			tempScaffoldSet->next = tempInsertSequenceHead->scaffoldSet;
			tempInsertSequenceHead->scaffoldSet = tempScaffoldSet;
		}

		if (orientation[startIndex] != preContigSequence->orientation && orientation[endIndex] != tempContigSequence->orientation && startIndex > endIndex) {

			ScaffoldSet* tempScaffoldSet = (ScaffoldSet*)malloc(sizeof(ScaffoldSet));
			tempScaffoldSet->contigSequence = NULL;
			tempScaffoldSet->sequenceCount = 1;
			tempScaffoldSet->contigNum = startIndex - endIndex + 1;
			tempScaffoldSet->next = NULL;
			ContigSequence* tempContigSequence0 = tempScaffoldSet->contigSequence;
			for (int i = startIndex; i >= endIndex; i--) {
				ContigSequence* tempContigSequence1 = (ContigSequence*)malloc(sizeof(ContigSequence));
				tempContigSequence1->index = contigIndex[i];
				tempContigSequence1->orientation = !orientation[i];
				if (i == endIndex) {
					tempContigSequence1->gapDistance = 0;
				}
				else {
					tempContigSequence1->gapDistance = distance[i - 1];
				}
				tempContigSequence1->next = NULL;
				if (tempContigSequence0 != NULL) {
					tempContigSequence0->next = tempContigSequence1;
				}
				else {
					tempScaffoldSet->contigSequence = tempContigSequence1;
				}
				tempContigSequence0 = tempContigSequence1;
			}
			tempScaffoldSet->next = tempInsertSequenceHead->scaffoldSet;
			tempInsertSequenceHead->scaffoldSet = tempScaffoldSet;

		}

	}

	ScaffoldSet* tempScaffoldSet = tempInsertSequenceHead->scaffoldSet;
	maxCount = -1;
	while (tempScaffoldSet != NULL) {
		if (tempScaffoldSet->sequenceCount > maxCount) {
			maxCount = tempScaffoldSet->sequenceCount;
		}
		tempScaffoldSet = tempScaffoldSet->next;
	}

	tempScaffoldSet = tempInsertSequenceHead->scaffoldSet;
	ContigSequence* tempContigSequence1 = NULL;
	int largeNum = 0;
	while (tempScaffoldSet != NULL) {
		if (tempScaffoldSet->sequenceCount == maxCount) {
			tempContigSequence1 = tempScaffoldSet->contigSequence;
			largeNum++;
		}
		tempScaffoldSet = tempScaffoldSet->next;
	}

	if ((maxCount > 1 || largeNum == 1) && tempContigSequence1 != NULL) {

		preContigSequence->gapDistance = tempContigSequence1->gapDistance;
		tempContigSequence1 = tempContigSequence1->next;
		ContigSequence* first = preContigSequence;
		while (tempContigSequence1 != NULL && tempContigSequence1->index != tempContigSequence->index) {
			ContigSequence* tempContigSequence2 = (ContigSequence*)malloc(sizeof(ContigSequence));
			printIndex[tempContigSequence1->index] = true;
			tempContigSequence2->index = tempContigSequence1->index;
			tempContigSequence2->orientation = tempContigSequence1->orientation;
			tempContigSequence2->gapDistance = tempContigSequence1->gapDistance;
			tempContigSequence2->next = tempContigSequence;

			first->next = tempContigSequence2;
			first = tempContigSequence2;
			tempContigSequence1 = tempContigSequence1->next;
		}
		first->next = tempContigSequence;
	}

	fclose(fp);

}



void InsertContigToSequence(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, bool* printIndex, char* file, char* line, int maxSize) {
	long int lineCount = 0;
	FILE* fp;
	if ((fp = fopen(file, "r")) == NULL) {
		printf("%s, does not exist!", file);	
		exit(0);
	}
	char* p;
	const char* split = ",";
	while ((fgets(line, maxSize, fp)) != NULL) {
		lineCount++;
	}
	fclose(fp);

	ScaffoldSet* tempScaffoldSet = scaffoldSetHead->scaffoldSet;
	int i = 0;
	while (tempScaffoldSet != NULL) {
		i++;
		if (tempScaffoldSet->contigNum < 2) {
			tempScaffoldSet = tempScaffoldSet->next;
			continue;
		}
		ContigSequence* tempContigSequence = tempScaffoldSet->contigSequence;
		ContigSequence* preContigSequence = NULL;
		if (tempContigSequence == NULL) {
			tempScaffoldSet = tempScaffoldSet->next;
			continue;
		}

		while (tempContigSequence != NULL) {
			if (preContigSequence != NULL) {

				InsertShortContigToSequence(tempContigSequence, preContigSequence, file, line, maxSize, printIndex, lineCount);
			}
			preContigSequence = tempContigSequence;
			tempContigSequence = tempContigSequence->next;
		}

		tempScaffoldSet = tempScaffoldSet->next;
	}



}


ScaffoldSetHead* GetScaffoldSet(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, SimpleResultHead* simpleResultHead, char* line, int maxSize, char* file, int readType) {
	bool* contigOrientation = (bool*)malloc(sizeof(bool) * scaffoldGraphHead->scaffoldGraphNodeCount);
	long int* contigPosition = (long int*)malloc(sizeof(long int) * scaffoldGraphHead->scaffoldGraphNodeCount);
	long int allContigLength = contigSetHead->allContigLength;

	DetermineOrientationOfContigs(scaffoldGraphHead->scaffoldGraph, scaffoldGraphHead->scaffoldGraphNodeCount, contigOrientation);

	IterativeDetermineOrderOfContigs(contigSetHead, scaffoldGraphHead->scaffoldGraph, scaffoldGraphHead->scaffoldGraphNodeCount, contigOrientation, contigPosition, allContigLength);

	ScaffoldGraphEdge* tempEdge = NULL;
	ScaffoldGraphEdge* temp = NULL;
	ScaffoldGraphEdge* t = NULL;
	ScaffoldGraphEdge* tempEdge1 = NULL;
	ScaffoldGraphEdge* temp1 = NULL;
	ScaffoldGraphEdge* t1 = NULL;
	double maxShortPer = -1;
	double maxLongPer = -1;
	int maxSIndex = -1;
	int maxLIndex = -1;
	int maxIndex = -1;
	bool maxOrientation = false;

	double maxWeight = -1;

	for (int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		if (scaffoldGraphHead->scaffoldGraph[i].outEdge != NULL) {
			tempEdge = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			double num = 0;
			double outshortNum = 0;
			double outlongNum = 0;
			while (tempEdge != NULL) {
				num++;
				outshortNum = outshortNum + tempEdge->shortReadCount;
				outlongNum = outlongNum + tempEdge->longReadCount;
				tempEdge = tempEdge->next;
			}

			temp = scaffoldGraphHead->scaffoldGraph[i].outEdge;
			maxShortPer = -1;
			maxLongPer = -1;
			maxSIndex = -1;
			maxLIndex = -1;
			maxIndex = temp->nodeIndex;
			maxOrientation = false;
			maxWeight = -1;

			while (temp != NULL) {
				int nodeIndex = temp->nodeIndex;
				int shortNum = temp->shortReadCount;
				int longNum = temp->longReadCount;
				bool orientation = temp->orientation;
				int tempWeight = temp->weight;

				double shortPer = double(shortNum) / outshortNum;
				double longPer = double(longNum) / outlongNum;
				if (tempWeight <= maxWeight) {
					t = temp->next;
					scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, temp->nodeIndex, orientation);
					if (orientation == true) {
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}
					else {
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}
					temp = t;
					continue;
				}
				else {
					if (maxShortPer != -1 && maxLongPer != -1) {
						scaffoldGraphHead->scaffoldGraph[i].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].outEdge, maxIndex, maxOrientation);
						if (maxOrientation == true) {
							scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge, i, maxOrientation);
						}
						else {
							scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge, i, maxOrientation);
						}

					}
					maxShortPer = shortPer;
					maxLongPer = longPer;
					maxOrientation = temp->orientation;
					maxIndex = temp->nodeIndex;
					maxWeight = temp->weight;
				}
				temp = temp->next;
			}
		}
		if (scaffoldGraphHead->scaffoldGraph[i].inEdge != NULL) {
			tempEdge1 = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			double num = 0;
			double inshortNum = 0;
			double inlongNum = 0;
			while (tempEdge1 != NULL) {
				num++;
				inshortNum = inshortNum + tempEdge1->shortReadCount;
				inlongNum = inlongNum + tempEdge1->longReadCount;
				tempEdge1 = tempEdge1->next;
			}

			temp1 = scaffoldGraphHead->scaffoldGraph[i].inEdge;
			maxShortPer = -1;
			maxLongPer = -1;
			maxSIndex = -1;
			maxLIndex = -1;
			maxIndex = temp1->nodeIndex;
			maxOrientation = false;
			maxWeight = -1;

			while (temp1 != NULL) {
				int nodeIndex = temp1->nodeIndex;
				int shortNum = temp1->shortReadCount;
				int longNum = temp1->longReadCount;
				bool orientation = temp1->orientation;
				int tempWeight = temp1->allReadCount;

				double shortPer = shortNum / inshortNum;
				double longPer = longNum / inlongNum;
				if (tempWeight <= maxWeight) {
					t1 = temp1->next;
					scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, temp1->nodeIndex, orientation);
					if (orientation == true) {
						scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].outEdge, i, orientation);
					}
					else {
						scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[nodeIndex].inEdge, i, orientation);
					}
					temp1 = t1;
					continue;
				}
				else {
					if (maxShortPer != -1 && maxLongPer != -1) {
						scaffoldGraphHead->scaffoldGraph[i].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[i].inEdge, maxIndex, maxOrientation);
						if (maxOrientation == true) {
							scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[maxIndex].outEdge, i, maxOrientation);
						}
						else {
							scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge = DeleteEdgeNode(scaffoldGraphHead->scaffoldGraph[maxIndex].inEdge, i, maxOrientation);
						}

					}
					maxShortPer = shortPer;
					maxLongPer = longPer;
					maxOrientation = temp1->orientation;
					maxIndex = temp1->nodeIndex;
					maxWeight = temp1->weight;
				}
				temp1 = temp1->next;
			}
		}
	}


	int i = 0;
	int j = 0;
	bool* printIndex = (bool*)malloc(sizeof(bool) * scaffoldGraphHead->scaffoldGraphNodeCount);
	int* sortNode = (int*)malloc(sizeof(int) * scaffoldGraphHead->scaffoldGraphNodeCount);
	int* tempLength = (int*)malloc(sizeof(int) * scaffoldGraphHead->scaffoldGraphNodeCount);
	for (i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		printIndex[i] = false;
		sortNode[i] = i;
		tempLength[i] = contigSetHead->contigSet[i].contigLength;
	}
	SortNode(sortNode, tempLength, 0, scaffoldGraphHead->scaffoldGraphNodeCount - 1);
	ScaffoldGraph* scaffoldGraph = (ScaffoldGraph*)malloc(sizeof(ScaffoldGraph));
	ScaffoldSetHead* scaffoldSetHead = (ScaffoldSetHead*)malloc(sizeof(ScaffoldSetHead));
	scaffoldSetHead->scaffoldSet = NULL;
	bool orientation = true;
	int a = 0;
	if (readType == 1) {
		a = 2000;
	}
	else if (readType == 2) {
		a = 1000;
	}

	for (long int p = 0; p < scaffoldGraphHead->scaffoldGraphNodeCount; p++) {
		i = sortNode[p];
		if (printIndex[i] == true || contigSetHead->contigSet[i].contigLength < a || (scaffoldGraphHead->scaffoldGraph[i].outEdge == NULL && scaffoldGraphHead->scaffoldGraph[i].inEdge == NULL)) {
			continue;
		}
		ScaffoldGraphEdge* temp = scaffoldGraphHead->scaffoldGraph[i].outEdge;

		ScaffoldSet* tempScaffoldSet = (ScaffoldSet*)malloc(sizeof(ScaffoldSet));
		tempScaffoldSet->next = NULL;
		if (scaffoldSetHead->scaffoldSet != NULL) {
			tempScaffoldSet->next = scaffoldSetHead->scaffoldSet;
		}
		scaffoldSetHead->scaffoldSet = tempScaffoldSet;

		ContigSequence* tempContigSequence = (ContigSequence*)malloc(sizeof(ContigSequence));
		tempContigSequence->index = i;
		tempContigSequence->orientation = 1;
		tempContigSequence->gapDistance = 0;
		tempContigSequence->next = NULL;
		scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequence;
		j = i;
		if (temp == NULL) {
			printIndex[i] = true;
		}
		orientation = 1;
		while (temp != NULL) {
			printIndex[i] = true;
			temp = GetOptimizeNodeIndex(scaffoldGraphHead, j, orientation, scaffoldSetHead->scaffoldSet->contigSequence, 1, printIndex);
			if (temp == NULL) {
				break;
			}

			j = temp->nodeIndex;
			if ((orientation == 1 && temp->orientation == 1) || (orientation == 0 && temp->orientation == 0)) {
				orientation = 1;
			}
			else {
				orientation = 0;
			}

			int cc = SearchScaffoldEdge(temp->nodeIndex, scaffoldSetHead->scaffoldSet->contigSequence);
			if (cc > 0) {
				printIndex[j] = true;
				break;
			}
			printIndex[j] = true;

			ContigSequence* tempContigSequence1 = (ContigSequence*)malloc(sizeof(ContigSequence));
			tempContigSequence1->index = temp->nodeIndex;
			tempContigSequence1->gapDistance = 0;
			tempContigSequence->gapDistance = temp->gapDistance;

			tempContigSequence->next = tempContigSequence1;
			tempContigSequence1->next = NULL;
			tempContigSequence = tempContigSequence1;
			tempContigSequence1 = NULL;
			tempContigSequence->orientation = orientation;

		}
		temp = scaffoldGraph[i].inEdge;
		if (temp == NULL) {
			continue;
		}
		orientation = 1;
		j = i;
		printIndex[i] = false;
		while (temp != NULL) {
			printIndex[j] = true;
			temp = GetOptimizeNodeIndex(scaffoldGraphHead, j, orientation, scaffoldSetHead->scaffoldSet->contigSequence, 0, printIndex);
			if (temp == NULL) {
				break;
			}
			j = temp->nodeIndex;

			if ((orientation == 1 && temp->orientation == 1) || (orientation == 0 && temp->orientation == 0)) {
				orientation = 1;
			}
			else {
				orientation = 0;
			}

			int cc = SearchScaffoldEdge(temp->nodeIndex, scaffoldSetHead->scaffoldSet->contigSequence);
			if (cc > 0) {
				printIndex[j] = true;
				break;
			}
			ContigSequence* tempContigSequence1 = (ContigSequence*)malloc(sizeof(ContigSequence));
			tempContigSequence1->index = temp->nodeIndex;
			tempContigSequence1->gapDistance = temp->gapDistance;
			tempContigSequence1->next = scaffoldSetHead->scaffoldSet->contigSequence;
			scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequence1;
			scaffoldSetHead->scaffoldSet->contigSequence->orientation = orientation;
		}
	}

	OptimizeScaffoldSetCongtigSequence(scaffoldSetHead->scaffoldSet, contigSetHead->contigCount);


	if (readType == 1) {
		InsertContigToSequence(contigSetHead, scaffoldSetHead, printIndex, file, line, maxSize);
	}
	
	for (i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		if (printIndex[i] != 1) {
			ScaffoldSet* tempScaffoldSet = (ScaffoldSet*)malloc(sizeof(ScaffoldSet));
			tempScaffoldSet->next = scaffoldSetHead->scaffoldSet;
			scaffoldSetHead->scaffoldSet = tempScaffoldSet;

			ContigSequence* tempContigSequence = (ContigSequence*)malloc(sizeof(ContigSequence));
			tempContigSequence->index = i;
			tempContigSequence->gapDistance = 0;
			tempContigSequence->orientation = 1;
			tempContigSequence->next = NULL;
			scaffoldSetHead->scaffoldSet->contigSequence = tempContigSequence;
		}
	}

	return scaffoldSetHead;
}

char* ReverseComplement(char* temp) {
	int len = strlen(temp);
	char* rcTemp = (char*)malloc(sizeof(char) * (len + 1));
	for (int i = 0; i < len; i++) {
		if (temp[i] == 'A' || temp[i] == 'a') {
			rcTemp[len - 1 - i] = 'T';
		}
		else if (temp[i] == 'T' || temp[i] == 't') {
			rcTemp[len - 1 - i] = 'A';
		}
		else if (temp[i] == 'G' || temp[i] == 'g') {
			rcTemp[len - 1 - i] = 'C';
		}
		else if (temp[i] == 'C' || temp[i] == 'c') {
			rcTemp[len - 1 - i] = 'G';
		}
		else if (temp[i] == 'N' || temp[i] == 'n') {
			rcTemp[len - 1 - i] = 'N';
		}
	}
	rcTemp[len] = '\0';
	return rcTemp;
}



void OutputScaffoldSet(ScaffoldSet* scaffoldSet, ContigSetHead* contigSetHead, char* result) {
	//cout << "begin cout" << endl;
	int i = 0;
	int j = 0;
	Contig* contigSet = contigSetHead->contigSet;
	int contigCount = contigSetHead->contigCount;

	fflush(stdout);
	setvbuf(stdout, NULL, _IONBF, 0);

	bool* printContigIndex = new bool[contigCount];
	for (i = 0; i < contigCount; i++) {
		printContigIndex[i] = false;
	}
	ScaffoldSet* tempScaffoldSet = scaffoldSet;
	scaffoldSet = tempScaffoldSet;

	char* scaffoldTagFileName = new char[1000];
	strcpy(scaffoldTagFileName, result);
	strcat(scaffoldTagFileName, "_tag.fa");
	ofstream ocoutTag;
	ocoutTag.open(scaffoldTagFileName);

	char* scaffoldSetFileName = new char[1000];
	strcpy(scaffoldSetFileName, result);
	strcat(scaffoldSetFileName, "_set.fa");
	ofstream ocout1;
	ocout1.open(scaffoldSetFileName);
	j = 0;

	int tempLength = 0;
	int tempGapDis = 0;

	while (tempScaffoldSet != NULL) {
		ContigSequence* tempContigSequence = tempScaffoldSet->contigSequence;
		if (tempContigSequence == NULL) {
			ocoutTag << endl;
			tempScaffoldSet = tempScaffoldSet->next;
			continue;
		}
		ocout1 << ">" << j << endl;

		int allLength = 0;
		tempLength = 0;
		tempGapDis = 0;

		while (tempContigSequence != NULL) {

			if (printContigIndex[tempContigSequence->index] == true) {
				ocoutTag << "--";
			}
			tempLength = tempGapDis + tempLength + strlen(contigSet[tempContigSequence->index].contig);
			printContigIndex[tempContigSequence->index] = true;
			ocoutTag << tempContigSequence->index << "(" << tempContigSequence->gapDistance << "--" << tempContigSequence->orientation << "--" << tempLength << "),";
			tempGapDis = tempContigSequence->gapDistance;

			if (tempContigSequence->orientation == 0) {
				char* tempChar1 = ReverseComplement(contigSet[tempContigSequence->index].contig);
				if (tempContigSequence->gapDistance < 0 && tempContigSequence->next != NULL) {
					int contigLength = strlen(tempChar1);
					if (contigLength + tempContigSequence->gapDistance < 0) {
						tempContigSequence = tempContigSequence->next;
						continue;
					}

					char* tempChar2 = (char*)malloc(sizeof(char) * (contigLength + tempContigSequence->gapDistance + 1));

					strncpy(tempChar2, tempChar1, contigLength + tempContigSequence->gapDistance);

					tempChar2[contigLength + tempContigSequence->gapDistance] = '\0';

					free(tempChar1);

					tempChar1 = tempChar2;
				}
				ocout1 << tempChar1;
				free(tempChar1);
			}
			else {
				if (tempContigSequence->gapDistance < 0 && tempContigSequence->next != NULL) {

					int contigLength = strlen(contigSet[tempContigSequence->index].contig);

					if (contigLength + tempContigSequence->gapDistance < 0) {
						tempContigSequence = tempContigSequence->next;
						continue;
					}

					char* tempChar2 = (char*)malloc(sizeof(char) * (contigLength + tempContigSequence->gapDistance + 1));
					strncpy(tempChar2, contigSet[tempContigSequence->index].contig, contigLength + tempContigSequence->gapDistance);

					tempChar2[contigLength + tempContigSequence->gapDistance] = '\0';

					ocout1 << tempChar2;
					free(tempChar2);
				}
				else {
					ocout1 << contigSet[tempContigSequence->index].contig;
				}
			}
			if (tempContigSequence->gapDistance >= 0 && tempContigSequence->next != NULL) {
				int cc = tempContigSequence->gapDistance;
				for (int tt = 0; tt < cc; tt++) {
					ocout1 << "N";
				}
			}
			tempContigSequence = tempContigSequence->next;
		}
		ocoutTag << endl;
		ocout1 << endl;
		j++;
		tempScaffoldSet = tempScaffoldSet->next;

	}

	ocoutTag << "----------------------------------------------------------------" << endl;
	for (i = 0; i < contigCount; i++) {
		if (printContigIndex[i] == false && contigSet[i].contig != NULL) {
			ocout1 << ">" << j << endl;
			ocout1 << contigSet[i].contig << endl;
			ocoutTag << i << "," << endl;
			j++;
		}
	}

	Contig* minContigSet = contigSetHead->minContigSet;
	int minContigCount = contigSetHead->minContigCount;

	for (i = 0; i < minContigCount; i++) {
		ocout1 << ">min-" << j << endl;
		ocout1 << minContigSet[i].contig << endl;
		ocoutTag << i << "," << endl;
		j++;
	}

}




#endif