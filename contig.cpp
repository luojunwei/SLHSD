#ifndef CONTIG_CPP_INCLUDED 
#define CONTIG_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "contig.h"

using namespace std;

ContigSetHead* GetContigSet(char* contigSetFile, int contigLengthThreshold, int insertSize) {

	ContigSetHead* contigSetHead = new ContigSetHead;
	contigSetHead->contigSet = NULL;
	contigSetHead->contigCount = 0;
	contigSetHead->allContigLength = 0;
	contigSetHead->minContigSet = 0;
	contigSetHead->minContigCount = 0;
	contigSetHead->minAllContigLength = 0;

	long int maxSize = 90000;

	char* contig = NULL;
	if (NULL == (contig = (char*)malloc(sizeof(char) * maxSize))) {
		perror("malloc error!");
		exit(1);
	}

	FILE* fp;
	if ((fp = fopen(contigSetFile, "r")) == NULL) {
		printf("%s, does not exist!", contigSetFile);
		exit(0);
	}

	while ((fgets(contig, maxSize, fp)) != NULL) {
		if (contig[0] == '>') {
			contigSetHead->contigCount++;
		}
	}
	fclose(fp);

	long int contigCount = 0;
	contigCount = contigSetHead->contigCount;

	contigSetHead->contigSet = new Contig[contigCount];

	for (long int i = 0; i < contigCount; i++) {
		contigSetHead->contigSet[i].contig = NULL;
		contigSetHead->contigSet[i].contigLength = 0;
		contigSetHead->contigSet[i].name = NULL;
		contigSetHead->contigSet[i].contigId = -1;
		contigSetHead->contigSet[i].repeativeContig = false;
		contigSetHead->contigSet[i].shortContig = false;
	}

	if ((fp = fopen(contigSetFile, "r")) == NULL) {
		printf("%s, does not exist!", contigSetFile);
		exit(0);
	}

	long int allocateLength = 0;
	long int contigIndex = -1;

	while ((fgets(contig, maxSize, fp)) != NULL) {
		if (contig[0] == '>') {
			if (strlen(contig) == maxSize - 1) {
				while ((fgets(contig, maxSize, fp)) != NULL) {
					if (strlen(contig) != maxSize - 1) {
						break;
					}
				}
			}
			if (contigIndex != -1) {
				contigSetHead->allContigLength = contigSetHead->allContigLength + contigSetHead->contigSet[contigIndex].contigLength;
			}

			contigIndex++;
			int len = strlen(contig);
			contigSetHead->contigSet[contigIndex].name = (char*)malloc(sizeof(char) * len);
			strncpy(contigSetHead->contigSet[contigIndex].name, contig + 1, len - 1);
			contigSetHead->contigSet[contigIndex].name[len - 2] = '\0';
			contigSetHead->contigSet[contigIndex].contigId = contigIndex;

			continue;
		}

		long int extendLength = strlen(contig);
		if (contig[extendLength - 1] == '\n' || contig[extendLength - 1] == '\t' || contig[extendLength - 1] == '\r') {
			extendLength--;
		}

		long int contigLength = 0;
		char* tempContig = NULL;
		if (contigSetHead->contigSet[contigIndex].contig != NULL) {
			if (contigSetHead->contigSet[contigIndex].contigLength + extendLength >= allocateLength) {
				contigLength = contigSetHead->contigSet[contigIndex].contigLength;
				contigSetHead->contigSet[contigIndex].contig = (char*)realloc(contigSetHead->contigSet[contigIndex].contig, allocateLength + maxSize + 1);

				allocateLength = allocateLength + maxSize + 1;

				strncpy(contigSetHead->contigSet[contigIndex].contig + contigLength, contig, extendLength); 
				contigSetHead->contigSet[contigIndex].contig[contigLength + extendLength] = '\0';
				contigSetHead->contigSet[contigIndex].contigLength = contigLength + extendLength;

			}
			else {
				strncpy(contigSetHead->contigSet[contigIndex].contig + contigSetHead->contigSet[contigIndex].contigLength, contig, extendLength);
				contigSetHead->contigSet[contigIndex].contig[contigSetHead->contigSet[contigIndex].contigLength + extendLength] = '\0';
				contigSetHead->contigSet[contigIndex].contigLength = contigSetHead->contigSet[contigIndex].contigLength + extendLength;
			}

		}
		else {
			contigSetHead->contigSet[contigIndex].contig = (char*)malloc(sizeof(char) * (maxSize + 1));
			strncpy(contigSetHead->contigSet[contigIndex].contig, contig, extendLength);
			contigSetHead->contigSet[contigIndex].contig[extendLength] = '\0';
			contigSetHead->contigSet[contigIndex].contigLength = extendLength;
			allocateLength = maxSize + 1;
		}

	}
	contigSetHead->allContigLength = contigSetHead->allContigLength + contigSetHead->contigSet[contigIndex].contigLength;
	fflush(fp);
	fclose(fp);

	if (insertSize < 1000) {
		contigLengthThreshold = 1500;
	}
	else {
		contigLengthThreshold = 1000;
	}
	int num = 0;
	for (long int i = 0; i < contigSetHead->contigCount; i++) {
		if (contigSetHead->contigSet[i].contigLength < contigLengthThreshold) {
			contigSetHead->contigSet[i].shortContig = true;
			num++;
		}
		else {
			contigSetHead->contigSet[i].shortContig = false;
		}
	}

	if (contigSetHead->contigCount <= 0) {
		cout << "contig count is zero,please check the contig file!" << endl;
		exit(0);
	}
	outputContigSet(contigSetHead, contigLengthThreshold);
	
	return contigSetHead;
}


void outputContigSet(ContigSetHead* contigSetHead, int contigLengthThreshold) {		
	char* file = (char*)malloc(sizeof(char) * 50);
	strcpy(file, "./SLHSD-OUTPut-Directory/contigSetforAll.fa");
	FILE* fp;
	if ((fp = fopen(file, "w")) == NULL) {
		cout << file << "does not exist!" << endl;
		exit(0);
	}

	for (long int i = 0; i < contigSetHead->contigCount; i++) {
		fprintf(fp, ">%ld--%ld\n", i, contigSetHead->contigSet[i].contigLength);
	}

	fclose(fp);
	free(file);
	file = NULL;
}



/*
void DetectRepeativeContigInSet(ContigSetHead* contigSetHead, char* bamFileName, float ratio) {

	long int contigIndex = 0;
	long int referenceIndex = -1;  //Ϊʲô���index��-1��ʼ
	string readName;
	string previousReadName = "a";
	//cout << previousReadName;

	BamReader bamReader;
	bamReader.Open(bamFileName);//bool Open(const std::string & filename);��BAM�ļ�
	BamAlignment alignment;

	long int contigCount = bamReader.GetReferenceCount();//int GetReferenceCount() const; ���زο����еĸ���
	cout << "contigCount=" <<contigCount << "---" << endl;  //3179��

	while (bamReader.GetNextAlignment(alignment)) {  //bool GetNextAlignment(BamAlignment& alignment);������һ�����õıȶ�
		readName = alignment.Name;  //���Name��BamAlignment�ļ��ж�����---read name
		//cout << "alignment=" << alignment << endl;
		//cout << "readName: " << readName << ",privoiusReadName:" << previousReadName << endl;  //readName��utg0-utg6250 ÿһ��readName�кü���read
		if (previousReadName != "a" && readName != previousReadName) {
			contigIndex++;  //contigIndex��0-3178
		}

		referenceIndex = alignment.RefID;  //�ο����е�ID
		//cout << "referenceIndex=" << referenceIndex << ",contigIndex=" << contigIndex << ",readName=" << readName << ",previousReadName=" << previousReadName << endl;
		if (referenceIndex == contigIndex) {	//����ο����е�id=contigIndex
			previousReadName = readName;
			continue;
		}
		//out << "111" << endl;
		const vector<CigarOp>& cigarData = alignment.CigarData;
		vector<CigarOp>::const_iterator cigarBegin = cigarData.begin();
		vector<CigarOp>::const_iterator cigarIter = cigarBegin;
		//cout << "222" << endl;
		double matchCount = 0;
		double unMatchCount = 0;
		for (; cigarIter != cigarEnd; ++cigarIter) {
			const CigarOp& op = (*cigarIter);
			if (op.Type == 'M') {
				matchCount = matchCount + op.Length;
			}
			else {
				unMatchCount = unMatchCount + op.Length;
			}
		}
		//cout << "333" << endl;
		uint16_t value;
		if (alignment.HasTag("NM")) {   //HasTag: returns true if alignment has a record for this tag name
			if (alignment.GetTag("NM", value)) {  //GetTag: retrieves tag data
				if (value != 0) {
					previousReadName = readName;
					continue;
				}
			}
		}
		//cout << "444" << endl;
		if (matchCount / (matchCount + unMatchCount) >= ratio && matchCount + unMatchCount > 0) {
			contigSetHead->contigSet[contigIndex].repeativeContig = true;
		}
		previousReadName = readName;

	}
}*/

#endif