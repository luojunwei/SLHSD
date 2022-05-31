#ifndef aligningFromBam_CPP_INCLUDED 
#define aligningFromBam_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include<vector>
#include <algorithm>  
#include <functional> 

#include "aligningFromBam.h"

using namespace std;


ReadMapPosition* GetReadInformation(ContigSetHead* contigSetHead, char* leftShortBamFile, char* rightShortBamFile, char* longBamFile, char* originalShortResult, char* originalLongResult, char* originalShortContig, int insertSize, int readType) {
	long int i = 0;
	long int j = 0;
	long int refLength = 0;
	long int refCount = 0;
	
	refCount = contigSetHead->contigCount;
	ReadMapPosition* readMapPosition = new ReadMapPosition[refCount];

	for (i = 0; i < contigSetHead->contigCount; i++) {
		refLength = contigSetHead->contigSet[i].contigLength;
		readMapPosition[i].leftShortIndex = new long int[refLength];
		readMapPosition[i].rightShortIndex = new long int[refLength];
		readMapPosition[i].leftLongIndex = new long int[refLength];
		readMapPosition[i].rightLongIndex = new long int[refLength];


		for (long int j = 0; j < refLength; j++) {
			readMapPosition[i].leftShortIndex[j] = 0;
			readMapPosition[i].rightShortIndex[j] = 0;
			readMapPosition[i].leftLongIndex[j] = 0;
			readMapPosition[i].rightLongIndex[j] = 0;
		}

		readMapPosition[i].leftShortCoverage = 0;
		readMapPosition[i].rightShortCoverage = 0;
		readMapPosition[i].leftLongCoverage = 0;
		readMapPosition[i].rightLongCoverage = 0;
	}
	BamReader bamReaderLeft;
	BamReader bamReaderRight;
	bamReaderLeft.Open(leftShortBamFile);
	bamReaderRight.Open(rightShortBamFile);
	BamAlignment alignmentLeft;
	BamAlignment alignmentRight;

	FILE* fp1;
	long int shortAlignCount = 0;
	if ((fp1 = fopen(originalShortResult, "w")) == NULL) {
		printf("%s, does not exist!", originalShortResult);
		exit(0);
	}

	int contigLengthIgnore = 2000;
	int minScore = contigSetHead->minAlignmentScore;
	bool token = false;

	long int leftid = -1;
	long int rightid = -1;

	if (readType == 1) {
		refLength = 0;
		while (bamReaderLeft.GetNextAlignment(alignmentLeft) && bamReaderRight.GetNextAlignment(alignmentRight)) {
			while ((alignmentLeft.AlignmentFlag & 0x900) != 0) {
				token = false;
				bamReaderLeft.GetNextAlignment(alignmentLeft);
				continue;
			}
			while ((alignmentRight.AlignmentFlag & 0x900) != 0) {
				token = false;
				bamReaderRight.GetNextAlignment(alignmentRight);
				continue;
			}

			refLength = contigSetHead->contigSet[alignmentLeft.RefID].contigLength;
			if (alignmentLeft.IsMapped()) {
				if (alignmentLeft.Position >= 0) {
					if (alignmentLeft.IsReverseStrand() == true) {
						readMapPosition[alignmentLeft.RefID].leftShortIndex[alignmentLeft.Position]++;
						readMapPosition[alignmentLeft.RefID].leftShortCoverage++;
					}
					if (alignmentLeft.IsReverseStrand() == false) {
						readMapPosition[alignmentLeft.RefID].rightShortIndex[alignmentLeft.Position]++;
						readMapPosition[alignmentLeft.RefID].rightShortCoverage++;
					}

				}
			}
			refLength = contigSetHead->contigSet[alignmentRight.RefID].contigLength;
			if (alignmentRight.IsMapped()) {
				if (alignmentRight.Position >= 0) {
					if (alignmentRight.IsReverseStrand() == true) {
						readMapPosition[alignmentRight.RefID].leftShortIndex[alignmentRight.Position]++;
						readMapPosition[alignmentRight.RefID].leftShortCoverage++;
					}
					if (alignmentRight.IsReverseStrand() == false) {
						readMapPosition[alignmentRight.RefID].rightShortIndex[alignmentRight.Position]++;
						readMapPosition[alignmentRight.RefID].rightShortCoverage++;
					}
				}
			}
			if (contigSetHead->contigSet[alignmentLeft.RefID].contigLength <= contigLengthIgnore || contigSetHead->contigSet[alignmentRight.RefID].contigLength <= contigLengthIgnore) {
				continue;
			}

			if (alignmentLeft.IsMapped() && alignmentRight.IsMapped() && alignmentLeft.RefID != alignmentRight.RefID && alignmentLeft.MapQuality > minScore && alignmentRight.MapQuality > minScore) {
				fprintf(fp1, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", alignmentLeft.RefID, alignmentRight.RefID, alignmentLeft.Position, alignmentRight.Position, !alignmentLeft.IsReverseStrand(), !alignmentRight.IsReverseStrand(), alignmentLeft.Length, alignmentRight.Length, alignmentLeft.MapQuality, alignmentRight.MapQuality);
				shortAlignCount++;
			}
		}
	}
	if (readType == 2) {

		refLength = 0;
		while (bamReaderLeft.GetNextAlignment(alignmentLeft) && bamReaderRight.GetNextAlignment(alignmentRight)) {
			while ((alignmentLeft.AlignmentFlag & 0x900) != 0) {
				token = false;
				bamReaderLeft.GetNextAlignment(alignmentLeft);
				continue;
			}
			while ((alignmentRight.AlignmentFlag & 0x900) != 0) {
				token = false;
				bamReaderRight.GetNextAlignment(alignmentRight);
				continue;
			}

			refLength = contigSetHead->contigSet[alignmentLeft.RefID].contigLength;

			if (alignmentLeft.IsMapped()) {
				if (alignmentLeft.Position >= 0) {
					if (refLength >= insertSize * 1.5) {
						if (alignmentLeft.Position <= insertSize * 1.5 && alignmentLeft.IsReverseStrand() == false) {
							readMapPosition[alignmentLeft.RefID].leftShortIndex[alignmentLeft.Position]++;
							readMapPosition[alignmentLeft.RefID].leftShortCoverage++;
						}
						if (alignmentLeft.Position >= refLength - insertSize * 1.5 && alignmentLeft.IsReverseStrand() == true) {
							readMapPosition[alignmentLeft.RefID].rightShortIndex[alignmentLeft.Position]++;
							readMapPosition[alignmentLeft.RefID].rightShortCoverage++;
						}
					}
					else if (refLength < insertSize * 1.5) {
						if (alignmentLeft.Position <= refLength * 0.7 && alignmentLeft.IsReverseStrand() == false) {
							readMapPosition[alignmentLeft.RefID].leftShortIndex[alignmentLeft.Position]++;
							readMapPosition[alignmentLeft.RefID].leftShortCoverage++;
						}
						if (alignmentLeft.Position > refLength * 0.7 && alignmentLeft.IsReverseStrand() == true) {
							readMapPosition[alignmentLeft.RefID].rightShortIndex[alignmentLeft.Position]++;
							readMapPosition[alignmentLeft.RefID].rightShortCoverage++;
						}
					}
				}
			}

			refLength = contigSetHead->contigSet[alignmentRight.RefID].contigLength;

			if (alignmentRight.IsMapped()) {
				if (alignmentRight.Position >= 0) {
					if (refLength >= insertSize * 1.5) {
						if (alignmentRight.Position <= insertSize * 1.5 && alignmentRight.IsReverseStrand() == false) {
							readMapPosition[alignmentRight.RefID].leftShortIndex[alignmentRight.Position]++;
							readMapPosition[alignmentRight.RefID].leftShortCoverage++;
						}
						if (alignmentRight.Position > refLength - insertSize * 1.5 && alignmentRight.IsReverseStrand() == true) {
							readMapPosition[alignmentRight.RefID].rightShortIndex[alignmentRight.Position]++;
							readMapPosition[alignmentRight.RefID].rightShortCoverage++;
						}
					}
					else if (refLength < insertSize * 1.5) {
						if (alignmentRight.Position <= refLength * 0.7 && alignmentRight.IsReverseStrand() == false) {
							readMapPosition[alignmentRight.RefID].leftShortIndex[alignmentRight.Position]++;
							readMapPosition[alignmentRight.RefID].leftShortCoverage++;
						}
						if (alignmentRight.Position > refLength * 0.7 && alignmentRight.IsReverseStrand() == true) {
							readMapPosition[alignmentRight.RefID].rightShortIndex[alignmentRight.Position]++;
							readMapPosition[alignmentRight.RefID].rightShortCoverage++;
						}
					}
				}
			}
			if (contigSetHead->contigSet[alignmentLeft.RefID].contigLength <= contigLengthIgnore || contigSetHead->contigSet[alignmentRight.RefID].contigLength <= contigLengthIgnore) {
				continue;
			}

			if (alignmentLeft.IsMapped() && alignmentRight.IsMapped() && alignmentLeft.RefID != alignmentRight.RefID && alignmentLeft.MapQuality > minScore && alignmentRight.MapQuality > minScore) {
				fprintf(fp1, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", alignmentLeft.RefID, alignmentRight.RefID, alignmentLeft.Position, alignmentRight.Position, !alignmentLeft.IsReverseStrand(), !alignmentRight.IsReverseStrand(), alignmentLeft.Length, alignmentRight.Length, alignmentLeft.MapQuality, alignmentRight.MapQuality);
				shortAlignCount++;
			}
		}
	}


	fflush(fp1);
	fclose(fp1);

	long int readLength = 0;
	long int previousReadLength = 0;
	long int readIndex = 0;

	FILE* fp3;
	if ((fp3 = fopen(originalShortContig, "w")) == NULL) {
		cout << originalShortContig << "does not exist!" << endl;
		exit(0);
	}

	FILE* fp2;
	if ((fp2 = fopen(originalLongResult, "w")) == NULL) {
		printf("%s, does not exist!", originalLongResult);
		exit(0);
	}

	string bamFileName = longBamFile;
	long int contigIndex = -1;
	long int previousContigIndex = -1;
	int readNameLength = 2000;
	string readName;
	string previousReadName = "a";

	BamReader bamReaderLong;
	bamReaderLong.Open(longBamFile);
	BamAlignment alignmentLong;
	long int contigCount = bamReaderLong.GetReferenceCount();
	long int minContigLengthAlign = 0;
	long int readLengthCutOff = 1000;

	LongResultHead* longresult = (LongResultHead*)malloc(sizeof(LongResultHead));
	longresult->allocateLongResultCount = 1000;
	longresult->longResultCount = 0;
	longresult->aligningShortContigResultCount = 0;
	longresult->longResult = (LongResult*)malloc(sizeof(LongResult) * longresult->allocateLongResultCount);
	for (i = 0; i < longresult->allocateLongResultCount; i++) {	
		longresult->longResult[i].readStartPosition = -1;
		longresult->longResult[i].readEndPosition = -1;
		longresult->longResult[i].contigStartPosition = -1;
		longresult->longResult[i].contigEndPosition = -1;
		longresult->longResult[i].overlapLength = -1;
		longresult->longResult[i].contigIndex = -1;
		longresult->longResult[i].overlapLength = -1;
		longresult->longResult[i].orientation = false;
		longresult->longResult[i].quality = 0;
	}


	i = 0;
	if (insertSize <= 1000) {
		contigSetHead->minAlignmentScore = 20;
	}

	while (bamReaderLong.GetNextAlignment(alignmentLong)) {
		if (alignmentLong.IsMapped()) {
			if (alignmentLong.Position >= 0) {

				if (alignmentLong.IsReverseStrand() == 1) {
					readMapPosition[alignmentLong.RefID].leftLongIndex[alignmentLong.Position]++;
					readMapPosition[alignmentLong.RefID].leftLongCoverage++;
				}
				if (alignmentLong.IsReverseStrand() == 0) {
					readMapPosition[alignmentLong.RefID].rightLongIndex[alignmentLong.Position]++;
					readMapPosition[alignmentLong.RefID].rightLongCoverage++;
				}


			}
		}

		readName = alignmentLong.Name;
		if (!alignmentLong.IsMapped() || alignmentLong.MapQuality < contigSetHead->minAlignmentScore || contigSetHead->contigSet[alignmentLong.RefID].contigLength < contigSetHead->contigLengthIgnore) {//
			previousReadName = readName;
			continue;
		}
		contigIndex = alignmentLong.RefID;
		readLength = alignmentLong.Length;

		if (contigSetHead->contigSet[contigIndex].contigLength < minContigLengthAlign || readLength < readLengthCutOff) {
			previousReadName = readName;
			continue;
		}
		if (previousReadName != "a" && readName != previousReadName) {
			if (longresult->longResultCount + longresult->aligningShortContigResultCount > 1) {
				outPutLongAlignResult(longresult, contigSetHead, fp2, fp3, readIndex, previousReadLength);
				readIndex++;
			}
			i = 0;
			longresult->longResultCount = 0;
			longresult->aligningShortContigResultCount = 0;
		}

		if (GetAlignLongResultOneLine(longresult, alignmentLong, contigSetHead, i) != false) {
			i++;
		}
		previousReadLength = readLength;
		previousReadName = readName;
	}
	fflush(fp2);
	fflush(fp3);
	fclose(fp2);
	fclose(fp3);

	return readMapPosition;
}

bool GetAlignLongResultOneLine(LongResultHead* longResultHead, BamAlignment alignment, ContigSetHead* contigSetHead, long int index) {
	int minAlignmentLength = contigSetHead->minAlignmentRevised;//150

	vector<int>clipSizes;
	vector<int>readPositions;
	vector<int>genomePositions;

	int readStartPosition = -1;
	int readEndPosition = -1;
	int contigStartPosition = -1;
	int contigEndPosition = -1;

	long int contigIndex = alignment.RefID;
	long int contigLength = contigSetHead->contigSet[contigIndex].contigLength;
	long int readLength = alignment.Length;

	if (alignment.GetSoftClips(clipSizes, readPositions, genomePositions) != true) {  
		readStartPosition = 0;
		readEndPosition = alignment.Length - 1;
		contigStartPosition = alignment.Position;  
		contigEndPosition = alignment.GetEndPosition() - 1;   
	}
	else {
		if (clipSizes.size() == 1) {
			if (clipSizes[0] == readPositions[0]) { 
				readStartPosition = readPositions[0];
				readEndPosition = alignment.Length - 1;
				contigStartPosition = genomePositions[0];
				contigEndPosition = alignment.GetEndPosition();
			}
			else {	
				readStartPosition = 0;
				readEndPosition = alignment.Length - clipSizes[0] - 1;
				contigStartPosition = alignment.Position;
				contigEndPosition = alignment.GetEndPosition() - 1;
			}
		}
		else {
			readStartPosition = readPositions[0];
			readEndPosition = alignment.Length - clipSizes[1] - 1;
			contigStartPosition = genomePositions[0];
			contigEndPosition = alignment.GetEndPosition() - 1;
		}
	}

	if (readStartPosition < contigStartPosition) {
		if (readStartPosition > minAlignmentLength) {
			return false;
		}
		longResultHead->longResult[index].readStartPosition = 0;
		longResultHead->longResult[index].contigStartPosition = contigStartPosition - readStartPosition;
	}
	else {
		if (contigStartPosition > minAlignmentLength) {
			return false;
		}
		longResultHead->longResult[index].readStartPosition = readStartPosition - contigStartPosition;
		longResultHead->longResult[index].contigStartPosition = 0;
	}


	if (readLength - readEndPosition < contigLength - contigEndPosition) {
		if (readLength - readEndPosition > minAlignmentLength) {
			return false;
		}
		longResultHead->longResult[index].readEndPosition = readLength - 1;
		longResultHead->longResult[index].contigEndPosition = contigEndPosition + (readLength - readEndPosition - 1);
	}
	else {
		if (contigLength - contigEndPosition > minAlignmentLength) {
			return false;
		}
		longResultHead->longResult[index].readEndPosition = readEndPosition + contigLength - contigEndPosition - 1;
		longResultHead->longResult[index].contigEndPosition = contigLength - 1;
	}

	if (alignment.IsReverseStrand() != false) {	
		int temp = longResultHead->longResult[index].readStartPosition;
		longResultHead->longResult[index].readStartPosition = readLength - longResultHead->longResult[index].readEndPosition - 1;
		longResultHead->longResult[index].readEndPosition = readLength - temp - 1;

	}

	longResultHead->longResult[index].contigIndex = contigIndex;
	longResultHead->longResult[index].overlapLength = longResultHead->longResult[index].contigEndPosition - longResultHead->longResult[index].contigStartPosition + 1;
	longResultHead->longResult[index].orientation = !alignment.IsReverseStrand();
	longResultHead->longResult[index].quality = alignment.MapQuality;

	if (longResultHead->longResult[index].overlapLength < contigSetHead->minOverlapLength) {
		return false;
	}

	clipSizes.clear();
	readPositions.clear();
	genomePositions.clear();

	if (contigSetHead->contigSet[contigIndex].shortContig == true) {
		longResultHead->aligningShortContigResultCount++;
	}
	else {
		longResultHead->longResultCount++;
	}
	
	return true;
}

void outPutLongAlignResult(LongResultHead* longResultHead, ContigSetHead* contigSetHead, FILE* fp, FILE* fp1, long int readIndex, long int readLength) {

	for (long int i = 0; i < longResultHead->longResultCount + longResultHead->aligningShortContigResultCount - 1; i++) {
		for (long int j = i + 1; j < longResultHead->longResultCount + longResultHead->aligningShortContigResultCount; j++) {
			if (longResultHead->longResult[i].contigIndex == longResultHead->longResult[j].contigIndex) {	
				return;
			}
		}
	}
	for (long int i = 0; i < longResultHead->longResultCount + longResultHead->aligningShortContigResultCount - 1; i++) {
		for (long int j = i + 1; j < longResultHead->longResultCount + longResultHead->aligningShortContigResultCount; j++) {
			if (longResultHead->longResult[i].readStartPosition > longResultHead->longResult[j].readStartPosition) {
				int temp = longResultHead->longResult[i].readStartPosition;
				longResultHead->longResult[i].readStartPosition = longResultHead->longResult[j].readStartPosition;
				longResultHead->longResult[j].readStartPosition = temp;

				temp = longResultHead->longResult[i].readEndPosition;
				longResultHead->longResult[i].readEndPosition = longResultHead->longResult[j].readEndPosition;
				longResultHead->longResult[j].readEndPosition = temp;

				temp = longResultHead->longResult[i].contigStartPosition;
				longResultHead->longResult[i].contigStartPosition = longResultHead->longResult[j].contigStartPosition;
				longResultHead->longResult[j].contigStartPosition = temp;

				temp = longResultHead->longResult[i].contigEndPosition;
				longResultHead->longResult[i].contigEndPosition = longResultHead->longResult[j].contigEndPosition;
				longResultHead->longResult[j].contigEndPosition = temp;

				temp = longResultHead->longResult[i].contigIndex;
				longResultHead->longResult[i].contigIndex = longResultHead->longResult[j].contigIndex;
				longResultHead->longResult[j].contigIndex = temp;

				temp = longResultHead->longResult[i].overlapLength;
				longResultHead->longResult[i].overlapLength = longResultHead->longResult[j].overlapLength;
				longResultHead->longResult[j].overlapLength = temp;

				temp = longResultHead->longResult[i].orientation;
				longResultHead->longResult[i].orientation = longResultHead->longResult[j].orientation;
				longResultHead->longResult[j].orientation = temp;

				temp = longResultHead->longResult[i].quality;
				longResultHead->longResult[i].quality = longResultHead->longResult[j].quality;
				longResultHead->longResult[j].quality = temp;

			}
		}
	}

	fprintf(fp, "%d,%ld", longResultHead->longResultCount, readLength);
	fprintf(fp1, "%d,%ld", longResultHead->longResultCount + longResultHead->aligningShortContigResultCount, readLength);


	for (long int i = 0; i < longResultHead->longResultCount + longResultHead->aligningShortContigResultCount; i++) {
		if (contigSetHead->contigSet[longResultHead->longResult[i].contigIndex].shortContig != true) {
			fprintf(fp, ",%d,%ld,%ld,%ld,%ld,%d,%d", longResultHead->longResult[i].contigIndex,
				longResultHead->longResult[i].readStartPosition, longResultHead->longResult[i].readEndPosition,
				longResultHead->longResult[i].contigStartPosition, longResultHead->longResult[i].contigEndPosition,
				longResultHead->longResult[i].overlapLength, longResultHead->longResult[i].orientation);
		}
		fprintf(fp1, ",%d,%ld,%ld,%ld,%ld,%d,%d", longResultHead->longResult[i].contigIndex,
			longResultHead->longResult[i].readStartPosition, longResultHead->longResult[i].readEndPosition,
			longResultHead->longResult[i].contigStartPosition, longResultHead->longResult[i].contigEndPosition,
			longResultHead->longResult[i].overlapLength, longResultHead->longResult[i].orientation);
	}
	fprintf(fp, ",\n");
	fprintf(fp1, ",\n");


}


AligningResultHead* GetAligningResultHead(ContigSetHead* contigSetHead, char* file, int readType, char* originalShortResult, int insertSize, int maxSize, char* line) {
	
	FILE* fp;
	if ((fp = fopen(originalShortResult, "r")) == NULL) {
		printf("%s, does not exist!", originalShortResult);
		exit(0);
	}

	char* p;
	const char* split = ",";
	long int i = 0;
	long int aligningCount = 0;
	while ((fgets(line, maxSize, fp)) != NULL) {
		aligningCount++;
	}
	fclose(fp);

	AligningResultHead* aligningResultHead = (AligningResultHead*)malloc(sizeof(AligningResultHead));
	aligningResultHead->aligningResultCount = aligningCount;
	aligningResultHead->aligningResult = (AligningResult*)malloc(sizeof(AligningResult) * aligningCount);

	if ((fp = fopen(originalShortResult, "r")) == NULL) {
		printf("%s, does not exist!", originalShortResult);
		exit(0);
	}
	long int temp = 0;
	bool t = true;
	double d = 0;
	while ((fgets(line, maxSize, fp)) != NULL) {

		p = strtok(line, split);
		aligningResultHead->aligningResult[i].leftContigIndex = atoi(p);

		p = strtok(NULL, split);
		aligningResultHead->aligningResult[i].rightContigIndex = atoi(p);

		p = strtok(NULL, split);
		aligningResultHead->aligningResult[i].leftPosition = atoi(p);

		p = strtok(NULL, split);
		aligningResultHead->aligningResult[i].rightPosition = atoi(p);

		p = strtok(NULL, split);
		aligningResultHead->aligningResult[i].leftOrientation = atoi(p);

		p = strtok(NULL, split);
		aligningResultHead->aligningResult[i].rightOrientation = atoi(p);

		p = strtok(NULL, split);
		aligningResultHead->aligningResult[i].leftLength = atoi(p);

		p = strtok(NULL, split);
		aligningResultHead->aligningResult[i].rightLength = atoi(p);

		d = 0;
		if (readType == 1) {
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				d = insertSize - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].leftContigIndex].contigLength - aligningResultHead->aligningResult[i].leftPosition) - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].rightContigIndex].contigLength - aligningResultHead->aligningResult[i].rightPosition);
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				d = insertSize - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].leftContigIndex].contigLength - aligningResultHead->aligningResult[i].leftPosition) - (aligningResultHead->aligningResult[i].rightPosition + aligningResultHead->aligningResult[i].rightLength);
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				d = insertSize - (aligningResultHead->aligningResult[i].leftPosition + aligningResultHead->aligningResult[i].leftLength) - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].rightContigIndex].contigLength - aligningResultHead->aligningResult[i].rightPosition);
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				d = insertSize - (aligningResultHead->aligningResult[i].leftLength + aligningResultHead->aligningResult[i].leftPosition) - (aligningResultHead->aligningResult[i].rightLength + aligningResultHead->aligningResult[i].rightPosition);
			}
		}
		if (readType == 2) {
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				d = insertSize - (aligningResultHead->aligningResult[i].leftLength + aligningResultHead->aligningResult[i].leftPosition) - (aligningResultHead->aligningResult[i].rightLength + aligningResultHead->aligningResult[i].rightPosition);
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				d = insertSize - (aligningResultHead->aligningResult[i].leftPosition + aligningResultHead->aligningResult[i].leftLength) - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].rightContigIndex].contigLength - aligningResultHead->aligningResult[i].rightPosition);
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				d = insertSize - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].leftContigIndex].contigLength - aligningResultHead->aligningResult[i].leftPosition) - (aligningResultHead->aligningResult[i].rightPosition + aligningResultHead->aligningResult[i].rightLength);
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				d = insertSize - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].leftContigIndex].contigLength - aligningResultHead->aligningResult[i].leftPosition) - (contigSetHead->contigSet[aligningResultHead->aligningResult[i].rightContigIndex].contigLength - aligningResultHead->aligningResult[i].rightPosition);
			}
		}

		if (aligningResultHead->aligningResult[i].leftContigIndex > aligningResultHead->aligningResult[i].rightContigIndex) {
			temp = aligningResultHead->aligningResult[i].leftContigIndex;
			aligningResultHead->aligningResult[i].leftContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			aligningResultHead->aligningResult[i].rightContigIndex = temp;

			temp = aligningResultHead->aligningResult[i].leftPosition;
			aligningResultHead->aligningResult[i].leftPosition = aligningResultHead->aligningResult[i].rightPosition;
			aligningResultHead->aligningResult[i].rightPosition = temp;

			t = aligningResultHead->aligningResult[i].leftOrientation;
			aligningResultHead->aligningResult[i].leftOrientation = aligningResultHead->aligningResult[i].rightOrientation;
			aligningResultHead->aligningResult[i].rightOrientation = t;

			temp = aligningResultHead->aligningResult[i].leftLength;
			aligningResultHead->aligningResult[i].leftLength = aligningResultHead->aligningResult[i].rightLength;
			aligningResultHead->aligningResult[i].rightLength = temp;

		}

		aligningResultHead->aligningResult[i].gapDiatanceForSR = d;

		i++;
	}

	fclose(fp);
	qsort(aligningResultHead->aligningResult, aligningResultHead->aligningResultCount, sizeof(aligningResultHead->aligningResult[0]), cmp);

	long int start = 0;
	long int end = -1;
	long int leftContigIndex = -1;
	long int rightContigIndex = -1;

	long int count = 0;
	long int m = 0;
	i = 0;
	for (i = 0; i < aligningResultHead->aligningResultCount + 1; i++) {
		if (aligningResultHead->aligningResult[i].leftContigIndex == leftContigIndex && aligningResultHead->aligningResult[i].rightContigIndex == rightContigIndex) {
			end++;
		}
		else {
			if (end - start + 1 > 0) {
				for (m = start; m <= end; m++) {
					if (abs(aligningResultHead->aligningResult[m].gapDiatanceForSR) > insertSize * (1 + 0.3)) {
						aligningResultHead->aligningResult[m].leftContigIndex = -1;
					}
				}
			}
			start = i;
			end = i;
			m = 0;
			leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
		}
	}
	count = 0;
	for (long int i = 0; i < aligningResultHead->aligningResultCount; i++) {
		if (aligningResultHead->aligningResult[i].leftContigIndex == -1) {
			count++;
			continue;
		}
		else {
			aligningResultHead->aligningResult[i - count].leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			aligningResultHead->aligningResult[i - count].rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			aligningResultHead->aligningResult[i - count].leftPosition = aligningResultHead->aligningResult[i].leftPosition;
			aligningResultHead->aligningResult[i - count].rightPosition = aligningResultHead->aligningResult[i].rightPosition;
			aligningResultHead->aligningResult[i - count].leftOrientation = aligningResultHead->aligningResult[i].leftOrientation;
			aligningResultHead->aligningResult[i - count].rightOrientation = aligningResultHead->aligningResult[i].rightOrientation;
			aligningResultHead->aligningResult[i - count].leftLength = aligningResultHead->aligningResult[i].leftLength;
			aligningResultHead->aligningResult[i - count].rightLength = aligningResultHead->aligningResult[i].rightLength;
			aligningResultHead->aligningResult[i - count].gapDiatanceForSR = aligningResultHead->aligningResult[i].gapDiatanceForSR;
		}
	}
	aligningResultHead->aligningResultCount = aligningResultHead->aligningResultCount - count;


	OptimaizeShortAlignResult(aligningResultHead);

	return aligningResultHead;

}

int cmp(const void* arg1, const void* arg2) {
	AligningResult* a1 = (AligningResult*)arg1;
	AligningResult* a2 = (AligningResult*)arg2;
	if (a1->leftContigIndex != a2->leftContigIndex) {
		return a1->leftContigIndex - a2->leftContigIndex;
	}
	else if (a1->rightContigIndex != a2->rightContigIndex && a1->leftContigIndex == a2->leftContigIndex) {
		return a1->rightContigIndex - a2->rightContigIndex;
	}
	else if (a1->leftPosition != a2->leftPosition && a1->rightContigIndex == a2->rightContigIndex && a1->leftContigIndex == a2->leftContigIndex) {
		return a1->leftPosition - a2->leftPosition;
	}
	else {
		return a1->rightPosition - a2->rightPosition;
	}
}

void OptimaizeShortAlignResult(AligningResultHead* aligningResultHead) {

	int* type = (int*)malloc(sizeof(int) * 4);
	long int leftContigIndex = -1;
	long int rightContigIndex = -1;
	long int count = 0;
	long int startIndex = -1;

	long int max = 0;
	long int maxIndex = -1;

	for (long int i = 0; i < aligningResultHead->aligningResultCount + 1; i++) {	
		if (aligningResultHead->aligningResult[i].leftContigIndex == leftContigIndex && aligningResultHead->aligningResult[i].rightContigIndex == rightContigIndex) {
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				type[0]++;
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				type[1]++;
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				type[2]++;
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				type[3]++;
			}
			count++;
		}
		else {
			if (leftContigIndex != rightContigIndex) {
				if (count > 0) {
					max = 0;
					maxIndex = -1;

					for (int j = 0; j <= 3; j++) {
						if (type[j] > max) {
							max = type[j];
							maxIndex = j;
						}
					}

					for (long int j = startIndex; j < startIndex + count; j++) {
						if (aligningResultHead->aligningResult[j].leftContigIndex == leftContigIndex && aligningResultHead->aligningResult[j].rightContigIndex == rightContigIndex) {
							if (aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 0) {
								if (type[0] != max) {
									aligningResultHead->aligningResult[j].leftContigIndex = -1;
								}
							}
							if (aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 1) {
								if (type[1] != max) {
									aligningResultHead->aligningResult[j].leftContigIndex = -1;
								}
							}
							if (aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 0) {
								if (type[2] != max) {
									aligningResultHead->aligningResult[j].leftContigIndex = -1;
								}
							}
							if (aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 1) {
								if (type[3] != max) {
									aligningResultHead->aligningResult[j].leftContigIndex = -1;
								}
							}
						}
					}
				}
			}

			type[0] = 0;
			type[1] = 0;
			type[2] = 0;
			type[3] = 0;
			leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			startIndex = i;
			count = 1;
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				type[0]++;
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				type[1]++;
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 0) {
				type[2]++;
			}
			if (aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 1) {
				type[3]++;
			}
		}
	}
	count = 0;
	for (long int i = 0; i < aligningResultHead->aligningResultCount; i++) {
		if (aligningResultHead->aligningResult[i].leftContigIndex == -1) {
			count++;
			continue;
		}
		else {
			aligningResultHead->aligningResult[i - count].leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			aligningResultHead->aligningResult[i - count].rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			aligningResultHead->aligningResult[i - count].leftPosition = aligningResultHead->aligningResult[i].leftPosition;
			aligningResultHead->aligningResult[i - count].rightPosition = aligningResultHead->aligningResult[i].rightPosition;
			aligningResultHead->aligningResult[i - count].leftOrientation = aligningResultHead->aligningResult[i].leftOrientation;
			aligningResultHead->aligningResult[i - count].rightOrientation = aligningResultHead->aligningResult[i].rightOrientation;
			aligningResultHead->aligningResult[i - count].leftLength = aligningResultHead->aligningResult[i].leftLength;
			aligningResultHead->aligningResult[i - count].rightLength = aligningResultHead->aligningResult[i].rightLength;
			aligningResultHead->aligningResult[i - count].gapDiatanceForSR = aligningResultHead->aligningResult[i].gapDiatanceForSR;
		}
	}
	aligningResultHead->aligningResultCount = aligningResultHead->aligningResultCount - count;

}

double CaculateShortDistance(AligningResultHead* aligningResultHead, SimpleResultHead* simpleResultHead, long int& count, long int leftIndex, long int rightIndex, long int startIndex, int* type, int max) {
	double dis = 0;
	int num = 0;
	for (long int j = startIndex; j < startIndex + count; j++) {
		if (aligningResultHead->aligningResult[j].leftContigIndex == leftIndex && aligningResultHead->aligningResult[j].rightContigIndex == rightIndex) {
			if (aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 0) {
				if (type[0] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[0] = 0;
					num++;
				}
			}
			if (aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 1) {
				if (type[1] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[1] = 0;
					num++;
				}
			}
			if (aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 0) {
				if (type[2] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[2] = 0;
					num++;
				}
			}
			if (aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 1) {
				if (type[3] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[3] = 0;
					num++;
				}
			}
		}
		if (aligningResultHead->aligningResult[j].leftContigIndex == rightIndex && aligningResultHead->aligningResult[j].rightContigIndex == leftIndex) {
			if (aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 0) {
				if (type[3] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[3] = 0;
					num++;
				}
			}
			if (aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 1) {
				if (type[2] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[2] = 0;
					num++;
				}
			}
			if (aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 0) {
				if (type[1] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[1] = 0;
					num++;
				}
			}
			if (aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 1) {
				if (type[0] != max) {
					aligningResultHead->aligningResult[j].leftContigIndex = -1;
					type[0] = 0;
					num++;
				}
			}
		}
	}

	for (long int i = startIndex; i < startIndex + count; i++) {
		if ((aligningResultHead->aligningResult[i].leftContigIndex == leftIndex && aligningResultHead->aligningResult[i].rightContigIndex == rightIndex)) {
			if (aligningResultHead->aligningResult[i].leftContigIndex != -1) {
				dis = dis + aligningResultHead->aligningResult[i].gapDiatanceForSR;
			}
		}
	}
	count = count - num;
	dis = (dis / (double)count);

	return dis;
}


void OptimizeLongRead(ContigSetHead* contigSetHead, char* originalfile, char* line, int maxSize, FILE* file) {	
	FILE* fp;
	if ((fp = fopen(originalfile, "r")) == NULL) {
		printf("%s,does not exist!", originalfile);
		exit(0);
	}
	const char* split = ",";
	char* p;

	int maxCount = 100;
	int aligningSingle[7 * maxCount];
	int lineNum = -1;
	while ((fgets(line, maxSize, fp)) != NULL) {
		
		lineNum++;
		p = strtok(line, split);	
		int count = atoi(p);
		if (count <= 1) {
			continue;
		}
		p = strtok(NULL, split);	
		int readLength = atoi(p);

		int a = 1;

		while (a <= count * 7) {
			p = strtok(NULL, split);
			aligningSingle[a - 1] = atoi(p);
			a++;
		}

		for (int i = 0; i < count - 1; i++) {
			for (int j = i + 1; j < count; j++) {
				if (aligningSingle[i * 7 + 1] == aligningSingle[j * 7 + 1] && aligningSingle[i * 7 + 1] == 0) {	
					if (aligningSingle[i * 7 + 5] > aligningSingle[j * 7 + 5]) {	
						aligningSingle[j * 7] = -1;
					}
					if (aligningSingle[i * 7 + 5] < aligningSingle[j * 7 + 5]) {
						aligningSingle[i * 7] = -1;
					}
					if (aligningSingle[i * 7 + 5] == aligningSingle[j * 7 + 5]) {
						aligningSingle[j * 7] = -1;
						aligningSingle[i * 7] = -1;
					}
				}
				if (aligningSingle[i * 7 + 2] == aligningSingle[j * 7 + 2] && aligningSingle[i * 7 + 2] == readLength - 1) {	
					if (aligningSingle[i * 7 + 5] > aligningSingle[j * 7 + 5]) {
						aligningSingle[j * 7] = -1;
					}
					if (aligningSingle[i * 7 + 5] < aligningSingle[j * 7 + 5]) {
						aligningSingle[i * 7] = -1;
					}
					if (aligningSingle[i * 7 + 5] == aligningSingle[j * 7 + 5]) {
						aligningSingle[j * 7] = -1;
						aligningSingle[i * 7] = -1;
					}
				}
			}
		}
		int realCount = count;
		int interval = 0;
		for (int i = 0; i < count; i++) {
			if (aligningSingle[i * 7] == -1) {
				interval++;
			}
			else {
				for (int j = 0; j < 7; j++) {
					aligningSingle[(i - interval) * 7 + j] = aligningSingle[i * 7 + j];
				}
			}
		}
		count = count - interval;
		realCount = count;

		int previousIndex = -1;
		fprintf(file, "%d,", realCount);
		int gapDistance = 0;
		for (int i = 0; i < count; i++) {
			if (aligningSingle[i * 7] == -1) {
				break;
			}
			if (i == count - 1) {	
				gapDistance = 0;
			}
			else {
				gapDistance = aligningSingle[(i + 1) * 7 + 1] - aligningSingle[i * 7 + 2];	
			}

			int j = i + 1;
			bool token = false;
			for (j = i + 1; j < count - 1; j++) {
				if (aligningSingle[j * 7] == -1) {
					gapDistance = gapDistance + contigSetHead->contigSet[aligningSingle[j * 7]].contigLength + aligningSingle[(j + 1) * 7 + 1] - aligningSingle[j * 7 + 2];
				}
				else {
					token = true;
					break;
				}
			}

			if (token == false) {
				if (aligningSingle[(count - 1) * 7 + 3] == -1) {
					gapDistance = 0;
					fprintf(file, "%d,%d,%d,%d,%d,", aligningSingle[i * 7], aligningSingle[i * 7 + 3], gapDistance, aligningSingle[i * 7 + 6], aligningSingle[i * 7 + 5]);//index  contigStartPosition  gap  ori  overlap
					break;
				}
			}

			fprintf(file, "%d,%d,%d,%d,%d,", aligningSingle[i * 7], aligningSingle[i * 7 + 3], gapDistance, aligningSingle[i * 7 + 6], aligningSingle[i * 7 + 5]);
			gapDistance = 0;
		}
		fprintf(file, "\n");

	}
	fclose(fp);

}

void AddEdgeInSimpleGraph(SimpleResultHead* simpleResultHead, long int simpleLongAllCount, int leftIndex, bool leftOrientation, int rightIndex, int rightOrientation, int distance, int minOverlapLength) {
	for (long int i = 0; i < simpleLongAllCount; i++) {
		if (simpleResultHead->simpleLongResult[i].leftIndex == -1) {
			simpleResultHead->simpleLongResult[i].leftIndex = leftIndex;
			simpleResultHead->simpleLongResult[i].rightIndex = rightIndex;
			simpleResultHead->simpleLongResult[i].gapSimpleLong = distance;
			simpleResultHead->simpleLongResult[i].overlapSimpleLong = minOverlapLength;

			if (leftOrientation == 0 && rightOrientation == 0) {
				simpleResultHead->simpleLongResult[i].count1++;
			}
			if (leftOrientation == 0 && rightOrientation == 1) {
				simpleResultHead->simpleLongResult[i].count2++;
			}
			if (leftOrientation == 1 && rightOrientation == 0) {
				simpleResultHead->simpleLongResult[i].count3++;
			}
			if (leftOrientation == 1 && rightOrientation == 1) {
				simpleResultHead->simpleLongResult[i].count4++;
			}
			break;
		}
		if (simpleResultHead->simpleLongResult[i].leftIndex == leftIndex && simpleResultHead->simpleLongResult[i].leftIndex == rightIndex) {
			if (leftOrientation == 0 && rightOrientation == 0) {
				simpleResultHead->simpleLongResult[i].count1++;
			}
			if (leftOrientation == 0 && rightOrientation == 1) {
				simpleResultHead->simpleLongResult[i].count2++;
			}
			if (leftOrientation == 1 && rightOrientation == 0) {
				simpleResultHead->simpleLongResult[i].count3++;
			}
			if (leftOrientation == 1 && rightOrientation == 1) {
				simpleResultHead->simpleLongResult[i].count4++;
			}
			break;
		}
		if (simpleResultHead->simpleLongResult[i].leftIndex == rightIndex && simpleResultHead->simpleLongResult[i].leftIndex == leftIndex) {
			if (leftOrientation == 0 && rightOrientation == 0) {
				simpleResultHead->simpleLongResult[i].count4++;
			}
			if (leftOrientation == 0 && rightOrientation == 1) {
				simpleResultHead->simpleLongResult[i].count2++;
			}
			if (leftOrientation == 1 && rightOrientation == 0) {
				simpleResultHead->simpleLongResult[i].count3++;
			}
			if (leftOrientation == 1 && rightOrientation == 1) {
				simpleResultHead->simpleLongResult[i].count1++;
			}
			break;
		}
	}
}

void AddEdgeInSimpleGraph1(SimpleLongResult* simpleLongResult, long int simpleLongAllCount, int leftIndex, bool leftOrientation, int rightIndex, int rightOrientation, int distance, int minOverlapLength) {
	for (long int i = 0; i < simpleLongAllCount; i++) {
		if (simpleLongResult[i].leftIndex == -1) {
			simpleLongResult[i].leftIndex = leftIndex;
			simpleLongResult[i].rightIndex = rightIndex;
			simpleLongResult[i].gapSimpleLong = distance;
			simpleLongResult[i].overlapSimpleLong = minOverlapLength;

			if (leftOrientation == 0 && rightOrientation == 0) {
				simpleLongResult[i].count1++;
			}
			if (leftOrientation == 0 && rightOrientation == 1) {
				simpleLongResult[i].count2++;
			}
			if (leftOrientation == 1 && rightOrientation == 0) {
				simpleLongResult[i].count3++;
			}
			if (leftOrientation == 1 && rightOrientation == 1) {
				simpleLongResult[i].count4++;
			}
			break;
		}
		if (simpleLongResult[i].leftIndex == leftIndex && simpleLongResult[i].leftIndex == rightIndex) {
			if (leftOrientation == 0 && rightOrientation == 0) {
				simpleLongResult[i].count1++;
			}
			if (leftOrientation == 0 && rightOrientation == 1) {
				simpleLongResult[i].count2++;
			}
			if (leftOrientation == 1 && rightOrientation == 0) {
				simpleLongResult[i].count3++;
			}
			if (leftOrientation == 1 && rightOrientation == 1) {
				simpleLongResult[i].count4++;
			}
			break;
		}
		if (simpleLongResult[i].leftIndex == rightIndex && simpleLongResult[i].leftIndex == leftIndex) {
			if (leftOrientation == 0 && rightOrientation == 0) {
				simpleLongResult[i].count4++;
			}
			if (leftOrientation == 0 && rightOrientation == 1) {
				simpleLongResult[i].count2++;
			}
			if (leftOrientation == 1 && rightOrientation == 0) {
				simpleLongResult[i].count3++;
			}
			if (leftOrientation == 1 && rightOrientation == 1) {
				simpleLongResult[i].count1++;
			}
			break;
		}
	}
}

SimpleLongResult* GetLineIndex1(ContigSetHead* contigSetHead, bool* lineIndex, char* file, char* line, int maxSize, long int& lineCount, long int& simpleGraphNodeCount) {
	FILE* fp;
	if ((fp = fopen(file, "r")) == NULL) {
		printf("%s, does not exist!", file);
		exit(0);
	}
	char* p;
	const char* split = ",";
	simpleGraphNodeCount = 0;
	lineCount = 0;

	while ((fgets(line, maxSize, fp)) != NULL) {
		lineCount++;
		p = strtok(line, split);

		int count = atoi(p);
		if (count <= 1) {
			continue;
		}
		else {
			simpleGraphNodeCount = simpleGraphNodeCount + count - 1;
		}
	}

	SimpleLongResult* simpleLongResult = (SimpleLongResult*)malloc(sizeof(SimpleLongResult) * simpleGraphNodeCount);
	for (long int i = 0; i < simpleGraphNodeCount; i++) {
		simpleLongResult[i].leftIndex = -1;
		simpleLongResult[i].rightIndex = -1;
		simpleLongResult[i].simpleAlignCount = -1;
		simpleLongResult[i].count1 = 0;
		simpleLongResult[i].count2 = 0;
		simpleLongResult[i].count3 = 0;
		simpleLongResult[i].count4 = 0;
		simpleLongResult[i].gapSimpleLong = -1;
		simpleLongResult[i].overlapSimpleLong = -1;
	}

	while ((fgets(line, maxSize, fp)) != NULL) {
		p = strtok(line, split);
		int count = atoi(p);
		if (count <= 1) {
			continue;
		}
		p = strtok(NULL, split);
		int startContigIndex = atoi(p);
		p = strtok(NULL, split);
		int distance = atoi(p);
		p = strtok(NULL, split);
		bool startOrientation = atoi(p);
		p = strtok(NULL, split);
		int startOverlapLength = atoi(p);

		int a = 2;
		int temp = 0;
		while (a <= count) {
			p = strtok(NULL, split);
			int endContigIndex = atoi(p);
			p = strtok(NULL, split);
			int distance1 = atoi(p);
			p = strtok(NULL, split);
			bool endOrientation = atoi(p);
			p = strtok(NULL, split);
			int endOverlapLength = atoi(p);

			int min = 0;
			if (startOverlapLength > endOverlapLength) {
				min = endOverlapLength;
			}
			else {
				min = startOverlapLength;
			}

			AddEdgeInSimpleGraph1(simpleLongResult, simpleGraphNodeCount, startContigIndex, startOrientation, endContigIndex, endOrientation, distance, min);
			distance = distance1;
			startContigIndex = endContigIndex;
			startOrientation = endOrientation;
			startOverlapLength = endOverlapLength;

			a++;
		}
	}
	fclose(fp);

	return simpleLongResult;
}

SimpleResultHead* OptimaizeLongAlign(ContigSetHead* contigSetHead, char* optimizeLongResult) {
	FILE* fp;
	if ((fp = fopen(optimizeLongResult, "r")) == NULL) {
		printf("%s, does not exist!", optimizeLongResult);
		exit(0);
	}

	int maxSize = 10000;
	char* line = (char*)malloc(sizeof(char) * maxSize);
	char* p;
	const char* split = ",";
	long int i = 0;
	long int longCount = 0;
	long int lineCount = 0;
	while ((fgets(line, maxSize, fp)) != NULL) {
		lineCount++;
		p = strtok(line, split);
		int count = atoi(p);
		if (count <= 1) {
			continue;
		}
		else {
			longCount = longCount + count - 1;
		}
	}
	fclose(fp);
	SimpleResultHead* simpleResultHead = (SimpleResultHead*)malloc(sizeof(SimpleResultHead));
	simpleResultHead->simpleLongLineCount = longCount;
	SimpleLongResult* simpleLongResult = (SimpleLongResult*)malloc(sizeof(SimpleLongResult) * longCount);
	simpleResultHead->simpleLongResult = simpleLongResult;
	for (long int i = 0; i < longCount; i++) {
		simpleResultHead->simpleLongResult[i].leftIndex = -1;
		simpleResultHead->simpleLongResult[i].rightIndex = -1;
		simpleResultHead->simpleLongResult[i].leftOritation = false;
		simpleResultHead->simpleLongResult[i].rightOritation = false;
		simpleResultHead->simpleLongResult[i].gapSimpleLong = -1;
		simpleResultHead->simpleLongResult[i].overlapSimpleLong = -1;
		simpleResultHead->simpleLongResult[i].leftPosition = -1;
		simpleResultHead->simpleLongResult[i].rightPosition = -1;
	}

	if ((fp = fopen(optimizeLongResult, "r")) == NULL) {
		printf("%s, does not exist!", optimizeLongResult);
		exit(0);
	}
	i = 0;
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
		int temp = 0;
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

					simpleResultHead->simpleLongResult[i].leftIndex = startContigIndex;
					simpleResultHead->simpleLongResult[i].rightIndex = endContigIndex;
					simpleResultHead->simpleLongResult[i].leftOritation = startOrientation;
					simpleResultHead->simpleLongResult[i].rightOritation = endOrientation;
					simpleResultHead->simpleLongResult[i].gapSimpleLong = distance;
					simpleResultHead->simpleLongResult[i].overlapSimpleLong = min;
					simpleResultHead->simpleLongResult[i].leftPosition = startPosition;
					simpleResultHead->simpleLongResult[i].rightPosition = endPosition;

					if (startContigIndex > endContigIndex) {
						simpleResultHead->simpleLongResult[i].leftIndex = endContigIndex;
						simpleResultHead->simpleLongResult[i].rightIndex = startContigIndex;
						simpleResultHead->simpleLongResult[i].gapSimpleLong = distance;
						simpleResultHead->simpleLongResult[i].overlapSimpleLong = min;
						simpleResultHead->simpleLongResult[i].leftPosition = endPosition;
						simpleResultHead->simpleLongResult[i].rightPosition = startPosition;
						if (startOrientation == 0 && endOrientation == 0) {
							simpleResultHead->simpleLongResult[i].leftOritation = 1;
							simpleResultHead->simpleLongResult[i].rightOritation = 1;
						}
						if (startOrientation == 1 && endOrientation == 1) {
							simpleResultHead->simpleLongResult[i].leftOritation = 0;
							simpleResultHead->simpleLongResult[i].rightOritation = 0;
						}
						if ((startOrientation == 0 && endOrientation == 1) || (startOrientation == 1 && endOrientation == 0)) {
							simpleResultHead->simpleLongResult[i].leftOritation = startOrientation;
							simpleResultHead->simpleLongResult[i].rightOritation = endOrientation;
						}
					}

					distance = distance1;
					startContigIndex = endContigIndex;
					startOrientation = endOrientation;
					startOverlapLength = endOverlapLength;
					startPosition = endPosition;
					i++;
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
	}
	fclose(fp);

	return simpleResultHead;
}

void OptimaizeLongAlignResult(ContigSetHead* contigSetHead, char* optimizeLongResult, SimpleResultHead* simpleResultHead) {

	FILE* fp;
	if ((fp = fopen(optimizeLongResult, "r")) == NULL) {
		printf("%s, does not exist!", optimizeLongResult);
		exit(0);
	}

	int maxSize = 10000;
	char* line = (char*)malloc(sizeof(char) * maxSize);
	char* p;
	const char* split = ",";
	long int simpleLongAllCount = 0;
	long int lineCount = 0;
	while ((fgets(line, maxSize, fp)) != NULL) {	
		lineCount++;
		p = strtok(line, split);
		int count = atoi(p);
		if (count <= 1) {
			continue;
		}
		else {
			simpleLongAllCount = simpleLongAllCount + count - 1;
		}
	}
	simpleResultHead->simpleLongLineCount = simpleLongAllCount;

	bool* lineIndex = (bool*)malloc(sizeof(bool) * lineCount);
	for (long int i = 0; i < lineCount; i++) {
		lineIndex[i] = false;
	}

	SimpleLongResult* simpleLongResult = (SimpleLongResult*)malloc(sizeof(SimpleLongResult) * simpleLongAllCount);
	simpleResultHead->simpleLongResult = simpleLongResult;
	for (long int i = 0; i < simpleLongAllCount; i++) {
		simpleResultHead->simpleLongResult[i].leftIndex = -1;
		simpleResultHead->simpleLongResult[i].rightIndex = -1;
		simpleResultHead->simpleLongResult[i].simpleAlignCount = -1;
		simpleResultHead->simpleLongResult[i].count1 = 0;
		simpleResultHead->simpleLongResult[i].count2 = 0;
		simpleResultHead->simpleLongResult[i].count3 = 0;
		simpleResultHead->simpleLongResult[i].count4 = 0;
		simpleResultHead->simpleLongResult[i].gapSimpleLong = -1;
		simpleResultHead->simpleLongResult[i].overlapSimpleLong = -1;
	}
	fclose(fp);
	if ((fp = fopen(optimizeLongResult, "r")) == NULL) {
		printf("%s, does not exist!", optimizeLongResult);
		exit(0);
	}
	while ((fgets(line, maxSize, fp)) != NULL) {
		p = strtok(line, split);
		int count = atoi(p);
		if (count <= 1) {
			continue;
		}
		p = strtok(NULL, split);
		int startContigIndex = atoi(p);
		p = strtok(NULL, split);
		int distance = atoi(p);
		p = strtok(NULL, split);
		bool startOrientation = atoi(p);
		p = strtok(NULL, split);
		int startOverlapLength = atoi(p);
		int a = 2;
		int temp = 0;
		while (a <= count) {
			p = strtok(NULL, split);
			int endContigIndex = atoi(p);
			p = strtok(NULL, split);
			int distance1 = atoi(p);
			p = strtok(NULL, split);
			bool endOrientation = atoi(p);
			p = strtok(NULL, split);
			int endOverlapLength = atoi(p);
			int min = 0;
			if (startOverlapLength > endOverlapLength) {
				min = endOverlapLength;
			}
			else {
				min = startOverlapLength;
			}
			AddEdgeInSimpleGraph(simpleResultHead, simpleLongAllCount, startContigIndex, startOrientation, endContigIndex, endOrientation, distance, min);
			distance = distance1;
			startContigIndex = endContigIndex;
			startOrientation = endOrientation;
			startOverlapLength = endOverlapLength;
			a++;
		}
	}
	fclose(fp);


	qsort(simpleResultHead->simpleLongResult, simpleResultHead->simpleLongLineCount, sizeof(simpleResultHead->simpleLongResult[0]), cmpl);

	deleteDuplicateInSimpleLong(simpleResultHead);

}

int cmpl(const void* arg1, const void* arg2) {
	SimpleLongResult* a1 = (SimpleLongResult*)arg1;
	SimpleLongResult* a2 = (SimpleLongResult*)arg2;
	//return a1->leftContigIndex - a2->leftContigIndex;
	if (a1->leftIndex != a2->leftIndex) {
		return a1->leftIndex - a2->leftIndex;
	}
	else {
		return a1->rightIndex - a2->rightIndex;
	}
}

SimpleResultHead* deleteDuplicateInSimpleLong(SimpleResultHead* simpleResultHead) {

	int num = 0;
	int count = 0;
	int i = 0;
	int t = 0;

	for (int j = 1; j < simpleResultHead->simpleLongLineCount + 1; j++) {
		if (simpleResultHead->simpleLongResult[i].leftIndex == simpleResultHead->simpleLongResult[j].leftIndex
			&& simpleResultHead->simpleLongResult[i].rightIndex == simpleResultHead->simpleLongResult[j].rightIndex) {
			t++;
			num++;
			mixOrientation(simpleResultHead, i, j);

		}
		else {
			cout << t << endl;
			simpleResultHead->simpleLongResult[i].simpleAlignCount = t + 1;
			simpleResultHead->simpleLongResult[++i] = simpleResultHead->simpleLongResult[j];
			t = 0;

		}
	}
	simpleResultHead->simpleLongLineCount = simpleResultHead->simpleLongLineCount - num;
	return simpleResultHead;
}

void mixOrientation(SimpleResultHead* simpleResultHead, int start, int end) {
	simpleResultHead->simpleLongResult[start].count1 += simpleResultHead->simpleLongResult[end].count1;
	simpleResultHead->simpleLongResult[start].count2 += simpleResultHead->simpleLongResult[end].count2;
	simpleResultHead->simpleLongResult[start].count3 += simpleResultHead->simpleLongResult[end].count3;
	simpleResultHead->simpleLongResult[start].count4 += simpleResultHead->simpleLongResult[end].count4;

	simpleResultHead->simpleLongResult[start].gapSimpleLong = int((simpleResultHead->simpleLongResult[start].gapSimpleLong + simpleResultHead->simpleLongResult[end].gapSimpleLong) / 2);

	simpleResultHead->simpleLongResult[start].overlapSimpleLong = min(simpleResultHead->simpleLongResult[start].overlapSimpleLong, simpleResultHead->simpleLongResult[end].overlapSimpleLong);
}





#endif