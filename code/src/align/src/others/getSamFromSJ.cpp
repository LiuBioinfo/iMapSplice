// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
//#include <omp.h>
#include "read_block_test.h"
#include "bwtmap_info.h"
#include "DoubleAnchorScore.h"
#include "sbndm.h"
#include "splice_info.h"
#include "jumpCode_info.h"

#define SPLICE_MIN_LENGTH 50
using namespace std;

typedef map<int, int> SpliceWeightMap; 
typedef map <int, SpliceWeightMap > SpliceJunctionHash; 
typedef SpliceJunctionHash::iterator SpliceJunctionHashIter;
typedef SpliceWeightMap::iterator SpliceWeightMapIter;

SpliceJunctionHash spliceJunction[22];
//SpliceJunctionHash spliceJunctionAcceptor[22]; //<spliceStartPos, spliceEndPos>
SpliceJunctionHashIter iter;
SpliceWeightMapIter weightMapIter;
int groundTruthJunctionNum = 0;
int groundTruthSpliceNum = 0;
int groundTruthDeletionNum = 0;
//int deletionNum = 0;
int tofindJunctionNum = 0;
int tofindSpliceNum = 0;
int tofindDeletionNum = 0;

int covertStringToInt(string chrName)
{
	if(chrName == "chr1")
		return 0;
	else if(chrName == "chr2")
		return 1;
	else if(chrName == "chr3")
		return 2;
	else if(chrName == "chr4")
		return 3;
	else if(chrName == "chr5")
		return 4;
	else if(chrName == "chr6")
		return 5;
	else if(chrName == "chr7")
		return 6;
	else if(chrName == "chr8")
		return 7;
	else if(chrName == "chr9")
		return 8;
	else if(chrName == "chr10")
		return 9;
	else if(chrName == "chr11")
		return 10;
	else if(chrName == "chr12")
		return 11;
	else if(chrName == "chr13")
		return 12;
	else if(chrName == "chr14")
		return 13;
	else if(chrName == "chr15")
		return 14;
	else if(chrName == "chr16")
		return 15;
	else if(chrName == "chr17")
		return 16;
	else if(chrName == "chr18")
		return 17;
	else if(chrName == "chr19")
		return 18;
	else if(chrName == "chrX")
		return 19;
	else if(chrName == "chrY")
		return 20;
	else if(chrName == "chrM")
		return 21;
	else
		return 22;
}

bool insertSpliceJunction2Hash(SpliceJunctionHash* spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
	int chrInt, int spliceStartPos, int spliceEndPos) 
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////insert to spliceJunction///////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	groundTruthJunctionNum ++;
	if((spliceEndPos - spliceStartPos) < SPLICE_MIN_LENGTH)
	{
		groundTruthDeletionNum ++;
		return true;
	}
	groundTruthSpliceNum ++;
	iter = spliceJunction[chrInt].find(spliceStartPos); // to find spliceStartPos in Hash
	if(iter == spliceJunction[chrInt].end()) 
	{
		SpliceWeightMap* newSpliceWeightMap = new SpliceWeightMap;
		(*newSpliceWeightMap).insert(pair<unsigned int, unsigned int> (spliceEndPos, 1));  // insert spliceEndPos and weight(1) into Map
		spliceJunction[chrInt].insert(pair<unsigned int, SpliceWeightMap> (spliceStartPos, (*newSpliceWeightMap)));
	}
	else
	{
		weightMapIter = (iter->second).find(spliceEndPos);
		if(weightMapIter == (iter->second).end())
		{
			(iter->second).insert(pair<unsigned int, unsigned int> (spliceEndPos, 1));
		}
		else
		{
			(weightMapIter->second) ++;
		}
	}
	return true;
}

bool findSpliceJunctionInHash(SpliceJunctionHash* spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
	int chrInt, int spliceStartPos, int spliceEndPos) 
{
	//tofindSpliceNum ++;
	if(chrInt > 21)
		return false;
	iter = spliceJunction[chrInt].find(spliceStartPos);
	if(iter == spliceJunction[chrInt].end()) 
	{
		return false;
	}
	else
	{
		weightMapIter = (iter->second).find(spliceEndPos);
		if(weightMapIter == (iter->second).end())
		{
			return false;
		}
		else
		{
			return true;	
		}
	}
}


int main(int argc, char**argv)
{
	if(argc < 3)
	{
		cout << "Executable <InputSam> <InputSJ> <OutputSAM>" << endl;
		exit(0);
	}
	
	char *inputSamFile = argv[1]; string inputSamFileStr = inputSamFile;
	char *inputSJfile = argv[2]; string inputSJfileStr = inputSJfile;
	char *outputSamFile = argv[3]; string outputSamFileStr = outputSamFile;
	
	ofstream OutputSamFile_ofs(outputSamFile);

	string outputSamFileStr_shortAnchor = outputSamFileStr + ".shortAnchor";
	string outputSamFileStr_indel = outputSamFileStr + ".indel";
	string outputSamFileStr_multiJunc = outputSamFileStr + ".multiJunc";
	string outputSamFileStr_singleJunc = outputSamFileStr + ".singleJunc";
	string outputSamFileStr_others = outputSamFileStr + ".others";

	ofstream OutputSamFile_shortAnchor_ofs(outputSamFileStr_shortAnchor.c_str());
	ofstream OutputSamFile_indel_ofs(outputSamFileStr_indel.c_str());
	ofstream OutputSamFile_multiJunc_ofs(outputSamFileStr_multiJunc.c_str());
	ofstream OutputSamFile_singleJunc_ofs(outputSamFileStr_singleJunc.c_str());
	ofstream OutputSamFile_others_ofs(outputSamFileStr_others.c_str());

	int shortAnchorJuncNum = 0;
	int indelJuncNum = 0;
	int multiJuncNum = 0;
	int singleJuncNum = 0;
	int othersJuncNum = 0;
	int allJuncNum = 0;

	FILE *fp_in1 = fopen(inputSJfile, "r");
	string entryString;
	int tabLocation1;
	int tabLocation2;
	int tabLocation3;
	char entry[500];
	int chrInt;
	int spliceStartPos;
	int spliceEndPos;
	string chrIntString;
	string spliceStartPosString;
	string spliceEndPosString;

	//int groundTruthSpliceNum = 0;
	//int toFindSpliceNum = 0;

	fgets(entry, sizeof(entry), fp_in1);
	while(!feof(fp_in1))
	{
		fgets(entry, sizeof(entry), fp_in1);
		//groundTruthSpliceNum ++;
		entryString = entry;
		tabLocation1 = entryString.find('\t', 0);
		tabLocation2 = entryString.find('\t', tabLocation1+1);
		tabLocation3 = entryString.find('\t', tabLocation2+1);
		chrIntString = entryString.substr(0, tabLocation1);
		spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
		spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
		chrInt = covertStringToInt(chrIntString);
		if(chrInt > 21)
			continue;
		spliceStartPos = atoi(spliceStartPosString.c_str());
		spliceEndPos = atoi(spliceEndPosString.c_str());
		bool insertion = insertSpliceJunction2Hash(spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
			chrInt, spliceStartPos, spliceEndPos);
		if(!insertion)
		{
			cout << "insertion failed " << endl;
			cout << chrInt << " " << spliceStartPos << " " << spliceEndPos << endl;
		}
	}
	
	FILE *fp;
	fp = fopen(inputSamFile, "r");
	char ch[500];
	char read_name[100], flag[3], rname[10],pos[20],mapq[3],cigar[50],rnext[10];
	char pnext[20], tlen[4],seq[200],qual[200],others[500];
	char NM[20],IH[20],HI[20],YI[20],YH[20],XS[20],XF[20],ZF[50];	
	string readNameStr, readSequenceStr, cigarStringStr;
	int M_loc_1 = 0, M_loc_2 = 0, M_loc_3 = 0;
	int N_loc_1 = 0, N_loc_2 = 0; 
	string smallExonLengthStr;	
	int smallExonLengthInt;
	
	while(!feof(fp))
	{   
		fgets(ch,sizeof(ch),fp);
		if(feof(fp))
			break;
			
		sscanf(ch,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,NM,IH,HI,YI,YH,XS,XF,ZF);    

		string chrNameStr = rname;
		int chrNameStrInt = covertStringToInt(chrNameStr);
		if(chrNameStrInt > 21)
			continue;

		int chrPosInt = atoi(pos);		
		string cigarStringStr = cigar;

		JumpCode_Info* jumpCodeInfo = new JumpCode_Info();
		jumpCodeInfo->jumpCodeStr2jumpCodeVec(cigarStringStr);
		jumpCodeInfo->jumpCodeVec2SJforCompareVec(chrPosInt);

		int SJvecSize = (jumpCodeInfo->SJforCompareVec).size();
		int JumpCodeSize = (jumpCodeInfo->cigarStringJumpCode).size();

		bool shortAnchorBool = false;
		bool indelBool = false;
		bool multiAnchorBool = false;
		bool singleAnchorBool = false;
		bool othersBool = false;		

		int headJumpCodeLength = (jumpCodeInfo->cigarStringJumpCode)[0].len;
		int tailJumpCodeLength = (jumpCodeInfo->cigarStringJumpCode)[JumpCodeSize-1].len;

		shortAnchorBool = ((headJumpCodeLength < 10)||(tailJumpCodeLength < 10));
		indelBool = ((cigarStringStr.find("I") != cigarStringStr.npos) 
						|| (cigarStringStr.find("D") != cigarStringStr.npos));
		multiAnchorBool = (SJvecSize > 1);
		singleAnchorBool = (SJvecSize == 1)&&(!shortAnchorBool)&&(!indelBool);
		othersBool = (!( shortAnchorBool || indelBool || multiAnchorBool || singleAnchorBool ));

		for(int tmpSJnum = 0; tmpSJnum < SJvecSize;
			tmpSJnum ++)
		{	

			int tmpSpliceStartPos = (jumpCodeInfo->SJforCompareVec)[tmpSJnum].first;
			int tmpSpliceEndPos = (jumpCodeInfo->SJforCompareVec)[tmpSJnum].second;
			bool findSplice = false;
			for(int tmp1 = 0; tmp1 <= 0; tmp1++)
			{
				if(findSplice)
					break;
				for(int tmp2 = 0; tmp2 <= 0; tmp2++)
				{
					findSplice = findSpliceJunctionInHash(spliceJunction, //SpliceJunctionHash* spliceJunctionAcceptor,
						chrNameStrInt, tmpSpliceStartPos+tmp1, tmpSpliceEndPos+tmp2);
					if(findSplice)
					{
						break;			
					}
				}
			}	
			if(findSplice)
			{
				allJuncNum ++;
				//OutputSamFile_ofs << ch << endl;
				OutputSamFile_ofs << "SJ: " << chrNameStr << " " 
					<< tmpSpliceStartPos << " " << tmpSpliceEndPos << endl 
					<< ch << endl;

				if(shortAnchorBool)
				{
					shortAnchorJuncNum ++;
					OutputSamFile_shortAnchor_ofs << "SJ: " << chrNameStr << " " 
						<< tmpSpliceStartPos << " " << tmpSpliceEndPos << endl 
						<< ch << endl;
				}
				if(indelBool)
				{
					indelJuncNum ++;
					OutputSamFile_indel_ofs << "SJ: " << chrNameStr << " " 
						<< tmpSpliceStartPos << " " << tmpSpliceEndPos << endl 
						<< ch << endl;
				}
				if(multiAnchorBool)
				{
					multiJuncNum ++;
					OutputSamFile_multiJunc_ofs << "SJ: " << chrNameStr << " " 
						<< tmpSpliceStartPos << " " << tmpSpliceEndPos << endl 
						<< ch << endl;
				}
				if(singleAnchorBool)
				{
					singleJuncNum ++;
					OutputSamFile_singleJunc_ofs << "SJ: " << chrNameStr << " " 
						<< tmpSpliceStartPos << " " << tmpSpliceEndPos << endl 
						<< ch << endl;				
				}
				if(othersBool)
				{
					othersJuncNum ++;
					OutputSamFile_others_ofs << "SJ: " << chrNameStr << " " 
						<< tmpSpliceStartPos << " " << tmpSpliceEndPos << endl 
						<< ch << endl;
				}

			}	
		}

	}	
	fclose(fp);

	cout << "shortAnchorJuncNum: " << shortAnchorJuncNum << endl
		<< "indelJuncNum: " << indelJuncNum << endl
		<< "multiJuncNum: " << multiJuncNum << endl
		<< "singleJuncNum: " << singleJuncNum << endl 
		<< "othersJuncNum: " << othersJuncNum << endl
		<< "allJuncNum: " << allJuncNum << endl;

	return 0;	
}

