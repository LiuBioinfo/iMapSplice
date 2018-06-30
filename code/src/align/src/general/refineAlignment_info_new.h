// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef REFALIGNMENT_INFO_NEW_H
#define REFALIGNMENT_INFO_NEW_H

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
#include <sstream>

using namespace std;

class Align2Refine_Info
{
private:
	int chrNameInt;
	int chrMapPos;
	vector<Jump_Code> alignJumpCodeVec;	
	vector< pair<int,int> > SJposPairVec;
	vector< int > SJindexInAlignInferJuncVec;
	string otherStr;
public:
};


class PeReadAlign2Refine_Info
{
private:
	PE_Read_Info peReadInfo;
	bool pairedOrNotBool;

	vector<Align2Refine_Info* > align2RefineInfoVec_Nor1;
	vector<Align2Refine_Info* > align2RefineInfoVec_Rcm1;
	vector<Align2Refine_Info* > align2RefineInfoVec_Nor2;
	vector<Align2Refine_Info* > align2RefineInfoVec_Rcm2;
public:
	PeReadAlign2Refine_Info()
	{

	}

	bool end1OrEnd2_bool(int tmpFlag)
	{
		if(tmpFlag & 0x40)
			return true;
		else
			return false;
	}

	bool forOrRev_bool(int tmpFlag)
	{
		if(tmpFlag & 0x10)
			return true;
		else
			return false;
	}

	bool bothEndsMapped(int tmpFlag)
	{
		if(tmpFlag & 0x2)
			return true;
		else
			return false;
	}

	void initiateWithAlignStrVec(vector<string>& peReadAlignmentStrVec)
	// only worsk for MPS3 results;
	// assumptions:
	//	1. all alignments paired
	//	2. alignments are pair-wisely presented ... each pair is listed together
	{
		string firstPeReadAlignmentStr = peReadAlignmentStrVec[0];
		string secondPeReadAlignmentStr = peReadAlignmentStrVec[1];
		this->parseAndInitiatePeReadInfo_pushBackFirstPeAlignment(
			firstPeReadAlignmentStr, secondPeReadAlignmentStr);
		for(int tmp = 2; tmp < peReadAlignmentStrVec.size(); tmp ++)
		{
			string tmpPeReadAlignmentStr_1 = peReadAlignmentStrVec[tmp];
			tmp++;
			string tmpPeReadAlignmentStr_2 = peReadAlignmentStrVec[tmp];
			this->parseAndPushBackPeAlignment(tmpPeReadAlignmentStr_1,
				tmpPeReadAlignmentStr_2);
		}
	}

};


#endif