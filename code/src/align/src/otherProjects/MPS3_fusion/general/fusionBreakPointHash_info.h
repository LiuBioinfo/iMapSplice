// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FUSIONBREAKPOINTHASH_INFO_H
#define FUSIONBREAKPOINTHASH_INFO_H

#include "fusionBreakPoint_info.h"

using namespace std;

typedef map<int, int> SecondBreakPoint2FusionBreakPointHashIndexMap;
typedef map<int, SecondBreakPoint2FusionBreakPointHashIndexMap> First2secondBreakPointMap;

typedef map<int, set<int> > BreakPointAreaHash;
typedef BreakPointAreaHash::iterator BreakPointAreaHashIter;

class FusionBreakPointHash_Info
{
private:
	vector < vector < First2secondBreakPointMap > > breakPointMapVecVec_for;
	vector < vector < First2secondBreakPointMap > > breakPointMapVecVec_rev;
	vector < FusionBreakPoint_Info* > breakPointInfoVec;

	vector < BreakPointAreaHash > breakPointAreaHashVec_start;
	vector < BreakPointAreaHash > breakPointAreaHashVec_end;

	int areaSize;
public:
	FusionBreakPointHash_Info()
	{
		areaSize = 1000;
	}

	void clearSupportNum()
	{
		for(int tmp = 0; tmp < breakPointInfoVec.size(); tmp++)
			breakPointInfoVec[tmp]->clearSupportNum();
	}

	void clearSupportNum_encompassing()
	{
		for(int tmp = 0; tmp < breakPointInfoVec.size(); tmp++)
			breakPointInfoVec[tmp]->clearSupportNum_encompassing();
	}

	void addSupportNum_withIndexInBreakPointInfoVec_encompassing(int toAddSupportNum, int index)
	{
		breakPointInfoVec[index]->addSupportNum_encompassing(toAddSupportNum);
	}

	void addSupportNum_withIndexInBreakPointInfoVec(int toAddSupportNum, int index)
	{
		breakPointInfoVec[index]->addSupportNum(toAddSupportNum);
	}

	void addSupportNum_withChrNameStrBreakPointPos_encompassing(
		int toAddSupportNum, string& tmpChrNameStr_gene1, string& tmpChrNameStr_gene2, 
		int tmpBreakPointPos_gene1, int tmpBreakPointPos_gene2, Index_Info* indexInfo)
	{
		int tmpChrNameInt_gene1 = indexInfo->convertStringToInt(tmpChrNameStr_gene1);
		int tmpChrNameInt_gene2 = indexInfo->convertStringToInt(tmpChrNameStr_gene2);
		int tmpIndexInFusionBreakPointInfoVec 
			= this->searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_gene1, tmpChrNameInt_gene2,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2);
		if(tmpIndexInFusionBreakPointInfoVec < 0)
			return;
		else
		{
			breakPointInfoVec[tmpIndexInFusionBreakPointInfoVec]->addSupportNum_encompassing(toAddSupportNum);
		}		
	}

	void updateAnchorSizeFromAnotherFusionBreakPointHashInfo(
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo)
	{
		int tmpBreakPointVecSize = breakPointInfoVec.size();
		for(int tmp = 0; tmp < tmpBreakPointVecSize; tmp++)
		{
			int tmpAnchorLength_1 = tmpFusionBreakPointHashInfo->returnAnchorLengthWithIndex_1(tmp);
			int tmpAnchorLength_2 = tmpFusionBreakPointHashInfo->returnAnchorLengthWithIndex_2(tmp);
			breakPointInfoVec[tmp]->updateAnchorSize(tmpAnchorLength_1, tmpAnchorLength_2);
		}
	}

	void updateXMfromAnotherFusionBreakPointHashInfo(
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo)
	{
		int tmpBreakPointVecSize = breakPointInfoVec.size();
		for(int tmp = 0; tmp < tmpBreakPointVecSize; tmp++)
		{
			int tmpXMmin = tmpFusionBreakPointHashInfo->returnXMminWithIndex(tmp);
			int tmpXMmax = tmpFusionBreakPointHashInfo->returnXMmaxWithIndex(tmp);
			breakPointInfoVec[tmp]->updateAnchorSize(tmpXMmin, tmpXMmax);
		}
	}

	void addSupportNumFromAnotherFusionBreakPointHashInfo(
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo)
	{
		int tmpBreakPointVecSize = breakPointInfoVec.size();
		for(int tmp = 0; tmp < tmpBreakPointVecSize; tmp++)
		{
			int tmpSupportNum = tmpFusionBreakPointHashInfo->returnSupportNumWithIndex(tmp);
			breakPointInfoVec[tmp]->addSupportNum(tmpSupportNum);
		}
	}

	void addSupportNumFromAnotherFusionBreakPointHashInfo_encompassing(
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo)
	{
		int tmpBreakPointVecSize = breakPointInfoVec.size();
		for(int tmp = 0; tmp < tmpBreakPointVecSize; tmp++)
		{
			int tmpSupportNum_encompassing 
				= tmpFusionBreakPointHashInfo->returnSupportNumWithIndex_encompassing(tmp);
			breakPointInfoVec[tmp]->addSupportNum_encompassing(tmpSupportNum_encompassing);
		}
	}

	void searchFusionBreakPointFromAreaHashWithinRange(
		int tmpChrNameInt, int leftPos, int rightPos,
		bool thisGene_1_or_2_bool, vector<int>& foundFusionBreakPointPosVec)
	{
		int leftPosAreaIndex = leftPos/areaSize;
		int rightPosAreaIndex = rightPos/areaSize;
		// cout << "leftPosAreaIndex: " << leftPosAreaIndex << endl;
		// cout << "rightPosAreaIndex: " << rightPosAreaIndex << endl;
		if(thisGene_1_or_2_bool)
		{	
			for(int tmpAreaIndex = leftPosAreaIndex;
				tmpAreaIndex <= rightPosAreaIndex; tmpAreaIndex++)
			{
				//cout << "tmpAreaIndex: " << tmpAreaIndex << endl;
				BreakPointAreaHashIter tmpAreaIter 
					= breakPointAreaHashVec_start[tmpChrNameInt].find(tmpAreaIndex);
				if(tmpAreaIter == breakPointAreaHashVec_start[tmpChrNameInt].end())
				{}
				else
				{
					for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin();
						tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
					{
						int tmpBreakPointPos = (*tmpSetIter);
						if((tmpBreakPointPos >= leftPos)&&(tmpBreakPointPos <= rightPos))
							foundFusionBreakPointPosVec.push_back(tmpBreakPointPos);
					}
				}
			}
		}
		else
		{
			for(int tmpAreaIndex = leftPosAreaIndex;
				tmpAreaIndex <= rightPosAreaIndex; tmpAreaIndex++)
			{
				BreakPointAreaHashIter tmpAreaIter 
					= breakPointAreaHashVec_end[tmpChrNameInt].find(tmpAreaIndex);
				if(tmpAreaIter == breakPointAreaHashVec_end[tmpChrNameInt].end())
				{}
				else
				{
					for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin();
						tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
					{
						int tmpBreakPointPos = (*tmpSetIter);
						if((tmpBreakPointPos >= leftPos)&&(tmpBreakPointPos <= rightPos))
							foundFusionBreakPointPosVec.push_back(tmpBreakPointPos);
					}
				}
			}
		}
	}

	void searchAndReturnFusionIndexVec_bothEnds(
		int tmpChrNameInt_gene1, int leftPos_gene1, int rightPos_gene1,
		int tmpChrNameInt_gene2, int leftPos_gene2, int rightPos_gene2, 
		vector<int>& tmpFusionIndexVec)
	{
		vector<int> tmpFusionIndexVec_searchWithRangeOfGene1;
		vector<int> tmpFusionIndexVec_searchWithRangeOfGene2;
		this->searchAndReturnFusionIndexVec_oneEnd(tmpChrNameInt_gene1, leftPos_gene1, 
			rightPos_gene1, true, tmpFusionIndexVec_searchWithRangeOfGene1);
		this->searchAndReturnFusionIndexVec_oneEnd(tmpChrNameInt_gene2, leftPos_gene2,
			rightPos_gene2, false, tmpFusionIndexVec_searchWithRangeOfGene2);
		for(int tmp = 0; tmp < tmpFusionIndexVec_searchWithRangeOfGene1.size(); tmp++)
		{
			int tmpIndex = tmpFusionIndexVec_searchWithRangeOfGene1[tmp];
			tmpFusionIndexVec.push_back(tmpIndex);
		}
		for(int tmp = 0; tmp < tmpFusionIndexVec_searchWithRangeOfGene2.size(); tmp++)
		{
			int tmpIndex = tmpFusionIndexVec_searchWithRangeOfGene2[tmp];
			int tmpFusionIndexVecSize = tmpFusionIndexVec.size();
			bool tmpIndex_existed_bool = false;
			for(int tmp2 = 0; tmp2 < tmpFusionIndexVecSize; tmp2 ++)
			{
				int tmpIndex_existed = tmpFusionIndexVec[tmp2];
				if(tmpIndex == tmpIndex_existed)
				{
					tmpIndex_existed_bool = true;
					break;
				}
			}
			if(!tmpIndex_existed_bool)
				tmpFusionIndexVec.push_back(tmpIndex);
		}
	}

	void searchAndReturnFusionIndexVec_oneEnd(int tmpChrNameInt, int leftPos, int rightPos,
		bool thisGene_1_or_2_bool, vector<int>& fusionIndexVec)
	{
		//cout << "start to do searchAndReturnFusionIndexVec_oneEnd ..." << endl;
		vector<int> foundFusionBreakPointPosVec;
		vector<int> theOtherGeneChrNameIntVec; 
		vector<int> theOtherGeneBreakPointVec;
		this->searchFusionBreakPointFromAreaHashWithinRange_returnFusionBreakPointPairVec(
			tmpChrNameInt, leftPos, rightPos, thisGene_1_or_2_bool, 
			foundFusionBreakPointPosVec, theOtherGeneChrNameIntVec, theOtherGeneBreakPointVec);
		//cout << "foundFusionBreakPointPosVec.size(): " << foundFusionBreakPointPosVec.size() << endl;
		for(int tmp = 0; tmp < foundFusionBreakPointPosVec.size(); tmp++)
		{
			int tmpChrNameInt_thisGene = tmpChrNameInt;
			int tmpBreakPointPos_thisGene = foundFusionBreakPointPosVec[tmp];
			int tmpChrNameInt_theOtherGene = theOtherGeneChrNameIntVec[tmp];
			int tmpBreakPointPos_theOtherGene = theOtherGeneBreakPointVec[tmp];
			int tmpFusionIndex;
			if(thisGene_1_or_2_bool)
				tmpFusionIndex = this->searchAndReturnIndexInBreakPointInfoVec_for(
					tmpChrNameInt_thisGene, tmpChrNameInt_theOtherGene,
					tmpBreakPointPos_thisGene, tmpBreakPointPos_theOtherGene);
			else
				tmpFusionIndex = this->searchAndReturnIndexInBreakPointInfoVec_for(
					tmpChrNameInt_theOtherGene, tmpChrNameInt_thisGene,
					tmpBreakPointPos_theOtherGene, tmpBreakPointPos_thisGene);
			if(tmpFusionIndex >= 0)
				fusionIndexVec.push_back(tmpFusionIndex);
			else
			{
				cout << "error in searchAndReturnFusionIndexVec, fusionIndex: " << tmpFusionIndex << endl;
				exit(1);

			}
		}
	}

	void searchFusionBreakPointFromAreaHashWithinRange_returnFusionBreakPointPairVec(
		int tmpChrNameInt, int leftPos, int rightPos,
		bool thisGene_1_or_2_bool, vector<int>& foundFusionBreakPointPosVec,
		vector<int>& theOtherGeneChrNameIntVec, vector<int>& theOtherGeneBreakPointVec)
	{
		vector<int> tmpFoundFusionBreakPointPosVec;
		this->searchFusionBreakPointFromAreaHashWithinRange(tmpChrNameInt, 
			leftPos, rightPos, thisGene_1_or_2_bool, tmpFoundFusionBreakPointPosVec);
		for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec.size(); tmp++)
		{
			int tmpFoundFusionBreakPointPos = tmpFoundFusionBreakPointPosVec[tmp];
			vector<int> tmpTheOtherGeneChrNameIntVec;
			vector<int> tmpTheOtherGeneBreakPointPosVec;
			this->returnTheOtherFusionGeneBreakPoint(tmpChrNameInt, tmpFoundFusionBreakPointPos,
				thisGene_1_or_2_bool, tmpTheOtherGeneChrNameIntVec, tmpTheOtherGeneBreakPointPosVec);
			for(int tmpFusionPair = 0; tmpFusionPair < tmpTheOtherGeneChrNameIntVec.size(); tmpFusionPair ++)
			{
				int tmpTheOtherGeneChrNameInt = tmpTheOtherGeneChrNameIntVec[tmpFusionPair];
				int tmpTheOtherGeneBreakPointPos = tmpTheOtherGeneBreakPointPosVec[tmpFusionPair];
				foundFusionBreakPointPosVec.push_back(tmpFoundFusionBreakPointPos);
				theOtherGeneChrNameIntVec.push_back(tmpTheOtherGeneChrNameInt);
				theOtherGeneBreakPointVec.push_back(tmpTheOtherGeneBreakPointPos);
			}
		}
	}

	void returnTheOtherFusionGeneBreakPoint(
		int tmpGeneChrNameInt, int tmpBreakPoint, bool thisGene_1_or_2_bool,
		vector<string>& thisGeneStrandVec, vector<int>& theOtherGeneChrNameIntVec,
		vector<int>& theOtherGeneBreakPointVec, vector<string>& theOtherGeneStrandVec)
	{
		// cout << "function returnTheOtherFusionGeneBreakPoint starts ......" << endl;
		// cout << "tmpGeneChrNameInt: " << tmpGeneChrNameInt << endl;
		// cout << "tmpBreakPoint: " << tmpBreakPoint << endl;
		int tmpChromNum = breakPointMapVecVec_for.size();
		//cout << "tmpGeneChrNameInt: " << tmpGeneChrNameInt << endl;
		//cout << "tmpBreakPoint: " << tmpBreakPoint << endl;
		if(thisGene_1_or_2_bool)
		{
			for(int tmpChr = 0; tmpChr < tmpChromNum; tmpChr++)
			{
				First2secondBreakPointMap::iterator tmpFirstMapIter
					= (breakPointMapVecVec_for[tmpGeneChrNameInt])[tmpChr].find(tmpBreakPoint);
				if(tmpFirstMapIter == (breakPointMapVecVec_for[tmpGeneChrNameInt])[tmpChr].end())// not found
				{
					//cout << "error not found breakPoint in searchTheOtherFusionBreakPoint from FusionBreakPointHash_Info.." << endl;
					//exit(1);
				}
				else // found 
				{
					for(SecondBreakPoint2FusionBreakPointHashIndexMap::iterator tmp2ndTypeMapIter
						= (tmpFirstMapIter->second).begin();
						tmp2ndTypeMapIter != (tmpFirstMapIter->second).end();
						tmp2ndTypeMapIter ++)
					{
						int tmpTheOtherGeneBreakPointPos = tmp2ndTypeMapIter->first;
						int tmpFusionIndex = tmp2ndTypeMapIter->second;
						string tmpFirstGeneStrand = breakPointInfoVec[tmpFusionIndex]->returnStrand_1();
						string tmpSecondGeneStrand = breakPointInfoVec[tmpFusionIndex]->returnStrand_2();
						thisGeneStrandVec.push_back(tmpFirstGeneStrand);
						theOtherGeneChrNameIntVec.push_back(tmpChr);
						theOtherGeneBreakPointVec.push_back(tmpTheOtherGeneBreakPointPos);
						theOtherGeneStrandVec.push_back(tmpSecondGeneStrand);
					}
				}
			}
		}
		else
		{
			for(int tmpChr = 0; tmpChr < tmpChromNum; tmpChr++)
			{
				//cout << "tmpChr: " << tmpChr << endl;
				First2secondBreakPointMap::iterator tmpFirstMapIter
					= (breakPointMapVecVec_rev[tmpGeneChrNameInt])[tmpChr].find(tmpBreakPoint);
				if(tmpFirstMapIter == (breakPointMapVecVec_rev[tmpGeneChrNameInt])[tmpChr].end())// not found
				{
					//cout << "breakPos in gene2 not found" << endl;
					//cout << "error not found breakPoint in searchTheOtherFusionBreakPoint from FusionBreakPointHash_Info.." << endl;
					//exit(1);
				}
				else // found 
				{
					for(SecondBreakPoint2FusionBreakPointHashIndexMap::iterator tmp2ndTypeMapIter
						= (tmpFirstMapIter->second).begin();
						tmp2ndTypeMapIter != (tmpFirstMapIter->second).end();
						tmp2ndTypeMapIter ++)
					{
						int tmpTheOtherGeneBreakPointPos = tmp2ndTypeMapIter->first;
						int tmpFusionIndex = tmp2ndTypeMapIter->second;
						string tmpFirstGeneStrand = breakPointInfoVec[tmpFusionIndex]->returnStrand_1();
						string tmpSecondGeneStrand = breakPointInfoVec[tmpFusionIndex]->returnStrand_2();
						thisGeneStrandVec.push_back(tmpSecondGeneStrand);
						theOtherGeneChrNameIntVec.push_back(tmpChr);
						theOtherGeneBreakPointVec.push_back(tmpTheOtherGeneBreakPointPos);
						theOtherGeneStrandVec.push_back(tmpFirstGeneStrand);
					}
				}
			}
		}
	}

	void returnTheOtherFusionGeneBreakPoint(
		int tmpGeneChrNameInt, int tmpBreakPoint, bool thisGene_1_or_2_bool,
		vector<int>& theOtherGeneChrNameIntVec, vector<int>& theOtherGeneBreakPointVec)
	{
		vector<string> tmpThisGeneStrandVec;
		vector<string> tmpTheOtherGeneStrandVec;
		this->returnTheOtherFusionGeneBreakPoint(tmpGeneChrNameInt, tmpBreakPoint, 
			thisGene_1_or_2_bool, tmpThisGeneStrandVec, theOtherGeneChrNameIntVec, 
			theOtherGeneBreakPointVec, tmpTheOtherGeneStrandVec);
	}

	void searchFusionBreakPointWithinRangeForBothTwoGenes(
		int tmpChrNameInt_gene1, int leftMostPos_gene1, int rightMostPos_gene1,
		int tmpChrNameInt_gene2, int leftMostPos_gene2, int rightMostPos_gene2,
		string& specificFusionStrand,
		vector<int>& foundFusionBreakPointPosVec_gene1,
		vector<int>& foundFusionBreakPointPosVec_gene2)
	{
		//cout << "function searchFusionBreakPointWithinRangeForBothTwoGenes starts ......" << endl;
		vector<int> tmpFoundFusionBreakPointPosVec_gene1;
		this->searchFusionBreakPointFromAreaHashWithinRange(
			tmpChrNameInt_gene1, leftMostPos_gene1, rightMostPos_gene1,
			true, tmpFoundFusionBreakPointPosVec_gene1);
		// cout << "tmpFoundFusionBreakPointPosVec_gene1.size(): " << tmpFoundFusionBreakPointPosVec_gene1.size() << endl;
		// for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1.size(); tmp++)
		// 	cout << "tmpFoundFusionBreakPointPosInGene1: " << tmpFoundFusionBreakPointPosVec_gene1[tmp] << endl;
		for(int tmp = 0; tmp < tmpFoundFusionBreakPointPosVec_gene1.size(); tmp++)
		{
			vector<string> tmpFoundGeneStrandVec_gene1;
			vector<int> tmpFoundGeneChrNameIntVec_gene2;
			vector<int> tmpFoundGeneBreakPointVec_gene2;
			vector<string> tmpFoundGeneStrandVec_gene2;
			this->returnTheOtherFusionGeneBreakPoint(tmpChrNameInt_gene1, 
				tmpFoundFusionBreakPointPosVec_gene1[tmp], true,
				tmpFoundGeneStrandVec_gene1, tmpFoundGeneChrNameIntVec_gene2,
				tmpFoundGeneBreakPointVec_gene2, tmpFoundGeneStrandVec_gene2);
			//cout << "tmpFoundGeneStrandVec_gene1.size(): " << tmpFoundGeneStrandVec_gene1.size() << endl;
			for(int tmp2 = 0; tmp2 < tmpFoundGeneStrandVec_gene1.size(); tmp2++)
			{	
				//cout << "tmp2: " << tmp2 << endl;
				int tmpFoundGeneChrNameInt_gene2 = tmpFoundGeneChrNameIntVec_gene2[tmp2];
				int tmpFoundGeneBreakPoint_gene2 = tmpFoundGeneBreakPointVec_gene2[tmp2];
				//cout << "tmpFoundGeneChrNameInt_gene2: " << tmpFoundGeneChrNameInt_gene2 << endl;
				//cout << "tmpFoundGeneBreakPoint_gene2: " << tmpFoundGeneBreakPoint_gene2 << endl; 
				string tmpFoundFusionStrand 
					= tmpFoundGeneStrandVec_gene1[tmp2] + tmpFoundGeneStrandVec_gene2[tmp2];
				//cout << "tmpFoundFusionStrand: " << tmpFoundFusionStrand << endl;
				if((tmpFoundFusionStrand == specificFusionStrand)
					&&(tmpFoundGeneChrNameInt_gene2 == tmpChrNameInt_gene2)
					&&(tmpFoundGeneBreakPoint_gene2 <= rightMostPos_gene2)
					&&(tmpFoundGeneBreakPoint_gene2 >= leftMostPos_gene2))
				{
					foundFusionBreakPointPosVec_gene1.push_back(tmpFoundFusionBreakPointPosVec_gene1[tmp]);
					foundFusionBreakPointPosVec_gene2.push_back(tmpFoundGeneBreakPoint_gene2);
				}
			}
		}

	}

	void returnTheOtherFusionBreakPoint_afterCheckingStrand(
		int tmpGeneChrNameInt, int tmpBreakPoint, bool thisGene_1_or_2_bool,
		string& thisGeneStrand, string& theOtherGeneStrand,
		vector<int>& theOtherGeneChrNameIntVec,
		vector<int>& theOtherGeneBreakPointVec)
	{
		vector<string> thisGeneStrandVec_inter;
		vector<int> theOtherGeneChrNameIntVec_inter;
		vector<int> theOtherGeneBreakPointVec_inter;
		vector<string> theOtherGeneStrandVec_inter;
		this->returnTheOtherFusionGeneBreakPoint(
			tmpGeneChrNameInt, tmpBreakPoint, thisGene_1_or_2_bool,
			thisGeneStrandVec_inter, theOtherGeneChrNameIntVec_inter,
			theOtherGeneBreakPointVec_inter, theOtherGeneStrandVec_inter);
		for(int tmp = 0; tmp < thisGeneStrandVec_inter.size(); tmp++)
		{
			string tmpThisGeneStrand = thisGeneStrandVec_inter[tmp];
			string tmpTheOtherGeneStrand = theOtherGeneStrandVec_inter[tmp];
			if((thisGeneStrand == tmpThisGeneStrand)&&
				(theOtherGeneStrand == tmpTheOtherGeneStrand))
			{
				theOtherGeneChrNameIntVec.push_back(theOtherGeneChrNameIntVec_inter[tmp]);
				theOtherGeneBreakPointVec.push_back(theOtherGeneBreakPointVec_inter[tmp]);
			}	
		}
	}

	void generateFusionBreakPointHashInfo_fromFuionJuncFileVec(
		vector<string>& tmpFusionBreakPointFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < tmpFusionBreakPointFileVec.size(); tmp++)
		{
			string tmpFusionBreakPointFile = tmpFusionBreakPointFileVec[tmp];
			this->generateFusionBreakPointHashInfo_fromFuionJuncFile(
				tmpFusionBreakPointFile, indexInfo);
		}
	}

	void generateFusionBreakPointHashInfo_fromNonStrandedFuionJuncFileVec(
		vector<string>& tmpFusionBreakPointFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < tmpFusionBreakPointFileVec.size(); tmp++)
		{
			string tmpFusionBreakPointFile = tmpFusionBreakPointFileVec[tmp];
			this->generateFusionBreakPointHashInfo_fromNonStrandedFuionJuncFile(
				tmpFusionBreakPointFile, indexInfo);
		}
	}

	void generateFusionBreakPointHashInfo_fromFuionJuncFile(
		string& tmpFusionBreakPointFile, Index_Info* indexInfo)
	{
		ifstream fusionJunc_ifs(tmpFusionBreakPointFile.c_str());
		while(!fusionJunc_ifs.eof())
		{
			string tmpFusionJuncStr;
			getline(fusionJunc_ifs, tmpFusionJuncStr);
			//cout << "tmpFusionJuncStr: " << tmpFusionJuncStr << endl;
			if(//fusionJunc_ifs.eof() || 
				(tmpFusionJuncStr == ""))
				break;
			//cout << "tmpFusionJuncStr: " << tmpFusionJuncStr << endl;
			int tmpChrNameInt_1, tmpChrNameInt_2;
			int tmpBreakPoint_1, tmpBreakPoint_2;
			string tmpStrand_1, tmpStrand_2;
			string tmpFlankString;
			int tmpAnchorLength_1, tmpAnchorLength_2;
			int tmpSupNum, tmpSupNumStr_encompassing;
			int tmpXMmin, tmpXMmax;
			string tmpFusionCaseStr;
			int tmpMapRangeMax_1, tmpMapRangeMax_2;
			//string tmpOtherStr;
			this->getFusionJuncInfoFromFusionJuncStr(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpFlankString,
				tmpAnchorLength_1, tmpAnchorLength_2,
				tmpSupNum, tmpSupNumStr_encompassing,
				tmpXMmin, tmpXMmax, //tmpOtherStr,
				tmpFusionCaseStr, tmpMapRangeMax_1, tmpMapRangeMax_2,				
				tmpFusionJuncStr, indexInfo);
			// cout << "tmpChrNameInt_1: " << tmpChrNameInt_1 << endl;
			// cout << "tmpChrNameInt_2: " << tmpChrNameInt_2 << endl;
			// cout << "tmpBreakPoint_1: " << tmpBreakPoint_1 << endl;
			// cout << "tmpBreakPoint_2: " << tmpBreakPoint_2 << endl;
			this->insertNewFusionBreakPoint(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpSupNum, tmpSupNumStr_encompassing, tmpFlankString, 
				tmpAnchorLength_1, tmpAnchorLength_2, tmpXMmin, tmpXMmax,
				tmpFusionCaseStr, tmpMapRangeMax_1, tmpMapRangeMax_2, indexInfo);
		}
		fusionJunc_ifs.close();
	}

	void generateFusionBreakPointHashInfo_fromNonStrandedFuionJuncFile(
		string& tmpFusionBreakPointFile, Index_Info* indexInfo)
	{
		ifstream fusionJunc_ifs(tmpFusionBreakPointFile.c_str());
		while(!fusionJunc_ifs.eof())
		{
			string tmpFusionJuncStr;
			getline(fusionJunc_ifs, tmpFusionJuncStr);
			//cout << "tmpFusionJuncStr: " << tmpFusionJuncStr << endl;
			if(//fusionJunc_ifs.eof() || 
				(tmpFusionJuncStr == ""))
				break;
			//cout << "tmpFusionJuncStr: " << tmpFusionJuncStr << endl;
			int tmpChrNameInt_1, tmpChrNameInt_2;
			int tmpBreakPoint_1, tmpBreakPoint_2;
			string tmpStrand_1, tmpStrand_2;
			string tmpFlankString;
			int tmpAnchorLength_1, tmpAnchorLength_2;
			int tmpSupNum, tmpSupNumStr_encompassing;
			int tmpXMmin, tmpXMmax;
			string tmpFusionCaseStr;
			int tmpMapRangeMax_1, tmpMapRangeMax_2;
			//string tmpOtherStr;
			this->getFusionJuncInfoFromFusionJuncStr(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpFlankString,
				tmpAnchorLength_1, tmpAnchorLength_2,
				tmpSupNum, tmpSupNumStr_encompassing,
				tmpXMmin, tmpXMmax, //tmpOtherStr,
				tmpFusionCaseStr, tmpMapRangeMax_1, tmpMapRangeMax_2,				
				tmpFusionJuncStr, indexInfo);
			// cout << "tmpChrNameInt_1: " << tmpChrNameInt_1 << endl;
			// cout << "tmpChrNameInt_2: " << tmpChrNameInt_2 << endl;
			// cout << "tmpBreakPoint_1: " << tmpBreakPoint_1 << endl;
			// cout << "tmpBreakPoint_2: " << tmpBreakPoint_2 << endl;
			if(!((tmpStrand_1 == "N")&&(tmpStrand_2 == "N")))
			{
				cout << "this is not nonStranded fusion junc" << endl;
				cout << "tmpStrand_1: " << tmpStrand_1 << endl;
				cout << "tmpStrand_2: " << tmpStrand_2 << endl;
				exit(1);
			}
			else
			{
				//cout << "stranded fusion junc found " << endl;
				if(tmpFusionCaseStr == "1,4,")
				{
					tmpStrand_1 = "+";
					tmpStrand_2 = "+";
				}
				else if(tmpFusionCaseStr == "2,5,")
				{
					tmpStrand_1 = "+";
					tmpStrand_2 = "+";
				}
				else if(tmpFusionCaseStr == "7,8,")
				{
					tmpStrand_1 = "+";
					tmpStrand_2 = "-";
				}
				else if(tmpFusionCaseStr == "10,11,")
				{
					tmpStrand_1 = "-";
					tmpStrand_2 = "+";
				}
				else
				{
					cout << "error in tmpFusionCaseStr: " << tmpFusionCaseStr << endl;
					exit(1);
				}
			}

			this->insertNewFusionBreakPoint(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpSupNum, tmpSupNumStr_encompassing, tmpFlankString, 
				tmpAnchorLength_1, tmpAnchorLength_2, tmpXMmin, tmpXMmax,
				tmpFusionCaseStr, tmpMapRangeMax_1, tmpMapRangeMax_2, indexInfo);
		}
		fusionJunc_ifs.close();
	}

	void getFusionJuncInfoFromFusionJuncStr(
		int& detectedChrNameInt_1, int& detectedChrNameInt_2,
		int& detectedBreakPoint_1, int& detectedBreakPoint_2,
		string& detectedFusionJuncStrand_1, string& detectedFusionJuncStrand_2,
		string& detectedFusionJuncFlankString,
		int& detectedAnchorLength_1, int& detectedAnchorLength_2,
		int& detectedSupNum, int& detectedSupNum_encompassing,
		int& detectedXMmin, int& detectedXMmax, //string& detectedFusionJuncOtherStr,
		string& tmpFusionCaseStr, int& tmpMapRangeMax_1, int& tmpMapRangeMax_2,
		string& tmpFusionBreakPointInfoStr, Index_Info* indexInfo)
	{
		vector<string> fusionJuncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 15; tmp++)
		{
			int tabLoc = tmpFusionBreakPointInfoStr.find("\t", startLoc);
			string tmpDetectedFusionJuncField = tmpFusionBreakPointInfoStr.substr(startLoc, tabLoc-startLoc);
			fusionJuncFieldVec.push_back(tmpDetectedFusionJuncField);
			startLoc = tabLoc + 1;
		}
		int nextTabLoc = tmpFusionBreakPointInfoStr.find("\t", startLoc);
		if(nextTabLoc == string::npos)
			fusionJuncFieldVec.push_back(tmpFusionBreakPointInfoStr.substr(startLoc));
		else
		 	fusionJuncFieldVec.push_back(tmpFusionBreakPointInfoStr.substr(startLoc, nextTabLoc-startLoc));
		string tmpChrNameStr_1 = fusionJuncFieldVec[0];
		string tmpChrNameStr_2 = fusionJuncFieldVec[1];
		string tmpChrPosStr_1 = fusionJuncFieldVec[2];
		string tmpChrPosStr_2 = fusionJuncFieldVec[3];
		detectedChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
		detectedChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
		detectedBreakPoint_1 = atoi(tmpChrPosStr_1.c_str());
		detectedBreakPoint_2 = atoi(tmpChrPosStr_2.c_str());

		detectedFusionJuncStrand_1 = fusionJuncFieldVec[4];
		detectedFusionJuncStrand_2 = fusionJuncFieldVec[5];
		detectedFusionJuncFlankString = fusionJuncFieldVec[6];
		string tmpAnchorLengthStr_1 = fusionJuncFieldVec[7];
		string tmpAnchorLengthStr_2 = fusionJuncFieldVec[8];
		detectedAnchorLength_1 = atoi(tmpAnchorLengthStr_1.c_str());
		detectedAnchorLength_2 = atoi(tmpAnchorLengthStr_2.c_str());

		string tmpSupNumStr = fusionJuncFieldVec[9];
		detectedSupNum = atoi(tmpSupNumStr.c_str());
		string tmpSupNumStr_encompassing = fusionJuncFieldVec[10];
		detectedSupNum_encompassing = atoi(tmpSupNumStr_encompassing.c_str());

		string tmpXMminStr = fusionJuncFieldVec[11];
		string tmpXMmaxStr = fusionJuncFieldVec[12];
		detectedXMmin = atoi(tmpXMminStr.c_str());
		detectedXMmax = atoi(tmpXMmaxStr.c_str());

		tmpFusionCaseStr = fusionJuncFieldVec[13];
		string tmpMapRangeMax_1_str = fusionJuncFieldVec[14];
		string tmpMapRangeMax_2_str = fusionJuncFieldVec[15];
		tmpMapRangeMax_1 = atoi(tmpMapRangeMax_1_str.c_str());
		tmpMapRangeMax_2 = atoi(tmpMapRangeMax_2_str.c_str());
	}

	void insertNewFusionBreakPoint_toBreakPointAreaHash(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2)
	{
		int breakPointAreaNO_start = (int)(tmpBreakPoint_1/areaSize);
		int breakPointAreaNO_end = (int)(tmpBreakPoint_2/areaSize);
		
		BreakPointAreaHashIter foundAreaIter;

		foundAreaIter 
			= breakPointAreaHashVec_start[tmpChrNameInt_1].find(breakPointAreaNO_start);
		if(foundAreaIter == breakPointAreaHashVec_start[tmpChrNameInt_1].end())
		{
			set<int> newBreakPointSet;
			newBreakPointSet.insert(tmpBreakPoint_1);
			breakPointAreaHashVec_start[tmpChrNameInt_1].insert(pair<int, set<int> >(
				breakPointAreaNO_start, newBreakPointSet));
		}
		else
		{
			if( (foundAreaIter->second).find(tmpBreakPoint_1) 
				== (foundAreaIter->second).end() )
				(foundAreaIter->second).insert(tmpBreakPoint_1);
			else
			{}
		}

		foundAreaIter 
			= breakPointAreaHashVec_end[tmpChrNameInt_2].find(breakPointAreaNO_end);
		if(foundAreaIter == breakPointAreaHashVec_end[tmpChrNameInt_2].end())
		{
			set<int> newBreakPointSet;
			newBreakPointSet.insert(tmpBreakPoint_2);
			breakPointAreaHashVec_end[tmpChrNameInt_2].insert(pair<int, set<int> >(
				breakPointAreaNO_end, newBreakPointSet));
		}
		else
		{
			if( (foundAreaIter->second).find(tmpBreakPoint_2) 
				== (foundAreaIter->second).end() )
				(foundAreaIter->second).insert(tmpBreakPoint_2);
			else
			{}
		}					
	}

	void insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, int tmpSupNum, int tmpSupNumStr_encompassing)
	{
		int tmpIndexInBreakPointInfoVec 
			= searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2);
		if(tmpIndexInBreakPointInfoVec < 0)
		{	
			int toAddBreakPointIndexInVec = breakPointInfoVec.size();
			FusionBreakPoint_Info* tmpFusionBreakPointInfo 
				= new FusionBreakPoint_Info();
			tmpFusionBreakPointInfo->initiateWithChrNameBreakPoint(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, 
				tmpStrand_1, tmpStrand_2, tmpSupNum, tmpSupNumStr_encompassing);
			breakPointInfoVec.push_back(tmpFusionBreakPointInfo);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, toAddBreakPointIndexInVec);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_rev(
				tmpChrNameInt_2, tmpChrNameInt_1,
				tmpBreakPoint_2, tmpBreakPoint_1, toAddBreakPointIndexInVec);
		}
		else
		{
			breakPointInfoVec[tmpIndexInBreakPointInfoVec]->addSupportNum(tmpSupNum);
		}
	}

	void insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, Index_Info* indexInfo)
	{
		int tmpIndexInBreakPointInfoVec 
			= searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2);
		if(tmpIndexInBreakPointInfoVec < 0)
		{	
			int toAddBreakPointIndexInVec = breakPointInfoVec.size();
			FusionBreakPoint_Info* tmpFusionBreakPointInfo 
				= new FusionBreakPoint_Info();
			tmpFusionBreakPointInfo->initiateWithChrNameBreakPoint(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, 
				tmpStrand_1, tmpStrand_2, 
				tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString,
				tmpAnchorLength_1, tmpAnchorLength_2);
			breakPointInfoVec.push_back(tmpFusionBreakPointInfo);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, toAddBreakPointIndexInVec);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_rev(
				tmpChrNameInt_2, tmpChrNameInt_1,
				tmpBreakPoint_2, tmpBreakPoint_1, toAddBreakPointIndexInVec);
		}
		else
		{
			breakPointInfoVec[tmpIndexInBreakPointInfoVec]->addSupportNum(tmpSupNum);
		}
	}

	/*
	void insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		int tmpXMmin, int tmpXMmax, Index_Info* indexInfo)
	{
		int tmpIndexInBreakPointInfoVec 
			= searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2);
		if(tmpIndexInBreakPointInfoVec < 0)
		{	
			int toAddBreakPointIndexInVec = breakPointInfoVec.size();
			FusionBreakPoint_Info* tmpFusionBreakPointInfo 
				= new FusionBreakPoint_Info();
			tmpFusionBreakPointInfo->initiateWithChrNameBreakPoint(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, 
				tmpStrand_1, tmpStrand_2, 
				tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString,
				tmpAnchorLength_1, tmpAnchorLength_2,
				tmpXMmin, tmpXMmax);
			breakPointInfoVec.push_back(tmpFusionBreakPointInfo);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, toAddBreakPointIndexInVec);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_rev(
				tmpChrNameInt_2, tmpChrNameInt_1,
				tmpBreakPoint_2, tmpBreakPoint_1, toAddBreakPointIndexInVec);
		}
		else
		{
			breakPointInfoVec[tmpIndexInBreakPointInfoVec]->addSupportNum(tmpSupNum);
		}
	}*/

	void insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		int tmpXMmin, int tmpXMmax, Index_Info* indexInfo)
	{
		int tmpIndexInBreakPointInfoVec 
			= searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2);
		if(tmpIndexInBreakPointInfoVec < 0)
		{	
			int toAddBreakPointIndexInVec = breakPointInfoVec.size();
			FusionBreakPoint_Info* tmpFusionBreakPointInfo 
				= new FusionBreakPoint_Info();
			tmpFusionBreakPointInfo->initiateWithChrNameBreakPoint(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, 
				tmpStrand_1, tmpStrand_2, 
				tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString,
				tmpAnchorLength_1, tmpAnchorLength_2,
				tmpXMmin, tmpXMmax);
			breakPointInfoVec.push_back(tmpFusionBreakPointInfo);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, toAddBreakPointIndexInVec);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_rev(
				tmpChrNameInt_2, tmpChrNameInt_1,
				tmpBreakPoint_2, tmpBreakPoint_1, toAddBreakPointIndexInVec);
		}
		else
		{
			breakPointInfoVec[tmpIndexInBreakPointInfoVec]->addSupportNum(tmpSupNum);
		}
	}

	void insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		int tmpXMmin, int tmpXMmax, string& tmpFusionCaseStr,
		int tmpMapRangeMax_1, int tmpMapRangeMax_2, Index_Info* indexInfo)
	{
		int tmpIndexInBreakPointInfoVec 
			= searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2);
		if(tmpIndexInBreakPointInfoVec < 0)
		{	
			int toAddBreakPointIndexInVec = breakPointInfoVec.size();
			FusionBreakPoint_Info* tmpFusionBreakPointInfo 
				= new FusionBreakPoint_Info();
			tmpFusionBreakPointInfo->initiateWithChrNameBreakPoint(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, 
				tmpStrand_1, tmpStrand_2, 
				tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString,
				tmpAnchorLength_1, tmpAnchorLength_2,
				tmpXMmin, tmpXMmax, tmpFusionCaseStr, tmpMapRangeMax_1, tmpMapRangeMax_2);
			breakPointInfoVec.push_back(tmpFusionBreakPointInfo);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_for(
				tmpChrNameInt_1, tmpChrNameInt_2,
				tmpBreakPoint_1, tmpBreakPoint_2, toAddBreakPointIndexInVec);
			this->insertNewFusionBreakPoint_toBreakPointMapVecVec_rev(
				tmpChrNameInt_2, tmpChrNameInt_1,
				tmpBreakPoint_2, tmpBreakPoint_1, toAddBreakPointIndexInVec);
		}
		else
		{
			breakPointInfoVec[tmpIndexInBreakPointInfoVec]->addSupportNum(tmpSupNum);
		}
	}

	void insertNewFusionBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, 
		int tmpSupNumStr_encompassing)
	{
		this->insertNewFusionBreakPoint_toBreakPointAreaHash(
				tmpChrNameInt_1, tmpChrNameInt_2, 
				tmpBreakPoint_1, tmpBreakPoint_2);
		this->insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
				tmpChrNameInt_1, tmpChrNameInt_2, 
				tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, 
				tmpSupNum, tmpSupNumStr_encompassing);
	}

	void insertNewFusionBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2, 
		int tmpBreakPoint_1, int tmpBreakPoint_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		Index_Info* indexInfo)
	{
		this->insertNewFusionBreakPoint_toBreakPointAreaHash(
				tmpChrNameInt_1, tmpChrNameInt_2, 
				tmpBreakPoint_1, tmpBreakPoint_2);
		this->insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
				tmpChrNameInt_1, tmpChrNameInt_2, 
				tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, 
				tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString, 
				tmpAnchorLength_1, tmpAnchorLength_2, indexInfo);
	}

	void insertNewFusionBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2, 
		int tmpBreakPoint_1, int tmpBreakPoint_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		int tmpXMmin, int tmpXMmax,
		Index_Info* indexInfo)
	{
		this->insertNewFusionBreakPoint_toBreakPointAreaHash(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2);
		this->insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString, tmpAnchorLength_1, tmpAnchorLength_2, 
				tmpXMmin, tmpXMmax, indexInfo);
	}	

	void insertNewFusionBreakPoint(
		string& tmpChrNameStr_1, string& tmpChrNameStr_2, 
		int tmpBreakPoint_1, int tmpBreakPoint_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		int tmpXMmin, int tmpXMmax, 
		Index_Info* indexInfo)
	{
		int tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1); 
		int tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
		this->insertNewFusionBreakPoint_toBreakPointAreaHash(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2);
		this->insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString, tmpAnchorLength_1, tmpAnchorLength_2, 
				tmpXMmin, tmpXMmax, indexInfo);
	}


	void insertNewFusionBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2, 
		int tmpBreakPoint_1, int tmpBreakPoint_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		int tmpXMmin, int tmpXMmax,
		string& tmpFusionCaseStr,
		int tmpMapRangeMax_1, int tmpMapRangeMax_2,
		Index_Info* indexInfo)
	{
		this->insertNewFusionBreakPoint_toBreakPointAreaHash(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2);
		this->insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString, tmpAnchorLength_1, tmpAnchorLength_2, 
				tmpXMmin, tmpXMmax, tmpFusionCaseStr, tmpMapRangeMax_1, tmpMapRangeMax_2, indexInfo);
	}	

	void insertNewFusionBreakPoint(
		string& tmpChrNameStr_1, string& tmpChrNameStr_2, 
		int tmpBreakPoint_1, int tmpBreakPoint_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing, string& tmpFlankString, 
		int tmpAnchorLength_1, int tmpAnchorLength_2, 
		int tmpXMmin, int tmpXMmax, 
		string& tmpFusionCaseStr,
		int tmpMapRangeMax_1, int tmpMapRangeMax_2,
		Index_Info* indexInfo)
	{
		int tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1); 
		int tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
		this->insertNewFusionBreakPoint_toBreakPointAreaHash(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2);
		this->insertNewFusionBreakPoint_toBreakPointMapVecVec_forAndRev(
				tmpChrNameInt_1, tmpChrNameInt_2, tmpBreakPoint_1, tmpBreakPoint_2,
				tmpStrand_1, tmpStrand_2, tmpSupNum, tmpSupNumStr_encompassing,
				tmpFlankString, tmpAnchorLength_1, tmpAnchorLength_2, 
				tmpXMmin, tmpXMmax, tmpFusionCaseStr, tmpMapRangeMax_1, tmpMapRangeMax_2, indexInfo);
	}

	void insertNewFusionBreakPoint(
		string& tmpChrNameStr_1, string& tmpChrNameStr_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2,
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupNum, int tmpSupNumStr_encompassing,
		Index_Info* indexInfo)
	{
		int tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
		int tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
		this->insertNewFusionBreakPoint(tmpChrNameInt_1, tmpChrNameInt_2,
			tmpBreakPoint_1, tmpBreakPoint_2, tmpStrand_1, tmpStrand_2, 
			tmpSupNum, tmpSupNumStr_encompassing);
	}	

	void initiateWithChromNum(int chromNum)
	{
		areaSize = 1000;
		for(int tmpChrom = 0; tmpChrom < chromNum; tmpChrom ++)
		{
			vector < First2secondBreakPointMap > tmpBreakPointMapVec_for;
			vector < First2secondBreakPointMap > tmpBreakPointMapVec_rev;
			for(int tmpChrom2 = 0; tmpChrom2 < chromNum; tmpChrom2 ++)
			{
				First2secondBreakPointMap tmpFirst2secondBreakPointMap_for;
				tmpBreakPointMapVec_for.push_back(tmpFirst2secondBreakPointMap_for);
				First2secondBreakPointMap tmpFirst2secondBreakPointMap_rev;
				tmpBreakPointMapVec_rev.push_back(tmpFirst2secondBreakPointMap_rev);
			}
			breakPointMapVecVec_for.push_back(tmpBreakPointMapVec_for);
			breakPointMapVecVec_rev.push_back(tmpBreakPointMapVec_rev);

			BreakPointAreaHash tmpBreakPointAreaHash_start;
			BreakPointAreaHash tmpBreakPointAreaHash_end;
			breakPointAreaHashVec_start.push_back(tmpBreakPointAreaHash_start);
			breakPointAreaHashVec_end.push_back(tmpBreakPointAreaHash_end);
		}
	}

	void detectAlterFusionJunc(int offset)
	{
		for(int tmp = 0; tmp < breakPointInfoVec.size(); tmp++)
		{
			int tmpChrNameInt_1 = breakPointInfoVec[tmp]->returnChrNameInt_1();
			int tmpChrNameInt_2 = breakPointInfoVec[tmp]->returnChrNameInt_2();
			int tmpBreakPoint_1 = breakPointInfoVec[tmp]->returnBreakPoint_1();
			int tmpBreakPoint_2 = breakPointInfoVec[tmp]->returnBreakPoint_2();
			int leftPos_1 = tmpBreakPoint_1 - offset;
			int rightPos_1 = tmpBreakPoint_1 + offset;
			int leftPos_2 = tmpBreakPoint_2 - offset;
			int rightPos_2 = tmpBreakPoint_2 + offset;
			//cout << endl << "tmpFusionJunc: " << endl;
			//cout << "tmpChrNameInt_1: " << tmpChrNameInt_1 << endl;
			//cout << "tmpBreakPoint_1: " << tmpBreakPoint_1 << endl;
			//cout << "tmpChrNameInt_2: " << tmpChrNameInt_2 << endl;
			//cout << "tmpBreakPoint_2: " << tmpBreakPoint_2 << endl;
			vector<int> tmpAlterFusionIndexVec_1stGeneShared;
			vector<int> tmpAlterFusionIndexVec_2ndGeneShared;
			this->searchAndReturnFusionIndexVec_oneEnd(tmpChrNameInt_1, 
				leftPos_1, rightPos_1, true, tmpAlterFusionIndexVec_1stGeneShared);
			this->searchAndReturnFusionIndexVec_oneEnd(tmpChrNameInt_2,
				leftPos_2, rightPos_2, false, tmpAlterFusionIndexVec_2ndGeneShared);
			//cout << "start to copyAlterFusionIndexVec_exceptThisFusionIndex ..." << endl;		
			breakPointInfoVec[tmp]->copyAlterFusionIndexVec_exceptThisFusionIndex(tmp,
				tmpAlterFusionIndexVec_1stGeneShared, tmpAlterFusionIndexVec_2ndGeneShared);
			//cout << "end of copyAlterFusionIndexVec_exceptThisFusionIndex ..." << endl;
		}
		//cout << "end of detecting alterFusionJunc ..." << endl;
	}

	string returnAlterFusionStr(int tmpFusionIndex, Index_Info* indexInfo)
	{
		string tmpAlterFusionStr = "";
		int alterFusionNum_1stEndShared = breakPointInfoVec[tmpFusionIndex]->returnAlterFusionNum_1stEndShared();
		int alterFusionNum_2ndEndShared = breakPointInfoVec[tmpFusionIndex]->returnAlterFusionNum_2ndEndShared();
		for(int tmp = 0; tmp < alterFusionNum_1stEndShared; tmp++)
		{
			int tmpFusionIndex_1stEndShared 
				= breakPointInfoVec[tmpFusionIndex]->returnAlterFusionIndex_1stEndShared(tmp);
			string tmpFusionInfoStr_forAlterFusionDetection
				= breakPointInfoVec[tmpFusionIndex_1stEndShared]->returnFusionInfoStr_forAlterFusionDetection(indexInfo);
			tmpAlterFusionStr += tmpFusionInfoStr_forAlterFusionDetection;
			tmpAlterFusionStr += ",";
		}
		for(int tmp = 0; tmp < alterFusionNum_2ndEndShared; tmp++)
		{
			int tmpFusionIndex_2ndEndShared 
				= breakPointInfoVec[tmpFusionIndex]->returnAlterFusionIndex_2ndEndShared(tmp); 
			string tmpFusionInfoStr_forAlterFusionDetection
				= breakPointInfoVec[tmpFusionIndex_2ndEndShared]->returnFusionInfoStr_forAlterFusionDetection(indexInfo);
			tmpAlterFusionStr += tmpFusionInfoStr_forAlterFusionDetection;
			tmpAlterFusionStr += ",";
		}
		return tmpAlterFusionStr;
	}

	void outputFusionBreakPointHashInfoStr(string& outputFilePath, Index_Info* indexInfo)
	{
		ofstream breakPoint_ofs(outputFilePath.c_str());
		for(int tmp = 0; tmp < breakPointInfoVec.size(); tmp++)
		{
			string tmpBreakPointStr 
				= breakPointInfoVec[tmp]->returnFusionBreakPointStr(indexInfo);
			breakPoint_ofs << tmpBreakPointStr << endl;
		}
		breakPoint_ofs.close();
	}

	void outputFusionBreakPointHashInfoStr_withAlterFusionJunc(string& outputFilePath, Index_Info* indexInfo)
	{
		ofstream breakPoint_ofs(outputFilePath.c_str());
		for(int tmp = 0; tmp < breakPointInfoVec.size(); tmp++)
		{
			string tmpBreakPointStr 
				= breakPointInfoVec[tmp]->returnFusionBreakPointStr(indexInfo);
			string tmpAlterFusionStr = "";
			tmpAlterFusionStr = this->returnAlterFusionStr(tmp, indexInfo);
			if(tmpAlterFusionStr == "")
				tmpAlterFusionStr = "noAlterFusion\tkept_AF";
			else
				tmpAlterFusionStr += "\tfilterOut_AF";
			string tmpBreakPointStr_withAlterFusion = tmpBreakPointStr + "\t" + tmpAlterFusionStr;
			breakPoint_ofs << tmpBreakPointStr_withAlterFusion << endl;
		}
		breakPoint_ofs.close();
	}

	void insertNewFusionBreakPoint_toBreakPointMapVecVec_for(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, int toAddBreakPointIndexInVec)
	{
		First2secondBreakPointMap::iterator tmpBreakPointMapIter
			= ((breakPointMapVecVec_for[tmpChrNameInt_1])[tmpChrNameInt_2]).find(tmpBreakPoint_1);
		if(tmpBreakPointMapIter == ((breakPointMapVecVec_for[tmpChrNameInt_1])[tmpChrNameInt_2]).end()) // not found
		{
			SecondBreakPoint2FusionBreakPointHashIndexMap tmp2ndTypeMap;
			tmp2ndTypeMap.insert(pair<int,int>(tmpBreakPoint_2, toAddBreakPointIndexInVec));
			((breakPointMapVecVec_for[tmpChrNameInt_1])[tmpChrNameInt_2]).insert(
				pair<int,SecondBreakPoint2FusionBreakPointHashIndexMap>
					(tmpBreakPoint_1, tmp2ndTypeMap));
		}
		else
		{
			SecondBreakPoint2FusionBreakPointHashIndexMap::iterator tmp2ndTypeMapIter
				= (tmpBreakPointMapIter->second).find(tmpBreakPoint_2);
			if(tmp2ndTypeMapIter == (tmpBreakPointMapIter->second).end())
			{
				(tmpBreakPointMapIter->second).insert(
					pair<int,int>(tmpBreakPoint_2, toAddBreakPointIndexInVec));
			}
			else
			{
				cout << "error in insertNewFusionBreakPoint_toBreakPointMapVecVec_for" << endl;
				exit(1);
			}
		}
	}
	
	void insertNewFusionBreakPoint_toBreakPointMapVecVec_rev(
		int tmpChrNameInt_2, int tmpChrNameInt_1,
		int tmpBreakPoint_2, int tmpBreakPoint_1, int toAddBreakPointIndexInVec)
	{
		First2secondBreakPointMap::iterator tmpBreakPointMapIter
			= ((breakPointMapVecVec_rev[tmpChrNameInt_2])[tmpChrNameInt_1]).find(tmpBreakPoint_2);
		if(tmpBreakPointMapIter == ((breakPointMapVecVec_rev[tmpChrNameInt_2])[tmpChrNameInt_1]).end()) // not found
		{
			SecondBreakPoint2FusionBreakPointHashIndexMap tmp2ndTypeMap;
			tmp2ndTypeMap.insert(pair<int,int>(tmpBreakPoint_1, toAddBreakPointIndexInVec));
			((breakPointMapVecVec_rev[tmpChrNameInt_2])[tmpChrNameInt_1]).insert(
				pair<int,SecondBreakPoint2FusionBreakPointHashIndexMap>
					(tmpBreakPoint_2, tmp2ndTypeMap));
		}
		else
		{
			SecondBreakPoint2FusionBreakPointHashIndexMap::iterator tmp2ndTypeMapIter
				= (tmpBreakPointMapIter->second).find(tmpBreakPoint_1);
			if(tmp2ndTypeMapIter == (tmpBreakPointMapIter->second).end())
			{
				(tmpBreakPointMapIter->second).insert(
					pair<int,int>(tmpBreakPoint_1, toAddBreakPointIndexInVec));
			}
			else
			{
				cout << "error in insertNewFusionBreakPoint_toBreakPointMapVecVec_rev" << endl;
				exit(1);
			}
		}
	}

	int searchAndReturnIndexInBreakPointInfoVec_for(int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2)
	{
		First2secondBreakPointMap::iterator tmpBreakPointMapIter
			= ((breakPointMapVecVec_for[tmpChrNameInt_1])[tmpChrNameInt_2]).find(tmpBreakPoint_1);
		if(tmpBreakPointMapIter == ((breakPointMapVecVec_for[tmpChrNameInt_1])[tmpChrNameInt_2]).end()) // not found
			return -1;
		else
		{
			SecondBreakPoint2FusionBreakPointHashIndexMap::iterator tmpInfoVecIndexMapIter
				= (tmpBreakPointMapIter->second).find(tmpBreakPoint_2);
			if(tmpInfoVecIndexMapIter == (tmpBreakPointMapIter->second).end())
			{
				return -1;
			}
			else
			{
				return (tmpInfoVecIndexMapIter->second);
			}
		}
	}

	int searchAndReturnIndexInBreakPointInfoVec_rev(int tmpChrNameInt_2, int tmpChrNameInt_1,
		int tmpBreakPoint_2, int tmpBreakPoint_1)
	{
		First2secondBreakPointMap::iterator tmpBreakPointMapIter
			= ((breakPointMapVecVec_rev[tmpChrNameInt_2])[tmpChrNameInt_1]).find(tmpBreakPoint_2);
		if(tmpBreakPointMapIter == ((breakPointMapVecVec_rev[tmpChrNameInt_2])[tmpChrNameInt_1]).end()) // not found
			return -1;
		else
		{
			SecondBreakPoint2FusionBreakPointHashIndexMap::iterator tmpInfoVecIndexMapIter
				= (tmpBreakPointMapIter->second).find(tmpBreakPoint_1);
			if(tmpInfoVecIndexMapIter == (tmpBreakPointMapIter->second).end())
			{
				return -1;
			}
			else
			{
				return (tmpInfoVecIndexMapIter->second);
			}
		}
	}

	void supportNumIncrementWithIndex(int index)
	{
		return breakPointInfoVec[index]->addSupportNum(1);
	}

	void supportNumIncrementWithIndex_encompassing(int index)
	{
		return breakPointInfoVec[index]->addSupportNum_encompassing(1);
	}

	int returnBreakPointInfoVecSize()
	{
		return breakPointInfoVec.size();
	}

	int returnChrNameIntWithIndex_1(int index)
	{
		return breakPointInfoVec[index]->returnChrNameInt_1();
	}

	int returnChrNameIntWithIndex_2(int index)
	{
		return breakPointInfoVec[index]->returnChrNameInt_2();
	}

	string returnChrNameStrWithIndex_1(int index, Index_Info* indexInfo)
	{
		int tmpChrNameInt_1 = breakPointInfoVec[index]->returnChrNameInt_1();
		return indexInfo->returnChrNameStr(tmpChrNameInt_1);
	}

	string returnChrNameStrWithIndex_2(int index, Index_Info* indexInfo)
	{
		int tmpChrNameInt_2 = breakPointInfoVec[index]->returnChrNameInt_2();
		return indexInfo->returnChrNameStr(tmpChrNameInt_2);
	}

	int returnBreakPointWithIndex_1(int index)
	{
		return breakPointInfoVec[index]->returnBreakPoint_1();
	}

	int returnBreakPointWithIndex_2(int index)
	{
		return breakPointInfoVec[index]->returnBreakPoint_2();
	}

	int returnSupportNumWithIndex(int index)
	{
		return breakPointInfoVec[index]->returnSupportNum();
	}

	int returnSupportNumWithIndex_encompassing(int index)
	{
		return breakPointInfoVec[index]->returnSupportNum_encompassing();
	}

	string returnStrandWithIndex_1(int index)
	{
		return breakPointInfoVec[index]->returnStrand_1();
	}

	string returnStrandWithIndex_2(int index)
	{
		return breakPointInfoVec[index]->returnStrand_2();
	}

	int returnAnchorLengthWithIndex_1(int index)
	{
		return breakPointInfoVec[index]->returnAnchorLength_1();
	}

	int returnAnchorLengthWithIndex_2(int index)
	{
		return breakPointInfoVec[index]->returnAnchorLength_2();
	}

	int returnXMminWithIndex(int index)
	{
		return breakPointInfoVec[index]->returnXMmin();
	}	

	int returnXMmaxWithIndex(int index)
	{
		return breakPointInfoVec[index]->returnXMmax();
	}	

	void memoryFree()
	{
		for(int tmp = 0; tmp < breakPointInfoVec.size(); tmp++)
		{
			delete breakPointInfoVec[tmp];
		}
	}
};

#endif