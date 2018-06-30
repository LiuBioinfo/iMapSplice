// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef PAIRING_INFO_H
#define PAIRING_INFO_H

using namespace std;

class Pairing_Info
{
private:
	vector<int> subSegGroupVec_Nor1_regionIndex;
	vector<int> subSegGroupVec_Rcm1_regionIndex;
	vector<int> subSegGroupVec_Nor2_regionIndex;
	vector<int> subSegGroupVec_Rcm2_regionIndex;

	vector<int> subSegGroupVec_Nor1_pathIndex_inPathVec_seg;
	vector<int> subSegGroupVec_Rcm1_pathIndex_inPathVec_seg;
	vector<int> subSegGroupVec_Nor2_pathIndex_inPathVec_seg;
	vector<int> subSegGroupVec_Rcm2_pathIndex_inPathVec_seg;	
	//set<int> subSegGroupVec_indexSet_Nor1;
	map<int, int> subSegGroupVec_indexMap_Rcm1; // map<regionIndex, index_subSegGroupVec_Rcm1>
	//set<int> subSegGroupVec_indexSet_Nor2;
	map<int, int> subSegGroupVec_indexMap_Rcm2;

	vector<bool> subSegGroupVec_Nor1_candiPairedBool;
	vector<bool> subSegGroupVec_Rcm1_candiPairedBool;
	vector<bool> subSegGroupVec_Nor2_candiPairedBool;
	vector<bool> subSegGroupVec_Rcm2_candiPairedBool;

	vector< pair<int, int> > pairIndexVec_inPathVec_seg_Nor1Rcm2; // index in subSegGroupVec_Nor1
	vector< pair<int, int> > pairIndexVec_inPathVec_seg_Nor2Rcm1;	

	vector< pair<int, int> > finalPairIndexVec_inFixedPathVec_Nor1Rcm2; // index in subSegGroupVec_Nor1
	vector< pair<int, int> > finalPairIndexVec_inFixedPathVec_Nor2Rcm1;	
public:
	Pairing_Info()
	{

	}

	void pushBack_subSegGroupVec_regionIndex(int regionIndex, int typeNO)
	{
		if(typeNO == 1)
			subSegGroupVec_Nor1_regionIndex.push_back(regionIndex);
		else if(typeNO == 2)
			subSegGroupVec_Rcm1_regionIndex.push_back(regionIndex);
		else if(typeNO == 3)
			subSegGroupVec_Nor2_regionIndex.push_back(regionIndex);
		else if(typeNO == 4)
			subSegGroupVec_Rcm2_regionIndex.push_back(regionIndex);
		else
		{
			cout << "error in pushBack_subSegGroupVec_index" << endl;
			exit(1);
		}
	}

	void pushBack_index_inPathInfo_PathVec_seg(int pathIndex, int typeNO)
	{
		if(typeNO == 1)
			subSegGroupVec_Nor1_pathIndex_inPathVec_seg.push_back(pathIndex);
		else if(typeNO == 2)
			subSegGroupVec_Nor1_pathIndex_inPathVec_seg.push_back(pathIndex);
		else if(typeNO == 3)
			subSegGroupVec_Nor1_pathIndex_inPathVec_seg.push_back(pathIndex);
		else if(typeNO == 4)
			subSegGroupVec_Nor1_pathIndex_inPathVec_seg.push_back(pathIndex);
		else
		{
			cout << "error in pushBack_subSegGroupVec_index" << endl;
			exit(1);
		}
	}

	void pushBack_subSegGroupVec_bool(int typeNO)
	{
		if(typeNO == 1)
			subSegGroupVec_Nor1_candiPairedBool.push_back(false);
		else if(typeNO == 2)
			subSegGroupVec_Rcm1_candiPairedBool.push_back(false);
		else if(typeNO == 3)
			subSegGroupVec_Nor2_candiPairedBool.push_back(false);
		else if(typeNO == 4)
			subSegGroupVec_Rcm2_candiPairedBool.push_back(false);
		else
		{
			cout << "error in pushBack_subSegGroupVec_index" << endl;
			exit(1);
		}		
	}

	void generate_subSegGroupVec_indexSet()
	{	
		for(int tmp = 0; tmp < subSegGroupVec_Rcm1_regionIndex.size(); tmp++)
		{
			int tmpRegionIndex = subSegGroupVec_Rcm1_regionIndex[tmp];
			subSegGroupVec_indexMap_Rcm1.insert(pair<int,int>(tmpRegionIndex, tmp));
		}
		for(int tmp = 0; tmp < subSegGroupVec_Rcm2_regionIndex.size(); tmp++)
		{
			int tmpRegionIndex = subSegGroupVec_Rcm2_regionIndex[tmp];
			subSegGroupVec_indexMap_Rcm2.insert(pair<int,int>(tmpRegionIndex, tmp));
		}								
	}

	void pair_subSegGroup()
	{
		for(int tmp = 0; tmp < subSegGroupVec_Nor1_regionIndex.size(); tmp++)
		{
			int tmpRegionIndex_Nor1 = subSegGroupVec_Nor1_regionIndex[tmp];
			Map<int,int>::iterator tmpMapIter;
			tmpMapIter = subSegGroupVec_indexMap_Rcm2.find(tmpRegionIndex_Nor1);
			if(tmpMapIter != subSegGroupVec_indexMap_Rcm2.end())
			{
				int targetRcm2Index = tmpMapIter->second;
				subSegGroupVec_Nor1_candiPairedBool[tmp] = true;
				subSegGroupVec_Rcm2_candiPairedBool[targetRcm2Index] = true;
				pairIndexVec_inPathVec_seg_Nor1Rcm2.push_back(pair<int,int>(tmp, targetRcm2Index));
			}
			tmpMapIter = subSegGroupVec_indexMap_Rcm2.find(tmpRegionIndex_Nor1+1);
			if(tmpMapIter != subSegGroupVec_indexMap_Rcm2.end())
			{
				int targetRcm2Index = tmpMapIter->second;
				subSegGroupVec_Nor1_candiPairedBool[tmp] = true;
				subSegGroupVec_Rcm2_candiPairedBool[targetRcm2Index] = true;
				pairIndexVec_inPathVec_seg_Nor1Rcm2.push_back(pair<int,int>(tmp, targetRcm1Index));
			}
		}

		for(int tmp = 0; tmp < subSegGroupVec_Nor2_regionIndex.size(); tmp++)
		{
			int tmpRegionIndex_Nor2 = subSegGroupVec_Nor2_regionIndex[tmp];
			Map<int,int>::iterator tmpMapIter;
			tmpMapIter = subSegGroupVec_indexMap_Rcm1.find(tmpRegionIndex_Nor2);
			if(tmpMapIter != subSegGroupVec_indexMap_Rcm1.end())
			{
				int targetRcm1Index = tmpMapIter->second;
				subSegGroupVec_Nor2_candiPairedBool[tmp] = true;
				subSegGroupVec_Rcm1_candiPairedBool[targetRcm1Index] = true;
				pairIndexVec_inPathVec_seg_Nor2Rcm1.push_back(pair<int,int>(tmp, (tmpMapIter->second)));
			}
			tmpMapIter = subSegGroupVec_indexMap_Rcm1.find(tmpRegionIndex_Nor2+1);
			if(tmpMapIter != subSegGroupVec_indexMap_Rcm1.end())
			{
				int targetRcm1Index = tmpMapIter->second;
				subSegGroupVec_Nor2_candiPairedBool[tmp] = true;
				subSegGroupVec_Rcm1_candiPairedBool[targetRcm1Index] = true;
				pairIndexVec_inPathVec_seg_Nor2Rcm1.push_back(pair<int,int>(tmp, (tmpMapIter->second)));
			}
		}		
	}

	bool returnSubSegGroupVec_candiPairedBool(int typeNO, int index)
	{
		if(typeNO == 1)
		{
			return subSegGroupVec_Nor1_candiPairedBool[index];
		}
		else if(typeNO == 2)
		{
			return subSegGroupVec_Rcm1_candiPairedBool[index];			
		}
		else if(typeNO == 3)
		{
			return subSegGroupVec_Nor2_candiPairedBool[index];			
		}
		else if(typeNO == 4)
		{
			return subSegGroupVec_Rcm_candiPairedBool[index];			
		}
		else
		{
			cout << "error in typeNO" << endl;
			exit(1);
		}						
	}
};

#endif