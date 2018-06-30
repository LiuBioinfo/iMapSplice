// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALIGNINFERJUNCTIONHASH_INFO_H
#define ALIGNINFERJUNCTIONHASH_INFO_H

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

#include "alignInfer_info.h"

using namespace std;

typedef map<int, int> AcceptorStartPos2AlignInferInfoMap;
typedef map<int, AcceptorStartPos2AlignInferInfoMap> AlignInferJunctionMap; 
// (<donerEnd, acceptorStart>, index_SJmapIndexVec) or (<acceptorStart, donerEnd>, index_SJmapIndexVec)

class AlignInferJunctionHash_Info
{
private:
	vector<AlignInferJunctionMap> alignInferJunctionMapVec;
	//vector<AlignInferJunctionMap> alignIfnerJunctionMapVec_rev;

	int currentAlignInferInfoVecIndex;
public:
	vector<AlignInfer_Info> alignInferInfoVec;

	AlignInferJunctionHash_Info()
	{
		currentAlignInferInfoVecIndex = 0;
	}

	void generateSpliceSiteSetVec(vector< set<int> >& spliceSiteSetVec_doner, 
		vector< set<int> >& spliceSiteSetVec_acceptor)
	{
		int tmpAlignInferInfoVecSize = alignInferInfoVec.size();
		for(int tmp = 0; tmp < tmpAlignInferInfoVecSize; tmp++)
		{	
			int tmpAlignInferInfo_chrNameInt = alignInferInfoVec[tmp].returnChrNameInt();
			int tmpAlignInferInfo_donerEndPosInt = alignInferInfoVec[tmp].returnDonerEndPos();
			int tmpAlignInferInfo_acceptorStartPosInt = alignInferInfoVec[tmp].returnAcceptorStartPos();
			spliceSiteSetVec_doner[tmpAlignInferInfo_chrNameInt].insert(tmpAlignInferInfo_donerEndPosInt);
			spliceSiteSetVec_acceptor[tmpAlignInferInfo_chrNameInt].insert(tmpAlignInferInfo_acceptorStartPosInt);
		}
	}

	string returnAlignInferInfo_flankString(int tmp, Index_Info* indexInfo)
	{
		return alignInferInfoVec[tmp].returnFlankString(indexInfo);
	}

	string returnAlignInferInfo_flankStringChanged(int tmp)
	{
		return alignInferInfoVec[tmp].returnFlankStringChanged();
	}

	int returnAnchorSizeMax_doner(int tmp)
	{
		return alignInferInfoVec[tmp].returnAnchorSizeMax_doner();
	}

	int returnAnchorSizeMax_acceptor(int tmp)
	{
		return alignInferInfoVec[tmp].returnAnchorSizeMax_acceptor();
	}

	int returnAlignInferInfo_chrNameInt(int index)
	{
		return alignInferInfoVec[index].returnChrNameInt();
	}

	int returnAlignInferInfo_donerAnchorSizeMax(int index)
	{
		return alignInferInfoVec[index].returnAnchorSizeMax_doner();
	}

	int returnAlignInferInfo_acceptorAnchorSizeMax(int index)
	{
		return alignInferInfoVec[index].returnAnchorSizeMax_acceptor();
	}

	int returnAlignInferInfo_XMmin(int index)
	{
		return alignInferInfoVec[index].returnXMmin();
	}

	int returnAlignInferInfo_XMmax(int index)
	{
		return alignInferInfoVec[index].returnXMmax();
	}

	int returnAlignInferInfo_donerEndPos(int index)
	{
		return alignInferInfoVec[index].returnDonerEndPos();
	}

	int returnAlignInferInfo_acceptorStartPos(int index)
	{
		return alignInferInfoVec[index].returnAcceptorStartPos();
	}

	int returnAlignInferInfo_supportNum(int index)
	{
		return alignInferInfoVec[index].returnSupportNum();
	}

	int returnAlignInferInfoVecSize()
	{
		return alignInferInfoVec.size();
	}

	// void countJuncNum_total_validIntronSize_invalidIntronSize(
	// 	int min_intron_size, int max_intron_size,
	// 	int& junc_total_num, 
	// 	int& junc_validIntronSize_num, int& junc_invalidIntronSize_num
	// 	// int& donerSpliceSite_total_num,
	// 	// int& donerSpliceSite_validIntronSize_num, int& donerSpliceSite_invalidIntronSize_num,
	// 	// int& acceptorSpliceSite_total_num,
	// 	// int& acceptorSpliceSite_validIntronSize_num, int& acceptorSpliceSite_invalidIntronSize_num
	// 	)
	// {
	// 	int alignInferJuncInfoVecSize = alignInferInfoVec.size();
	// 	int validJuncNum = 0;
	// 	int invalidJuncNum = 0;
	// 	for(int tmp = 0; tmp < alignInferJuncInfoVecSize; tmp++)
	// 	{
	// 		int tmpJunc_donerSpliceSite = this->returnAlignInferInfo_donerEndPos();
	// 		int tmpJunc_acceptorSpliceSite = this->returnAlignInferInfo_acceptorStartPos();
	// 		int tmpJunc_intronSize = tmpJunc_acceptorSpliceSite - tmpJunc_donerSpliceSite - 1;
	// 		if((tmpJunc_intronSize >= min_intron_size)&&(tmpJunc_intronSize <= max_intron_size))
	// 			validJuncNum ++;
	// 		else
	// 			invalidJuncNum ++;
	// 	}
	// 	junc_total_num = alignInferJuncInfoVecSize;
	// 	junc_validIntronSize_num = validJuncNum;
	// 	junc_invalidIntronSize_num = invalidJuncNum;
	// }

	int searchAndReturnAlignInferJuncHashSupNum(int tmpChrInt,
		int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		int tmpIndex = this->searchAndReturnAlignInferInfoVecIndex(
			tmpChrInt, tmpDonerEndPos, tmpAcceptorStartPos);
		if(tmpIndex < 0)
			return 0;
		else
			return this->returnAlignInferInfo_supportNum(tmpIndex);
	}

	void compared2otherAlignInferJuncHashVecAndReturnSJsupNumVec(
		vector< AlignInferJunctionHash_Info* >& juncHashVec,
		vector<int>& chrNameIntVecInThisJuncHash,
		vector<int>& donerSpliceSiteVecInThisJuncHash,
		vector<int>& acceptorSpliceSiteVecInThisJuncHash,
		vector<int>& supNumInThisJuncHash,
		vector< vector<int> >& supNumVecInJuncHashVec2compare,
		vector<bool>& sharedByAllOrNotBoolVec,
		int min_intron_size, int max_intron_size,
		int& junc_total_num, int& junc_valid_num, int& junc_invalid_num,
		Index_Info* indexInfo)
	{
		//cout << "compared2otherAlignInferJuncHashVecAndReturnSJsupNumVec starts ......" << endl;
		int tmpValidJuncNum = 0;
		int tmpInvalidJuncNum = 0;
		int alignInferJuncInfoVecSize = alignInferInfoVec.size();
		for(int tmp = 0; tmp < alignInferJuncInfoVecSize; tmp++)
		{
			//cout << "tmp in alignInferJuncInfoVec: " << tmp << endl;
			int tmpJunc_chrNameInt = this->returnAlignInferInfo_chrNameInt(tmp);
			string tmpJunc_chrNameStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
			int tmpJunc_donerSpliceSite = this->returnAlignInferInfo_donerEndPos(tmp);
			int tmpJunc_acceptorSpliceSite = this->returnAlignInferInfo_acceptorStartPos(tmp);
			int tmpJunc_supportNum = this->returnAlignInferInfo_supportNum(tmp);
			int tmpJunc_intronSize = tmpJunc_acceptorSpliceSite - tmpJunc_donerSpliceSite - 1;
			//cout << "tmpJunc_chrNameStr: " << tmpJunc_chrNameStr << endl;
			//cout << "tmpJunc_donerSpliceSite: " << tmpJunc_donerSpliceSite << endl;
			//cout << "tmpJunc_acceptorSpliceSite: " << tmpJunc_acceptorSpliceSite << endl;
			//cout << "tmpJunc_supportNum: " << tmpJunc_supportNum << endl;
			if((tmpJunc_intronSize >= min_intron_size)&&(tmpJunc_intronSize <= max_intron_size))
			{
				//cout << "junction valid ..." << endl;
				tmpValidJuncNum ++;
				chrNameIntVecInThisJuncHash.push_back(tmpJunc_chrNameInt);
				donerSpliceSiteVecInThisJuncHash.push_back(tmpJunc_donerSpliceSite);
				acceptorSpliceSiteVecInThisJuncHash.push_back(tmpJunc_acceptorSpliceSite);
				supNumInThisJuncHash.push_back(tmpJunc_supportNum);
				bool tmpSharedByAllSJfileBool = true;
				vector<int> tmp2compareJuncSupNumVec;
				//cout << "juncHashVec.size(): " << juncHashVec.size() << endl;
				for(int tmpAlignInferJuncHash = 0; tmpAlignInferJuncHash < juncHashVec.size(); 
					tmpAlignInferJuncHash ++)
				{
					int tmp2compareJuncSupNum 
						= juncHashVec[tmpAlignInferJuncHash]->searchAndReturnAlignInferJuncHashSupNum(
							tmpJunc_chrNameInt, tmpJunc_donerSpliceSite, tmpJunc_acceptorSpliceSite);
					tmp2compareJuncSupNumVec.push_back(tmp2compareJuncSupNum);
					if(tmp2compareJuncSupNum == 0)
						tmpSharedByAllSJfileBool = false;
				}
				supNumVecInJuncHashVec2compare.push_back(tmp2compareJuncSupNumVec);
				if(tmpSharedByAllSJfileBool)
					sharedByAllOrNotBoolVec.push_back(true);
				else
					sharedByAllOrNotBoolVec.push_back(false);
			}
			else
				tmpInvalidJuncNum ++;
		}

		junc_total_num = alignInferInfoVec.size();
		junc_valid_num = tmpValidJuncNum;
		junc_invalid_num = tmpInvalidJuncNum;
	}

	void countTotalValidInvalidJuncNum_compareWithAnotherAlignInferJuncHash_juncWise(
		AlignInferJunctionHash_Info* juncHash_compared2, int offset,
		int min_intron_size, int max_intron_size,
		int& junc_total_num, int& junc_valid_num, int& junc_invalid_num,
		int& foundValidJuncNum_inThisJuncHash_withinOffset, 
		int& unfoundValidJuncNum_inThisJuncHash_withinOffset,
		string& foundSJ_inJuncFile, string& unfoundSJ_inJuncFile,
		Index_Info* indexInfo)
	{
		ofstream foundSJ_ofs(foundSJ_inJuncFile.c_str());
		ofstream unfoundSJ_ofs(unfoundSJ_inJuncFile.c_str());
		int alignInferJuncInfoVecSize = alignInferInfoVec.size();
		int validJuncNum = 0;
		int invalidJuncNum = 0;
		int validJuncNum_found = 0;
		int validJuncNum_unfound = 0;
		for(int tmp = 0; tmp < alignInferJuncInfoVecSize; tmp++)
		{
			int tmpJunc_chrNameInt = this->returnAlignInferInfo_chrNameInt(tmp);
			string tmpJunc_chrNameStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
			int tmpJunc_donerSpliceSite = this->returnAlignInferInfo_donerEndPos(tmp);
			int tmpJunc_acceptorSpliceSite = this->returnAlignInferInfo_acceptorStartPos(tmp);
			int tmpJunc_supportNum = this->returnAlignInferInfo_supportNum(tmp);
			int tmpJunc_intronSize = tmpJunc_acceptorSpliceSite - tmpJunc_donerSpliceSite - 1;
			if((tmpJunc_intronSize >= min_intron_size)&&(tmpJunc_intronSize <= max_intron_size))
			{
				validJuncNum ++;
				bool foundInJuncHashCompared2 
					= juncHash_compared2->foundInAlignInferJunctionHash(
						tmpJunc_chrNameInt, tmpJunc_donerSpliceSite, 
						tmpJunc_acceptorSpliceSite, offset);
				if(foundInJuncHashCompared2)
				{
					validJuncNum_found ++;
					foundSJ_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_donerSpliceSite
						<< "\t" << tmpJunc_acceptorSpliceSite << "\tJUNC_ID\t" << tmpJunc_supportNum << endl;
				}
				else
				{
					validJuncNum_unfound ++;
					unfoundSJ_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_donerSpliceSite
						<< "\t" << tmpJunc_acceptorSpliceSite  << "\tJUNC_ID\t" << tmpJunc_supportNum << endl;
				}
			}
			else
				invalidJuncNum ++;
		}

		junc_total_num = alignInferJuncInfoVecSize;
		junc_valid_num = validJuncNum;
		junc_invalid_num = invalidJuncNum;
		foundValidJuncNum_inThisJuncHash_withinOffset 
			= validJuncNum_found;
		unfoundValidJuncNum_inThisJuncHash_withinOffset
			= validJuncNum_unfound;
		foundSJ_ofs.close();
		unfoundSJ_ofs.close();
	}

	void countTotalValidInvalidJuncNum_compareWithAnotherAlignInferJuncHash_juncWise_outputChrNamePosOnly(
		AlignInferJunctionHash_Info* juncHash_compared2, int offset,
		int min_intron_size, int max_intron_size,
		int& junc_total_num, int& junc_valid_num, int& junc_invalid_num,
		int& foundValidJuncNum_inThisJuncHash_withinOffset, 
		int& unfoundValidJuncNum_inThisJuncHash_withinOffset,
		string& foundSJ_inJuncFile, string& unfoundSJ_inJuncFile,
		Index_Info* indexInfo)
	{
		ofstream foundSJ_ofs(foundSJ_inJuncFile.c_str());
		ofstream unfoundSJ_ofs(unfoundSJ_inJuncFile.c_str());
		int alignInferJuncInfoVecSize = alignInferInfoVec.size();
		int validJuncNum = 0;
		int invalidJuncNum = 0;
		int validJuncNum_found = 0;
		int validJuncNum_unfound = 0;
		for(int tmp = 0; tmp < alignInferJuncInfoVecSize; tmp++)
		{
			int tmpJunc_chrNameInt = this->returnAlignInferInfo_chrNameInt(tmp);
			string tmpJunc_chrNameStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
			int tmpJunc_donerSpliceSite = this->returnAlignInferInfo_donerEndPos(tmp);
			int tmpJunc_acceptorSpliceSite = this->returnAlignInferInfo_acceptorStartPos(tmp);
			//int tmpJunc_supportNum = this->returnAlignInferInfo_supportNum(tmp);
			int tmpJunc_intronSize = tmpJunc_acceptorSpliceSite - tmpJunc_donerSpliceSite - 1;
			if((tmpJunc_intronSize >= min_intron_size)&&(tmpJunc_intronSize <= max_intron_size))
			{
				validJuncNum ++;
				bool foundInJuncHashCompared2 
					= juncHash_compared2->foundInAlignInferJunctionHash(
						tmpJunc_chrNameInt, tmpJunc_donerSpliceSite, 
						tmpJunc_acceptorSpliceSite, offset);
				if(foundInJuncHashCompared2)
				{
					validJuncNum_found ++;
					foundSJ_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_donerSpliceSite
						<< "\t" << tmpJunc_acceptorSpliceSite << "\tJUNC_ID"  << endl;
				}
				else
				{
					validJuncNum_unfound ++;
					unfoundSJ_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_donerSpliceSite
						<< "\t" << tmpJunc_acceptorSpliceSite  << "\tJUNC_ID" << endl;
				}
			}
			else
				invalidJuncNum ++;
		}

		junc_total_num = alignInferJuncInfoVecSize;
		junc_valid_num = validJuncNum;
		junc_invalid_num = invalidJuncNum;
		foundValidJuncNum_inThisJuncHash_withinOffset 
			= validJuncNum_found;
		unfoundValidJuncNum_inThisJuncHash_withinOffset
			= validJuncNum_unfound;
		foundSJ_ofs.close();
		unfoundSJ_ofs.close();
	}


	void countTotalValidInvalidJuncNum_compareWithAnotherAlignInferJuncHash_juncWise_outputSupNumAnchorSize(
		AlignInferJunctionHash_Info* juncHash_compared2, int offset,
		int min_intron_size, int max_intron_size,
		int& junc_total_num, int& junc_valid_num, int& junc_invalid_num,
		int& foundValidJuncNum_inThisJuncHash_withinOffset, 
		int& unfoundValidJuncNum_inThisJuncHash_withinOffset,
		string& foundSJ_inJuncFile, string& unfoundSJ_inJuncFile,
		Index_Info* indexInfo)
	{
		ofstream foundSJ_ofs(foundSJ_inJuncFile.c_str());
		ofstream unfoundSJ_ofs(unfoundSJ_inJuncFile.c_str());
		int alignInferJuncInfoVecSize = alignInferInfoVec.size();
		int validJuncNum = 0;
		int invalidJuncNum = 0;
		int validJuncNum_found = 0;
		int validJuncNum_unfound = 0;
		for(int tmp = 0; tmp < alignInferJuncInfoVecSize; tmp++)
		{
			int tmpJunc_chrNameInt = this->returnAlignInferInfo_chrNameInt(tmp);
			string tmpJunc_chrNameStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
			int tmpJunc_donerSpliceSite = this->returnAlignInferInfo_donerEndPos(tmp);
			int tmpJunc_acceptorSpliceSite = this->returnAlignInferInfo_acceptorStartPos(tmp);
			int tmpJunc_supportNum = this->returnAlignInferInfo_supportNum(tmp);
			int tmpJunc_donerAnchorSize = this->returnAlignInferInfo_donerAnchorSizeMax(tmp);
			int tmpJunc_acceptorAnchorSize = this->returnAlignInferInfo_acceptorAnchorSizeMax(tmp);
			int tmpJunc_intronSize = tmpJunc_acceptorSpliceSite - tmpJunc_donerSpliceSite - 1;
			if((tmpJunc_intronSize >= min_intron_size)&&(tmpJunc_intronSize <= max_intron_size))
			{
				validJuncNum ++;
				bool foundInJuncHashCompared2 
					= juncHash_compared2->foundInAlignInferJunctionHash(
						tmpJunc_chrNameInt, tmpJunc_donerSpliceSite, 
						tmpJunc_acceptorSpliceSite, offset);
				if(foundInJuncHashCompared2)
				{
					validJuncNum_found ++;
					foundSJ_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_donerSpliceSite
						<< "\t" << tmpJunc_acceptorSpliceSite << "\tJUNC_ID\t" << tmpJunc_supportNum << "\t" 
						<< tmpJunc_donerAnchorSize << "\t" << tmpJunc_acceptorAnchorSize << endl;
				}
				else
				{
					validJuncNum_unfound ++;
					unfoundSJ_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_donerSpliceSite
						<< "\t" << tmpJunc_acceptorSpliceSite  << "\tJUNC_ID\t" << tmpJunc_supportNum << "\t" 
						<< tmpJunc_donerAnchorSize << "\t" << tmpJunc_acceptorAnchorSize << endl;
				}
			}
			else
				invalidJuncNum ++;
		}

		junc_total_num = alignInferJuncInfoVecSize;
		junc_valid_num = validJuncNum;
		junc_invalid_num = invalidJuncNum;
		foundValidJuncNum_inThisJuncHash_withinOffset 
			= validJuncNum_found;
		unfoundValidJuncNum_inThisJuncHash_withinOffset
			= validJuncNum_unfound;
		foundSJ_ofs.close();
		unfoundSJ_ofs.close();
	}

	int searchAndReturnJuncSupNumMax_withOffset(int tmpChrInt,
		int tmpDonerEndPos, int tmpAcceptorStartPos, int offset)
	{
		int tmpJuncSupNumMax = 0;
		for(int tmp = 0-offset; tmp <= offset; tmp++)
		{
			for(int tmp2 = 0-offset; tmp2 <= offset; tmp2++)
			{
				int tmpNewDonerEndPos = tmpDonerEndPos + tmp;
				int tmpNewAcceptorStartPos = tmpAcceptorStartPos + tmp2;
				int foundIndex = this->searchAndReturnAlignInferInfoVecIndex(
					tmpChrInt, tmpNewDonerEndPos, tmpNewAcceptorStartPos);
				if(foundIndex != -1)
				{
					int tmpJuncSupNum = alignInferInfoVec[foundIndex].returnSupportNum();
					if(tmpJuncSupNum >= tmpJuncSupNumMax)
						tmpJuncSupNumMax = tmpJuncSupNum;
				}
			}
		}
		return tmpJuncSupNumMax;
	}

	void compareWithAnotherAlignInferJuncHash_juncWise_varySupNumMaxInTheOtherJuncHash(
		AlignInferJunctionHash_Info* theOtherJuncHash, int offset, 
		int min_intron_size_int, int max_intron_size_int, int sup_num_thres_max,
		vector<int>& foundJuncNumVec_inThisJuncHash_withinOffset_varySupNum, 
		Index_Info* indexInfo)
	{
		int alignInferJuncInfoVecSize = alignInferInfoVec.size();
		for(int tmp = 0; tmp < alignInferJuncInfoVecSize; tmp++)
		{
			//cout << "tmp: " << tmp+1 << endl;
			int tmpJunc_chrNameInt = this->returnAlignInferInfo_chrNameInt(tmp);
			//string tmpJunc_chrNameStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
			int tmpJunc_donerSpliceSite = this->returnAlignInferInfo_donerEndPos(tmp);
			int tmpJunc_acceptorSpliceSite = this->returnAlignInferInfo_acceptorStartPos(tmp);
			int tmpJunc_intronSize = tmpJunc_acceptorSpliceSite - tmpJunc_donerSpliceSite - 1;
			if((tmpJunc_intronSize >= min_intron_size_int)&&(tmpJunc_intronSize <= max_intron_size_int))
			{
				int tmpJuncSupNumMax_withOffset_inTheOtherJuncHash
					= theOtherJuncHash->searchAndReturnJuncSupNumMax_withOffset(
						tmpJunc_chrNameInt, tmpJunc_donerSpliceSite, 
						tmpJunc_acceptorSpliceSite, offset);
				if(tmpJuncSupNumMax_withOffset_inTheOtherJuncHash >= sup_num_thres_max)
					tmpJuncSupNumMax_withOffset_inTheOtherJuncHash = sup_num_thres_max;
				if(tmpJuncSupNumMax_withOffset_inTheOtherJuncHash >= 1)
				{
					for(int tmpSupportNum = 1; tmpSupportNum <= tmpJuncSupNumMax_withOffset_inTheOtherJuncHash; tmpSupportNum++)
						foundJuncNumVec_inThisJuncHash_withinOffset_varySupNum[tmpSupportNum-1]++;
				}
				else
				{}
			}
			else
			{}
		}
	}

	void mergeWithAnotherAlignInferJuncHash_chrNamePosOnly(
		AlignInferJunctionHash_Info* tmpAlignInferJuncHashInfo, Index_Info* indexInfo)
	{
		int tmpAlignInferInfoVecSize = tmpAlignInferJuncHashInfo->returnAlignInferInfoVecSize();
		for(int tmp = 0; tmp < tmpAlignInferInfoVecSize; tmp++)
		{
			int mapChrNameInt = tmpAlignInferJuncHashInfo->returnAlignInferInfo_chrNameInt(tmp);
			int tmpSJposDonerEnd = tmpAlignInferJuncHashInfo->returnAlignInferInfo_donerEndPos(tmp);
			int tmpSJposAcceptorStart = tmpAlignInferJuncHashInfo->returnAlignInferInfo_acceptorStartPos(tmp);
			
			AlignInferJunctionMap::iterator alignInferJuncMapIter 
				= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
			if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
			{
				//cout << "not found" << endl;
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePosOnly(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, indexInfo);
				//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
				AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
				tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

				alignInferJunctionMapVec[mapChrNameInt].insert(
					pair< int, AcceptorStartPos2AlignInferInfoMap> (
						tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else // old donerEndPos found, search for acceptorStartPos
			{
				AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
					= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
				//cout << "found" << endl;
				if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
				{
					AlignInfer_Info tmpAlignInferInfo;
					tmpAlignInferInfo.initiateAlignInferInfo_chrNamePosOnly(
						mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, indexInfo);				
					(alignInferJuncMapIter->second).insert(pair<int, int>(
						tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
					alignInferInfoVec.push_back(tmpAlignInferInfo);
					currentAlignInferInfoVecIndex ++;
				}
				else
				{	
					// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
					// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
					// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
					//alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
				}
			}			
		}
	}

	void mergeWithAnotherAlignInferJuncHash_chrNamePos_supportNum(
		AlignInferJunctionHash_Info* tmpAlignInferJuncHashInfo, Index_Info* indexInfo)
	{
		int tmpAlignInferInfoVecSize = tmpAlignInferJuncHashInfo->returnAlignInferInfoVecSize();
		for(int tmp = 0; tmp < tmpAlignInferInfoVecSize; tmp++)
		{
			int mapChrNameInt = tmpAlignInferJuncHashInfo->returnAlignInferInfo_chrNameInt(tmp);
			int tmpSJposDonerEnd = tmpAlignInferJuncHashInfo->returnAlignInferInfo_donerEndPos(tmp);
			int tmpSJposAcceptorStart = tmpAlignInferJuncHashInfo->returnAlignInferInfo_acceptorStartPos(tmp);
			int tmpSJsupportNum = tmpAlignInferJuncHashInfo->returnAlignInferInfo_supportNum(tmp);

			AlignInferJunctionMap::iterator alignInferJuncMapIter 
				= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
			if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
			{
				//cout << "not found" << endl;
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
					tmpSJsupportNum, indexInfo);
				//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
				AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
				tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

				alignInferJunctionMapVec[mapChrNameInt].insert(
					pair< int, AcceptorStartPos2AlignInferInfoMap> (
						tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else // old donerEndPos found, search for acceptorStartPos
			{
				AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
					= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
				//cout << "found" << endl;
				if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
				{
					AlignInfer_Info tmpAlignInferInfo;
					tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum(
						mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
						tmpSJsupportNum, indexInfo);				
					(alignInferJuncMapIter->second).insert(pair<int, int>(
						tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
					alignInferInfoVec.push_back(tmpAlignInferInfo);
					currentAlignInferInfoVecIndex ++;
				}
				else
				{	
					// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
					// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
					// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
					alignInferInfoVec[acceptorStartPosMapIter->second].addSupportNum(
						tmpSJsupportNum);
				}
			}			
		}
	}

	void mergeWithAnotherAlignInferJuncHash_chrNamePos_supportNum_anchorSize(
		AlignInferJunctionHash_Info* tmpAlignInferJuncHashInfo, Index_Info* indexInfo)
	{
		int tmpAlignInferInfoVecSize = tmpAlignInferJuncHashInfo->returnAlignInferInfoVecSize();
		for(int tmp = 0; tmp < tmpAlignInferInfoVecSize; tmp++)
		{
			int mapChrNameInt = tmpAlignInferJuncHashInfo->returnAlignInferInfo_chrNameInt(tmp);
			int tmpSJposDonerEnd = tmpAlignInferJuncHashInfo->returnAlignInferInfo_donerEndPos(tmp);
			int tmpSJposAcceptorStart = tmpAlignInferJuncHashInfo->returnAlignInferInfo_acceptorStartPos(tmp);
			int tmpSJsupportNum = tmpAlignInferJuncHashInfo->returnAlignInferInfo_supportNum(tmp);
			int tmpSJdonerAnchorSizeMax = tmpAlignInferJuncHashInfo->returnAlignInferInfo_donerAnchorSizeMax(tmp);
			int tmpSJacceptorAnchorSizeMax = tmpAlignInferJuncHashInfo->returnAlignInferInfo_acceptorAnchorSizeMax(tmp);

			AlignInferJunctionMap::iterator alignInferJuncMapIter 
				= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
			if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
			{
				//cout << "not found" << endl;
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
					tmpSJsupportNum, tmpSJdonerAnchorSizeMax, tmpSJacceptorAnchorSizeMax, indexInfo);
				//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
				AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
				tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

				alignInferJunctionMapVec[mapChrNameInt].insert(
					pair< int, AcceptorStartPos2AlignInferInfoMap> (
						tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else // old donerEndPos found, search for acceptorStartPos
			{
				AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
					= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
				//cout << "found" << endl;
				if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
				{
					AlignInfer_Info tmpAlignInferInfo;
					tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize(
						mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
						tmpSJsupportNum, tmpSJdonerAnchorSizeMax, tmpSJacceptorAnchorSizeMax, indexInfo);				
					(alignInferJuncMapIter->second).insert(pair<int, int>(
						tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
					alignInferInfoVec.push_back(tmpAlignInferInfo);
					currentAlignInferInfoVecIndex ++;
				}
				else
				{	
					// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
					// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
					// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
					alignInferInfoVec[acceptorStartPosMapIter->second].addSupportNum(
						tmpSJsupportNum);
					alignInferInfoVec[acceptorStartPosMapIter->second].updateAnchorSize(
						tmpSJdonerAnchorSizeMax, tmpSJacceptorAnchorSizeMax);
				}
			}			
		}
	}	

	void mergeWithAnotherAlignInferJuncHash_chrNamePos_supportNum_anchorSize_XM(
		AlignInferJunctionHash_Info* tmpAlignInferJuncHashInfo, Index_Info* indexInfo)
	{
		int tmpAlignInferInfoVecSize = tmpAlignInferJuncHashInfo->returnAlignInferInfoVecSize();
		for(int tmp = 0; tmp < tmpAlignInferInfoVecSize; tmp++)
		{
			int mapChrNameInt = tmpAlignInferJuncHashInfo->returnAlignInferInfo_chrNameInt(tmp);
			int tmpSJposDonerEnd = tmpAlignInferJuncHashInfo->returnAlignInferInfo_donerEndPos(tmp);
			int tmpSJposAcceptorStart = tmpAlignInferJuncHashInfo->returnAlignInferInfo_acceptorStartPos(tmp);
			int tmpSJsupportNum = tmpAlignInferJuncHashInfo->returnAlignInferInfo_supportNum(tmp);
			int tmpSJdonerAnchorSizeMax = tmpAlignInferJuncHashInfo->returnAlignInferInfo_donerAnchorSizeMax(tmp);
			int tmpSJacceptorAnchorSizeMax = tmpAlignInferJuncHashInfo->returnAlignInferInfo_acceptorAnchorSizeMax(tmp);
			int tmpSJminXM = tmpAlignInferJuncHashInfo->returnAlignInferInfo_XMmin(tmp);
			int tmpSJmaxXM = tmpAlignInferJuncHashInfo->returnAlignInferInfo_XMmax(tmp);
			AlignInferJunctionMap::iterator alignInferJuncMapIter 
				= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
			if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
			{
				//cout << "not found" << endl;
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize_XM(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
					tmpSJsupportNum, tmpSJdonerAnchorSizeMax, tmpSJacceptorAnchorSizeMax, 
					tmpSJminXM, tmpSJmaxXM, indexInfo);
				//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
				AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
				tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

				alignInferJunctionMapVec[mapChrNameInt].insert(
					pair< int, AcceptorStartPos2AlignInferInfoMap> (
						tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else // old donerEndPos found, search for acceptorStartPos
			{
				AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
					= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
				//cout << "found" << endl;
				if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
				{
					AlignInfer_Info tmpAlignInferInfo;
					tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize_XM(
						mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
						tmpSJsupportNum, tmpSJdonerAnchorSizeMax, tmpSJacceptorAnchorSizeMax, 
						tmpSJminXM, tmpSJmaxXM, indexInfo);				
					(alignInferJuncMapIter->second).insert(pair<int, int>(
						tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
					alignInferInfoVec.push_back(tmpAlignInferInfo);
					currentAlignInferInfoVecIndex ++;
				}
				else
				{	
					// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
					// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
					// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
					alignInferInfoVec[acceptorStartPosMapIter->second].addSupportNum(
						tmpSJsupportNum);
					alignInferInfoVec[acceptorStartPosMapIter->second].updateAnchorSize(
						tmpSJdonerAnchorSizeMax, tmpSJacceptorAnchorSizeMax);
					alignInferInfoVec[acceptorStartPosMapIter->second].updateXM(
						tmpSJminXM, tmpSJmaxXM);
				}
			}			
		}
	}	

	void refineMultiAlignment_SE(
		string& multi_sam_file_ori,
		string& multi_final_sam_file,
		string& multi_corrected_sam_file_original,
		string& multi_corrected_sam_file, Index_Info* indexInfo)
	{
		ifstream multi_sam_ori_ifs(multi_sam_file_ori.c_str());
		ofstream multi_final_sam_ofs(multi_final_sam_file.c_str());
		ofstream multi_corrected_sam_ori_ofs(multi_corrected_sam_file_original.c_str());
		ofstream multi_corrected_sam_ofs(multi_corrected_sam_file.c_str());
		while(!multi_sam_ori_ifs.eof())
		{
			string tmpLineStr;
			getline(multi_sam_ori_ifs, tmpLineStr);
			if((multi_sam_ori_ifs.eof())||(tmpLineStr == ""))
				break;
			multi_final_sam_ofs << tmpLineStr << endl;
		}
		multi_sam_ori_ifs.close();
		multi_final_sam_ofs.close();
		multi_corrected_sam_ori_ofs.close();
		multi_corrected_sam_ofs.close();
	}

	void refineMultiAlignment_PE(
		string& multi_sam_file_ori,
		string& multi_final_sam_file,
		string& multi_corrected_sam_file_original,
		string& multi_corrected_sam_file, Index_Info* indexInfo)
	{
		ifstream multi_sam_ori_ifs(multi_sam_file_ori.c_str());
		ofstream multi_final_sam_ofs(multi_final_sam_file.c_str());
		ofstream multi_corrected_sam_ori_ofs(multi_corrected_sam_file_original.c_str());
		ofstream multi_corrected_sam_ofs(multi_corrected_sam_file.c_str());
		while(!multi_sam_ori_ifs.eof())
		{
			string tmpLineStr;
			getline(multi_sam_ori_ifs, tmpLineStr);
			if((multi_sam_ori_ifs.eof())||(tmpLineStr == ""))
				break;
			multi_final_sam_ofs << tmpLineStr << endl;
		}
		multi_sam_ori_ifs.close();
		multi_final_sam_ofs.close();
		multi_corrected_sam_ori_ofs.close();
		multi_corrected_sam_ofs.close();
	}

	void correct_output_SAM_withLowSupportLowConfidenceSJ_SE(RefineSAM_SE_Vec_Info* refineSAMvecInfo_SE,
		ofstream& unique_final_sam_ofs, 
		ofstream& unique_corrected_sam_original_ofs, ofstream& unique_corrected_sam_ofs, Index_Info* indexInfo)
	{
		for(int tmpRefineSAMindex = 0; tmpRefineSAMindex < refineSAMvecInfo_SE->returnRefineSAMvecSize(); tmpRefineSAMindex++)
		{
			if(refineSAMvecInfo_SE->point2NULLbool(tmpRefineSAMindex))
				continue;
			bool lowSupportLowConfidenceSJ_exists_bool = false;
			for(int tmpSJindex = 0; tmpSJindex < refineSAMvecInfo_SE->returnRefineSAMinfoSJvecSize(tmpRefineSAMindex); tmpSJindex++)
			{
				int tmpSJindexInAlignInferJuncHash = refineSAMvecInfo_SE->returnIndexInAlignInferJuncHash(tmpRefineSAMindex, tmpSJindex);
				bool SJvalidOrNotBool = alignInferInfoVec[tmpSJindexInAlignInferJuncHash].validSJ();
				if(!SJvalidOrNotBool)
				{
					lowSupportLowConfidenceSJ_exists_bool = true;
					break;
				}
			}
			if(lowSupportLowConfidenceSJ_exists_bool)
			{
				unique_final_sam_ofs << refineSAMvecInfo_SE->returnUnmapSAMstr(tmpRefineSAMindex) << endl;
				unique_corrected_sam_original_ofs << refineSAMvecInfo_SE->returnOriSAMstr(tmpRefineSAMindex) << endl;
				unique_corrected_sam_ofs << refineSAMvecInfo_SE->returnUnmapSAMstr(tmpRefineSAMindex) << endl; 
			}
			else
			{
				unique_final_sam_ofs << refineSAMvecInfo_SE->returnOriSAMstr(tmpRefineSAMindex) << endl;
			}
		}
	}

	void correct_output_SAM_withLowSupportLowConfidenceSJ_PE(RefineSAM_PE_Vec_Info* refineSAMvecInfo_PE,
		ofstream& unique_final_sam_ofs, 
		ofstream& unique_corrected_sam_original_ofs, ofstream& unique_corrected_sam_ofs, Index_Info* indexInfo)
	{
		for(int tmpRefineSAMindex = 0; tmpRefineSAMindex < refineSAMvecInfo_PE->returnRefineSAMvecSize(); tmpRefineSAMindex++)
		{
			if(refineSAMvecInfo_PE->point2NULLbool(tmpRefineSAMindex))
				continue;
			bool lowSupportLowConfidenceSJ_exists_bool = false;
			for(int tmpSJindex = 0; tmpSJindex < refineSAMvecInfo_PE->returnRefineSAMinfoSJvecSize(tmpRefineSAMindex); tmpSJindex++)
			{
				int tmpSJindexInAlignInferJuncHash = refineSAMvecInfo_PE->returnIndexInAlignInferJuncHash(tmpRefineSAMindex, tmpSJindex);
				bool SJvalidOrNotBool = alignInferInfoVec[tmpSJindexInAlignInferJuncHash].validSJ();
				if(!SJvalidOrNotBool)
				{
					lowSupportLowConfidenceSJ_exists_bool = true;
					break;
				}
			}
			if(lowSupportLowConfidenceSJ_exists_bool)
			{
				unique_final_sam_ofs << refineSAMvecInfo_PE->returnUnmapSAMstr(tmpRefineSAMindex) << endl;
				unique_corrected_sam_original_ofs << refineSAMvecInfo_PE->returnOriSAMstr(tmpRefineSAMindex) << endl;
				unique_corrected_sam_ofs << refineSAMvecInfo_PE->returnUnmapSAMstr(tmpRefineSAMindex) << endl; 
			}
			else
			{
				unique_final_sam_ofs << refineSAMvecInfo_PE->returnOriSAMstr(tmpRefineSAMindex) << endl;
			}
		}
	}

	void parseSAMfileVec_SE(vector<string>& inputSAMfileVec, 
		Index_Info* indexInfo, int maxReadBaseNumInPathStructure,
		string& multiSAMstr, string& unmapSAMstr, ofstream& uniqueSAM_ofs,
		string& headerSectionFileStr, RefineSAM_SE_Vec_Info* refineSAMvecInfo)
	{
		ofstream headerSection_ofs(headerSectionFileStr.c_str());
		ofstream multiSAM_ofs(multiSAMstr.c_str());
		ofstream unmapSAM_ofs(unmapSAMstr.c_str());
		for(int tmp = 0; tmp < inputSAMfileVec.size(); tmp++)
		{
			cout << "tmp SAM file: " << inputSAMfileVec[tmp] << endl;
			parseSAMfile_SE(
				inputSAMfileVec[tmp], indexInfo, maxReadBaseNumInPathStructure,
				multiSAM_ofs, unmapSAM_ofs, uniqueSAM_ofs, headerSection_ofs, refineSAMvecInfo);
		}
		multiSAM_ofs.close();
		unmapSAM_ofs.close();
		headerSection_ofs.close();
	}	

	void insertJuncFromAlignmentFileVec_storeRefineSAMinfoWithLowSupSJ_filterMultiUnmap_SE(
		vector<string>& inputSAMfileVec, Index_Info* indexInfo, int maxReadBaseNumInPathStructure,
		string& multiSAMstr, string& unmapSAMstr, ofstream& uniqueSAM_ofs,
		string& headerSectionFileStr,
		RefineSAM_SE_Vec_Info* refineSAMvecInfo)
	{
		ofstream headerSection_ofs(headerSectionFileStr.c_str());
		ofstream multiSAM_ofs(multiSAMstr.c_str());
		ofstream unmapSAM_ofs(unmapSAMstr.c_str());
		for(int tmp = 0; tmp < inputSAMfileVec.size(); tmp++)
		{
			cout << "tmp SAM file: " << inputSAMfileVec[tmp] << endl;
			insertJuncFromAlignmentFile_storeRefineSAMinfoWithLowSupSJ_filterMultiUnmap_SE(
				inputSAMfileVec[tmp], indexInfo, maxReadBaseNumInPathStructure,
				multiSAM_ofs, unmapSAM_ofs, uniqueSAM_ofs, headerSection_ofs, refineSAMvecInfo);
		}
		multiSAM_ofs.close();
		unmapSAM_ofs.close();
		headerSection_ofs.close();
	}

	void insertJuncFromAlignmentFileVec_storeRefineSAMinfoWithLowSupSJ_filterMultiUnmap_PE(
		vector<string>& inputSAMfileVec, Index_Info* indexInfo, int maxReadBaseNumInPathStructure,
		string& multiSAMstr, string& unmapSAMstr, ofstream& uniqueSAM_ofs,
		string& headerSectionFileStr,
		RefineSAM_PE_Vec_Info* refineSAMvecInfo)
	{
		ofstream headerSection_ofs(headerSectionFileStr.c_str());
		ofstream multiSAM_ofs(multiSAMstr.c_str());
		ofstream unmapSAM_ofs(unmapSAMstr.c_str());
		for(int tmp = 0; tmp < inputSAMfileVec.size(); tmp++)
		{
			cout << "tmp SAM file: " << inputSAMfileVec[tmp] << endl;
			insertJuncFromAlignmentFile_storeRefineSAMinfoWithLowSupSJ_filterMultiUnmap_PE(
				inputSAMfileVec[tmp], indexInfo, maxReadBaseNumInPathStructure,
				multiSAM_ofs, unmapSAM_ofs, uniqueSAM_ofs, headerSection_ofs, refineSAMvecInfo);
		}
		multiSAM_ofs.close();
		unmapSAM_ofs.close();
		headerSection_ofs.close();
	}	

	void parseSAMfile_SE(
		string& alignmentFile, Index_Info* indexInfo, int maxReadBaseNumInPathStructure,
		ofstream& multiSAM_ofs, ofstream& unmapSAM_ofs, ofstream& uniqueSAM_ofs, 
		ofstream& headerSection_ofs, RefineSAM_SE_Vec_Info* refineSAMvecInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(!sam_ifs.eof())
		{
			string tmpLineStr;
			getline(sam_ifs, tmpLineStr);
			if((sam_ifs.eof())||(tmpLineStr == ""))
				break;
			if(tmpLineStr.at(0) == '@')
			{
				headerSection_ofs << tmpLineStr << endl;
				continue;
			}
			RefineSAM_SE_Info* refineSamSEinfo = new RefineSAM_SE_Info();
			//int multi_or_unmap_type = UNIQUE_SAM_TYPE;
			int multi_or_unmap_type = refineSamSEinfo->checkMultiUnmap_initiate_SE(
				tmpLineStr, indexInfo, maxReadBaseNumInPathStructure);
			if(multi_or_unmap_type == MULTI_SAM_TYPE) // MULTI
			{
				multiSAM_ofs << tmpLineStr << endl;
				delete refineSamSEinfo;
			}
			else if(multi_or_unmap_type == UNMAP_SAM_TYPE) // UNMAP
			{
				unmapSAM_ofs << tmpLineStr << endl;
				delete refineSamSEinfo;
			}
			else // UNIQUE
			{
				uniqueSAM_ofs << tmpLineStr << endl;
				delete refineSamSEinfo;
			}
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_storeRefineSAMinfoWithLowSupSJ_filterMultiUnmap_SE(
		string& alignmentFile, Index_Info* indexInfo, int maxReadBaseNumInPathStructure,
		ofstream& multiSAM_ofs, ofstream& unmapSAM_ofs, ofstream& uniqueSAM_ofs, ofstream& headerSection_ofs,
		RefineSAM_SE_Vec_Info* refineSAMvecInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(!sam_ifs.eof())
		{
			string tmpLineStr;
			getline(sam_ifs, tmpLineStr);
			if((sam_ifs.eof())||(tmpLineStr == ""))
				break;
			if(tmpLineStr.at(0) == '@')
			{
				headerSection_ofs << tmpLineStr << endl;
				continue;
			}
			RefineSAM_SE_Info* refineSamSEinfo = new RefineSAM_SE_Info();
			int multi_or_unmap_type = refineSamSEinfo->checkMultiUnmap_initiate_SE(
				tmpLineStr, indexInfo, maxReadBaseNumInPathStructure);
			if(multi_or_unmap_type == MULTI_SAM_TYPE) // MULTI
			{
				multiSAM_ofs << tmpLineStr << endl;
				delete refineSamSEinfo;
			}
			else if(multi_or_unmap_type == UNMAP_SAM_TYPE) // UNMAP
			{
				unmapSAM_ofs << tmpLineStr << endl;
				delete refineSamSEinfo;
			}
			else
			{
				//cout << "unique ...." << endl;
				this->checkSJindexInAlignInferJuncHash_addNewSJ_updateOldSJ_SE(
					refineSAMvecInfo, refineSamSEinfo, indexInfo);
				//cout << "start to check withLowSupSJ_bool " << endl;
				bool withLowSupSJ_bool = this->withLowSJbool_refineSAMinfo_SE(refineSamSEinfo);
				//cout << "withLowSupSJ_bool: " << withLowSupSJ_bool << endl;
				if(withLowSupSJ_bool)
				{
					//cout << "start to do pushBackRefineSAMseInfo " << endl;
					refineSAMvecInfo->pushBackRefineSAMseInfo(refineSamSEinfo);
					//cout << "end of pushBackRefineSAMseInfo " << endl;
				}
				else
				{	
					uniqueSAM_ofs << tmpLineStr << endl;
					delete refineSamSEinfo;
				}
			}
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_storeRefineSAMinfoWithLowSupSJ_filterMultiUnmap_PE(
		string& alignmentFile, Index_Info* indexInfo, int maxReadBaseNumInPathStructure,
		ofstream& multiSAM_ofs, ofstream& unmapSAM_ofs, ofstream& uniqueSAM_ofs, ofstream& headerSection_ofs,
		RefineSAM_PE_Vec_Info* refineSAMvecInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(!sam_ifs.eof())
		{
			string tmpLineStr_1;
			getline(sam_ifs, tmpLineStr_1);
			if((sam_ifs.eof())||(tmpLineStr_1 == ""))
				break;
			if(tmpLineStr_1.at(0) == '@')
			{
				headerSection_ofs << tmpLineStr_1 << endl;
				continue;
			}
			string tmpLineStr_2;
			getline(sam_ifs, tmpLineStr_2);
			RefineSAM_PE_Info* refineSamPEinfo = new RefineSAM_PE_Info();
			int multi_or_unmap_type = refineSamPEinfo->checkMultiUnmap_initiate_PE(
				tmpLineStr_1, tmpLineStr_2, indexInfo, maxReadBaseNumInPathStructure);
			if(multi_or_unmap_type == MULTI_SAM_TYPE) // MULTI
			{
				multiSAM_ofs << tmpLineStr_1 << endl;
				multiSAM_ofs << tmpLineStr_2 << endl;
				delete refineSamPEinfo;
			}
			else if(multi_or_unmap_type == UNMAP_SAM_TYPE) // UNMAP
			{
				unmapSAM_ofs << tmpLineStr_1 << endl;
				unmapSAM_ofs << tmpLineStr_2 << endl;
				delete refineSamPEinfo;
			}
			else
			{
				//cout << "unique ...." << endl;
				this->checkSJindexInAlignInferJuncHash_addNewSJ_updateOldSJ_PE(
					refineSAMvecInfo, refineSamPEinfo, indexInfo);
				//cout << "start to check withLowSupSJ_bool " << endl;
				bool withLowSupSJ_bool = this->withLowSJbool_refineSAMinfo_PE(refineSamPEinfo);
				//cout << "withLowSupSJ_bool: " << withLowSupSJ_bool << endl;
				if(withLowSupSJ_bool)
				{
					//cout << "start to do pushBackRefineSAMseInfo " << endl;
					refineSAMvecInfo->pushBackRefineSAMseInfo(refineSamPEinfo);
					//cout << "end of pushBackRefineSAMseInfo " << endl;
				}
				else
				{	
					uniqueSAM_ofs << tmpLineStr_1 << endl;
					uniqueSAM_ofs << tmpLineStr_2 << endl;
					delete refineSamPEinfo;
				}
			}
		}
		sam_ifs.close();
	}

	bool withLowSJbool_refineSAMinfo_SE(RefineSAM_SE_Info* tmpRefineSAMinfo)
	{
		//cout << "start to do withLowSJbool_refineSAMinfo..." << endl;
		int tmpSJvecSizeInRefineSAMinfo = tmpRefineSAMinfo->returnSJvecSize();
		//cout << "tmpSJvecSizeInRefineSAMinfo: " << tmpSJvecSizeInRefineSAMinfo << endl;
		for(int tmp = 0; tmp < tmpSJvecSizeInRefineSAMinfo; tmp++)
		{
			//cout << "tmpSJindex: " << tmp << endl;
			int tmpSJindexInAlignInferJuncHash_inRefineSAMinfo
				= tmpRefineSAMinfo->returnSJindexInAlignInferJuncHash(tmp);
			//cout << "tmpSJindexInAlignInferJuncHash_inRefineSAMinfo: " << tmpSJindexInAlignInferJuncHash_inRefineSAMinfo << endl;
			//cout << "tmpAlignInferJuncVecSize: " << alignInferInfoVec.size() << endl;
			bool tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool // SJ support num no higher than STORED_SUPPORT_READ_NUM_MAX 
				= alignInferInfoVec[tmpSJindexInAlignInferJuncHash_inRefineSAMinfo].supportNumNoHigherThanStoredReadNumMaxBool();
			//cout << "tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool: " << tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool << endl;
			if(tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool)
				return true;
		}
		return false;
	}

	bool withLowSJbool_refineSAMinfo_PE(RefineSAM_PE_Info* tmpRefineSAMinfo)
	{
		//cout << "start to do withLowSJbool_refineSAMinfo..." << endl;
		int tmpSJvecSizeInRefineSAMinfo = tmpRefineSAMinfo->returnSJvecSize();
		//cout << "tmpSJvecSizeInRefineSAMinfo: " << tmpSJvecSizeInRefineSAMinfo << endl;
		for(int tmp = 0; tmp < tmpSJvecSizeInRefineSAMinfo; tmp++)
		{
			//cout << "tmpSJindex: " << tmp << endl;
			int tmpSJindexInAlignInferJuncHash_inRefineSAMinfo
				= tmpRefineSAMinfo->returnSJindexInAlignInferJuncHash(tmp);
			//cout << "tmpSJindexInAlignInferJuncHash_inRefineSAMinfo: " << tmpSJindexInAlignInferJuncHash_inRefineSAMinfo << endl;
			//cout << "tmpAlignInferJuncVecSize: " << alignInferInfoVec.size() << endl;
			bool tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool // SJ support num no higher than STORED_SUPPORT_READ_NUM_MAX 
				= alignInferInfoVec[tmpSJindexInAlignInferJuncHash_inRefineSAMinfo].supportNumNoHigherThanStoredReadNumMaxBool();
			//cout << "tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool: " << tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool << endl;
			if(tmpSJ_supportNumNoHigherThanMAXrefineSAM2store_bool)
				return true;
		}
		return false;
	}

	void checkSJindexInAlignInferJuncHash_addNewSJ_updateOldSJ_SE(
		RefineSAM_SE_Vec_Info* refineSamSEvecInfo,
		RefineSAM_SE_Info* refineSamSEinfo, Index_Info* indexInfo)
	{
		//cout << "checkSJindexInAlignInferJuncHash_addNewSJ_updateOldSJ_SE starts ..." << endl;
		int tmpSJ_chrNameInt = refineSamSEinfo->returnChrNameInt();
		//cout << "tmpSJ_chrNameInt: " << tmpSJ_chrNameInt << endl;
		//cout << "tmpSJ num: " << refineSamSEinfo->returnSJvecSize() << endl;
		for(int tmp = 0; tmp < refineSamSEinfo->returnSJvecSize(); tmp++)
		{
			//cout << "tmpSJindex: " << tmp << endl;
			int tmpSJ_donerEndPos = refineSamSEinfo->returnSJ_donerEndPos(tmp);
			int tmpSJ_acceptorStartPos = refineSamSEinfo->returnSJ_acceptorStartPos(tmp);
			//cout << "tmpSJ_donerEndPos: " << tmpSJ_donerEndPos << endl;
			//cout << "tmpSJ_acceptorStartPos: " << tmpSJ_acceptorStartPos << endl;
			int tmpSJindexInAlignInferJuncHash 
				= this->searchAndReturnAlignInferInfoVecIndex(tmpSJ_chrNameInt,
					tmpSJ_donerEndPos, tmpSJ_acceptorStartPos);
			// if((tmpSJ_donerEndPos == 10045)&&(tmpSJ_acceptorStartPos == 10808))
			// {
			// 	cout << endl << "SJ: 10045 ~ 10808 detected" << endl;
			// 	cout << "tmpSJ_chrNameInt: " << tmpSJ_chrNameInt << endl;
			// 	cout << "tmpSJindexInAlignInferJuncHash: " << tmpSJindexInAlignInferJuncHash << endl;
			// 	cout << "alignmentSAMstr: " << endl << refineSamSEinfo->returnOriSAMstr() << endl;
 		// 	}

			//cout << "tmpSJindexInAlignInferJuncHash: " << tmpSJindexInAlignInferJuncHash << endl;
			if(tmpSJindexInAlignInferJuncHash < 0) // new SJ
			{
				// add new SJ 2 alignInferJuncHash
				//cout << "tmpAlignInferJuncVecSize: " << alignInferInfoVec.size() << endl;
				int newSJindexInAlignInferJuncHash = currentAlignInferInfoVecIndex;
				//cout << "newSJindexInAlignInferJuncHash: " << newSJindexInAlignInferJuncHash << endl;
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_withRefineSAMinfo_SE(
					refineSamSEinfo, tmp, indexInfo);
				tmpAlignInferInfo.pushBackNewSupportRefineSAMindex_SE(refineSamSEvecInfo);
				//cout << "end of initiate " << endl;

				AlignInferJunctionMap::iterator alignInferJuncMapIter 
					= alignInferJunctionMapVec[tmpSJ_chrNameInt].find(tmpSJ_donerEndPos);
				if(alignInferJuncMapIter == (alignInferJunctionMapVec[tmpSJ_chrNameInt]).end())
				{
					AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
					tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
						tmpSJ_acceptorStartPos, newSJindexInAlignInferJuncHash));
					alignInferJunctionMapVec[tmpSJ_chrNameInt].insert(
						pair< int, AcceptorStartPos2AlignInferInfoMap> (
							tmpSJ_donerEndPos, tmpAcceptorStartPos2AlignInferInfoMap));
				}
				else
				{
					(alignInferJuncMapIter->second).insert(pair<int, int>(
						tmpSJ_acceptorStartPos, newSJindexInAlignInferJuncHash));					
				}
				//cout << "end of insertAlignInferInfo 2 map" << endl;					
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				//cout << "end of pushBackRefineSAMseInfo" << endl;
				// pushBack alignInferJuncHashIndex
				refineSamSEinfo->pushBackAlignInferJuncIndex(newSJindexInAlignInferJuncHash);
				//cout << "end of pushBackAlignInferJuncIndex " << endl;
				//cout << "tmpAlignInferJuncVecSize: " << alignInferInfoVec.size() << endl;
				currentAlignInferInfoVecIndex ++;
			}
			else
			{
				//cout << "start to update old SJs" << endl;
				int tmpSJsupportNum_old = alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnSupportNum();
				//cout << "tmpSJsupportNum_old: " << tmpSJsupportNum_old << endl;
				//cout << "start to do updateOldAlignInferJunc_withRefineSAMinfo_SE" << endl;
				alignInferInfoVec[tmpSJindexInAlignInferJuncHash].updateOldAlignInferJunc_withRefineSAMinfo_SE(
						refineSamSEinfo, tmp, indexInfo);
				//cout << "end of updateOldAlignInferJunc_withRefineSAMinfo_SE" << endl;
				if(tmpSJsupportNum_old < STORED_SUPPORT_READ_NUM_MAX) // very low support
				{
					alignInferInfoVec[tmpSJindexInAlignInferJuncHash].pushBackNewSupportRefineSAMindex_SE(refineSamSEvecInfo);					
				}
				else if(tmpSJsupportNum_old == STORED_SUPPORT_READ_NUM_MAX) // full support reads to store, release 
				{
					this->releaseStoredSupportRefineSAMinfoIfpossible_SE(
						tmpSJindexInAlignInferJuncHash, refineSamSEvecInfo);
				}
				else // High support num
				{}
				// pushBack alignInferJuncHashIndex
				refineSamSEinfo->pushBackAlignInferJuncIndex(tmpSJindexInAlignInferJuncHash);
			}
		}
	}

	void checkSJindexInAlignInferJuncHash_addNewSJ_updateOldSJ_PE(
		RefineSAM_PE_Vec_Info* refineSamPEvecInfo,
		RefineSAM_PE_Info* refineSamPEinfo, Index_Info* indexInfo)
	{
		//cout << "checkSJindexInAlignInferJuncHash_addNewSJ_updateOldSJ_SE starts ..." << endl;
		int tmpSJ_chrNameInt = refineSamPEinfo->returnChrNameInt();
		//cout << "tmpSJ_chrNameInt: " << tmpSJ_chrNameInt << endl;
		//cout << "tmpSJ num: " << refineSamPEinfo->returnSJvecSize() << endl;
		for(int tmp = 0; tmp < refineSamPEinfo->returnSJvecSize(); tmp++)
		{
			//cout << "tmpSJindex: " << tmp << endl;
			int tmpSJ_donerEndPos = refineSamPEinfo->returnSJ_donerEndPos(tmp);
			int tmpSJ_acceptorStartPos = refineSamPEinfo->returnSJ_acceptorStartPos(tmp);
			//cout << "tmpSJ_donerEndPos: " << tmpSJ_donerEndPos << endl;
			//cout << "tmpSJ_acceptorStartPos: " << tmpSJ_acceptorStartPos << endl;
			int tmpSJindexInAlignInferJuncHash 
				= this->searchAndReturnAlignInferInfoVecIndex(tmpSJ_chrNameInt,
					tmpSJ_donerEndPos, tmpSJ_acceptorStartPos);
			//cout << "tmpSJindexInAlignInferJuncHash: " << tmpSJindexInAlignInferJuncHash << endl;
			if(tmpSJindexInAlignInferJuncHash < 0) // new SJ
			{
				// add new SJ 2 alignInferJuncHash
				//cout << "tmpAlignInferJuncVecSize: " << alignInferInfoVec.size() << endl;
				int newSJindexInAlignInferJuncHash = currentAlignInferInfoVecIndex;
				//cout << "newSJindexInAlignInferJuncHash: " << newSJindexInAlignInferJuncHash << endl;
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_withRefineSAMinfo_PE(
					refineSamPEinfo, tmp, indexInfo);
				tmpAlignInferInfo.pushBackNewSupportRefineSAMindex_PE(refineSamPEvecInfo);
				//cout << "end of initiate " << endl;

				AlignInferJunctionMap::iterator alignInferJuncMapIter 
					= alignInferJunctionMapVec[tmpSJ_chrNameInt].find(tmpSJ_donerEndPos);
				if(alignInferJuncMapIter == (alignInferJunctionMapVec[tmpSJ_chrNameInt]).end())
				{					
					AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
					tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
						tmpSJ_acceptorStartPos, newSJindexInAlignInferJuncHash));
					alignInferJunctionMapVec[tmpSJ_chrNameInt].insert(
						pair< int, AcceptorStartPos2AlignInferInfoMap> (
							tmpSJ_donerEndPos, tmpAcceptorStartPos2AlignInferInfoMap));
				}
				else
				{
					(alignInferJuncMapIter->second).insert(pair<int, int>(
						tmpSJ_acceptorStartPos, newSJindexInAlignInferJuncHash));					
				}
				//cout << "end of insertAlignInferInfo 2 map" << endl;	
				
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				//cout << "end of pushBackRefineSAMseInfo" << endl;
				// pushBack alignInferJuncHashIndex
				refineSamPEinfo->pushBackAlignInferJuncIndex(newSJindexInAlignInferJuncHash);
				//cout << "end of pushBackAlignInferJuncIndex " << endl;
				//cout << "tmpAlignInferJuncVecSize: " << alignInferInfoVec.size() << endl;
				currentAlignInferInfoVecIndex ++;
			}
			else
			{
				//cout << "start to update old SJs" << endl;
				int tmpSJsupportNum_old = alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnSupportNum();
				//cout << "tmpSJsupportNum_old: " << tmpSJsupportNum_old << endl;
				//cout << "start to do updateOldAlignInferJunc_withRefineSAMinfo_SE" << endl;
				alignInferInfoVec[tmpSJindexInAlignInferJuncHash].updateOldAlignInferJunc_withRefineSAMinfo_PE(
						refineSamPEinfo, tmp, indexInfo);
				//cout << "end of updateOldAlignInferJunc_withRefineSAMinfo_SE" << endl;
				if(tmpSJsupportNum_old < STORED_SUPPORT_READ_NUM_MAX) // very low support
				{
					alignInferInfoVec[tmpSJindexInAlignInferJuncHash].pushBackNewSupportRefineSAMindex_PE(refineSamPEvecInfo);					
				}
				else if(tmpSJsupportNum_old == STORED_SUPPORT_READ_NUM_MAX) // full support reads to store, release 
				{
					this->releaseStoredSupportRefineSAMinfoIfpossible_PE(
						tmpSJindexInAlignInferJuncHash, refineSamPEvecInfo);
				}
				else // High support num
				{}
				// pushBack alignInferJuncHashIndex
				refineSamPEinfo->pushBackAlignInferJuncIndex(tmpSJindexInAlignInferJuncHash);
			}
		}
	}

	void releaseStoredSupportRefineSAMinfoIfpossible_SE(int tmpSJindexInAlignInferJuncHash,
		RefineSAM_SE_Vec_Info* refineSamSEvecInfo)
	{
		// cout << "tmpSJindexInAlignInferJuncHash: " << tmpSJindexInAlignInferJuncHash << endl;
		// cout << "start to do releaseStoredSupportRefineSAMinfoIfpossible_SE" << endl;
		// cout << "tmpSupportNum: " << alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnSupportNum() << endl;
		// cout << "tmpStoredRefineSAMinfoVecSize: " << alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnStoredRefineSAMinfoVecSize() << endl;
		for(int tmpSupportReadIndex = 0; 
			tmpSupportReadIndex < alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnStoredRefineSAMinfoVecSize(); 
			tmpSupportReadIndex++)
		{
			//cout << "tmpSupportReadIndex: " << tmpSupportReadIndex << endl;
			int tmpIndexInRefineSAMinfoVec 
				= alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnRefineSAMinfoIndex(tmpSupportReadIndex);
			//cout << "tmpIndexInRefineSAMinfoVec: " << tmpIndexInRefineSAMinfoVec << endl;
			int tmpRefineSAMinfoSJvecSize 
				= refineSamSEvecInfo->returnRefineSAMinfoSJvecSize(tmpIndexInRefineSAMinfoVec);
			//cout << "tmpRefineSAMinfoSJvecSize: " << tmpRefineSAMinfoSJvecSize << endl;
			bool allSJwithHighSupportNumBool = true;
			for(int tmp = 0; tmp < tmpRefineSAMinfoSJvecSize; tmp++)
			{
				//cout << "tmp: " << tmp << endl;
				int tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash
					= refineSamSEvecInfo->returnIndexInAlignInferJuncHash(
						tmpIndexInRefineSAMinfoVec, tmp);
				//cout << "tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash: " << tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash << endl; 
				if(tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash == tmpSJindexInAlignInferJuncHash)
					continue;
				int tmpRefineSAMinfoSJvec_tmpSJsupportNum 
					= alignInferInfoVec[tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash].returnSupportNum();
				//cout << "tmpRefineSAMinfoSJvec_tmpSJsupportNum: " << tmpRefineSAMinfoSJvec_tmpSJsupportNum  << endl;
				if(tmpRefineSAMinfoSJvec_tmpSJsupportNum <= STORED_SUPPORT_READ_NUM_MAX)
				{
					allSJwithHighSupportNumBool = false;
					break;
				}
			}
			// cout << "end of checking tmpRefineSAMinfo SJs" << endl;
			// if(tmpIndexInRefineSAMinfoVec == 225)
			// {
			// 	cout << "alignInferJuncHash:" << 73 << " " << alignInferInfoVec[73].returnSupportNum() << endl;
			// 	cout << "alignInferJuncHash:" << 171 << " " << alignInferInfoVec[171].returnSupportNum() << endl;
			// 	cout << "alignInferJuncHash:" << 172 << " " << alignInferInfoVec[172].returnSupportNum() << endl;
			// }
			// all SJ with high support numbers
			if(allSJwithHighSupportNumBool)
				refineSamSEvecInfo->deleteRefineSAMinfo(tmpIndexInRefineSAMinfoVec);
			else
			{}
		}
	}

	void releaseStoredSupportRefineSAMinfoIfpossible_PE(int tmpSJindexInAlignInferJuncHash,
		RefineSAM_PE_Vec_Info* refineSamPEvecInfo)
	{
		// cout << "tmpSJindexInAlignInferJuncHash: " << tmpSJindexInAlignInferJuncHash << endl;
		// cout << "start to do releaseStoredSupportRefineSAMinfoIfpossible_SE" << endl;
		// cout << "tmpSupportNum: " << alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnSupportNum() << endl;
		// cout << "tmpStoredRefineSAMinfoVecSize: " << alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnStoredRefineSAMinfoVecSize() << endl;
		for(int tmpSupportReadIndex = 0; 
			tmpSupportReadIndex < alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnStoredRefineSAMinfoVecSize(); 
			tmpSupportReadIndex++)
		{
			//cout << "tmpSupportReadIndex: " << tmpSupportReadIndex << endl;
			int tmpIndexInRefineSAMinfoVec 
				= alignInferInfoVec[tmpSJindexInAlignInferJuncHash].returnRefineSAMinfoIndex(tmpSupportReadIndex);
			//cout << "tmpIndexInRefineSAMinfoVec: " << tmpIndexInRefineSAMinfoVec << endl;
			int tmpRefineSAMinfoSJvecSize 
				= refineSamPEvecInfo->returnRefineSAMinfoSJvecSize(tmpIndexInRefineSAMinfoVec);
			//cout << "tmpRefineSAMinfoSJvecSize: " << tmpRefineSAMinfoSJvecSize << endl;
			bool allSJwithHighSupportNumBool = true;
			for(int tmp = 0; tmp < tmpRefineSAMinfoSJvecSize; tmp++)
			{
				//cout << "tmp: " << tmp << endl;
				int tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash
					= refineSamPEvecInfo->returnIndexInAlignInferJuncHash(
						tmpIndexInRefineSAMinfoVec, tmp);
				//cout << "tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash: " << tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash << endl; 
				if(tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash == tmpSJindexInAlignInferJuncHash)
					continue;
				int tmpRefineSAMinfoSJvec_tmpSJsupportNum 
					= alignInferInfoVec[tmpRefineSAMinfoSJvec_tmpSJindexInAlignInferJuncHash].returnSupportNum();
				//cout << "tmpRefineSAMinfoSJvec_tmpSJsupportNum: " << tmpRefineSAMinfoSJvec_tmpSJsupportNum  << endl;
				if(tmpRefineSAMinfoSJvec_tmpSJsupportNum <= STORED_SUPPORT_READ_NUM_MAX)
				{
					allSJwithHighSupportNumBool = false;
					break;
				}
			}
			// cout << "end of checking tmpRefineSAMinfo SJs" << endl;
			// if(tmpIndexInRefineSAMinfoVec == 225)
			// {
			// 	cout << "alignInferJuncHash:" << 73 << " " << alignInferInfoVec[73].returnSupportNum() << endl;
			// 	cout << "alignInferJuncHash:" << 171 << " " << alignInferInfoVec[171].returnSupportNum() << endl;
			// 	cout << "alignInferJuncHash:" << 172 << " " << alignInferInfoVec[172].returnSupportNum() << endl;
			// }
			// all SJ with high support numbers
			if(allSJwithHighSupportNumBool)
				refineSamPEvecInfo->deleteRefineSAMinfo(tmpIndexInRefineSAMinfoVec);
			else
			{}
		}
	}

	bool searchAndReturnSupNumInAlignInferJuncHash(int& supNumInThisJuncHash, 
		int tmpChrNameInt, int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		int foundIndex = this->searchAndReturnAlignInferInfoVecIndex(
			tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
		if(foundIndex < 0)
			return false;
		else
		{
			supNumInThisJuncHash = this->returnAlignInferInfo_supportNum(foundIndex);
			return true;
		}
	}

	bool searchAndReturnSupNumInAlignInferJuncHash(int& supNumInThisJuncHash, 
		string& tmpChrName, int tmpDonerEndPos, int tmpAcceptorStartPos, 
		Index_Info* indexInfo)
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
			return false;
		bool foundOrNotBool = this->searchAndReturnSupNumInAlignInferJuncHash(
			supNumInThisJuncHash, tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
		return foundOrNotBool;
	}

	bool foundInAlignInferJunctionHash(int tmpChrInt,
		int tmpDonerEndPos, int tmpAcceptorStartPos, int offset)
	{
		//bool found_bool = false;
		for(int tmp = 0-offset; tmp <= offset; tmp++)
		{
			for(int tmp2 = 0-offset; tmp2 <= offset; tmp2++)
			{
				int tmpNewDonerEndPos = tmpDonerEndPos + tmp;
				int tmpNewAcceptorStartPos = tmpAcceptorStartPos + tmp2;
				int foundIndex = this->searchAndReturnAlignInferInfoVecIndex(
					tmpChrInt, tmpNewDonerEndPos, tmpNewAcceptorStartPos);
				if(foundIndex != -1)
					return true;
			}
		}
		return false;
	}

	bool exactTheSame2someAlterSpliceSite_doner(int index)
	{
		return alignInferInfoVec[index].exactTheSame2someAlterSpliceSite_doner();
	}

	bool exactTheSame2someAlterSpliceSite_acceptor(int index)
	{
		return alignInferInfoVec[index].exactTheSame2someAlterSpliceSite_acceptor();
	}

	bool SJcanBeExtendedAtAcceptorSpliceSite(int index)
	{
		return alignInferInfoVec[index].canBeExtended_doner_bool();
	}
	
	bool SJcanBeExtendedAtDonerSpliceSite(int index)
	{
		return alignInferInfoVec[index].canBeExtended_acceptor_bool();
	}

	bool SJfoundAndInvalid(int tmpChrInt, int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		int SJindex = this->searchAndReturnAlignInferInfoVecIndex(
			tmpChrInt, tmpDonerEndPos, tmpAcceptorStartPos);
		if(SJindex < 0)
			return false;
		else
		{
			return (this->SJinvalid(SJindex));
		}
	}

	bool SJinvalid(int index)
	{
		bool tmpSJvalid = alignInferInfoVec[index].validSJ();
		return (!tmpSJvalid);
	}

	bool SJexistInAlignInferJuncHash(int tmpChrInt, int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		int tmpIndexInAlignerInferInfoVec = this->searchAndReturnAlignInferInfoVecIndex(
			tmpChrInt, tmpDonerEndPos, tmpAcceptorStartPos);
		if(tmpIndexInAlignerInferInfoVec < 0)
			return false;
		else
			return true;
	}

	bool someJuncInVecExistInAlignInferJuncHash(int tmpChrInt,
		vector< pair<int, int> >& tmpSJposPairVec)
	{
		for(int tmp = 0; tmp < tmpSJposPairVec.size(); tmp++)
		{
			//int tmpChrInt = tmpChrIntVec[tmp];
			int tmpDonerEndPos = tmpSJposPairVec[tmp].first;
			int tmpAcceptorStartPos = tmpSJposPairVec[tmp].second;
			bool tmpSJexistInAlignInferJuncHashBool = this->SJexistInAlignInferJuncHash(
				tmpChrInt, tmpDonerEndPos, tmpAcceptorStartPos);
			if(tmpSJexistInAlignInferJuncHashBool)
				return true;
		}
		return false;
	}

	int searchAndReturnAlignInferInfoVecIndex(int tmpChrInt,
		int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		//cout << "searchAndReturnAlignInferInfoVecIndex starts ......" << endl;
		AlignInferJunctionMap::iterator alignInferJunctionMapIter
			= alignInferJunctionMapVec[tmpChrInt].find(tmpDonerEndPos);
		if(alignInferJunctionMapIter != alignInferJunctionMapVec[tmpChrInt].end())
		{
			//cout << "tmpDonerEndPos found ......" << endl;
			AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPos2AlignInferInfoMapIter
				= (alignInferJunctionMapIter->second).find(tmpAcceptorStartPos);
			if(acceptorStartPos2AlignInferInfoMapIter != (alignInferJunctionMapIter->second).end())
			{
				//cout << "tmpAcceptorStartPos found ....." << endl;
				return (acceptorStartPos2AlignInferInfoMapIter->second);
			}
			else
			{
				//cout << "tmpAcceptorStartPos not found ....." << endl;
				//cout << "not found in searchAndReturnAlignInferInfoVecIndex 2nd level" << endl;
				return -1;
				//exit(1);					
			}
		}
		else
		{
			//cout << "tmpDonerEndPos not found ..." << endl;
			//cout << "not found in searchAndReturnAlignInferInfoVecIndex 1st level" << endl;
			return -1;
			//exit(1);			
		}
	}

	bool generateInferedUnfixedHeadPath_alignInferHash(
		vector<Jump_Code>& tmpInferedPathJumpCodeVec, 
		vector<int>& tmpInferedPathMismatchPosVec, 
		vector<char>& tmpInferedPathMismatchCharVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedHeadLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr)
	{
		//cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\ngenerateInferedUnfixedHeadPath_alignInferHash starts ..." << endl;
		int foundIndex = this->searchAndReturnAlignInferInfoVecIndex(tmpChromNameInt,
			tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);
		if(foundIndex < 0)
			return false;
		//cout << "foundIndex: " << foundIndex << endl;
		bool tmpInferedPathGeneratedBool
			= alignInferInfoVec[foundIndex].generateInferedUnfixedHeadPath_alignInfer(
				tmpInferedPathJumpCodeVec, 
				tmpInferedPathMismatchPosVec, 
				tmpInferedPathMismatchCharVec,
				indexInfo, tmpChromNameInt, 
				unfixedHeadLength, readSeqWithDirection,
				tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
				tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);
		//cout << "tmpInferedPathJumpCodeVecSize: " << tmpInferedPathJumpCodeVec.size() << endl
		//	<< "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl ; 
		return tmpInferedPathGeneratedBool;
	}

	
	bool generateInferedUnfixedTailPath_alignInferHash(
		vector<Jump_Code>& tmpInferedPathJumpCodeVec, 
		vector<int>& tmpInferedPathMismatchPosVec, 
		vector<char>& tmpInferedPathMismatchCharVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedTailLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr)
	{
		int foundIndex = this->searchAndReturnAlignInferInfoVecIndex(tmpChromNameInt,
			tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);
		//cout << "foundIndex: " << foundIndex << endl;
		if(foundIndex < 0)
			return false;
		bool tmpInferedPathGeneratedBool
			= alignInferInfoVec[foundIndex].generateInferedUnfixedTailPath_alignInfer(
				tmpInferedPathJumpCodeVec, 
				tmpInferedPathMismatchPosVec, 
				tmpInferedPathMismatchCharVec,
				indexInfo, tmpChromNameInt, 
				unfixedTailLength, readSeqWithDirection,
				tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
				tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);
		//cout << "tmpInferedPathGeneratedBool: " << tmpInferedPathGeneratedBool << endl;
		return tmpInferedPathGeneratedBool;
	}

	/*
	bool extendHeadWithAlignInferInfo_backward(int tmpChrInt,
		int tmpDonerEndPos, int tmpAcceptorStartPos, const string& readSeq_inProcess,
		int unfixedHeadSeqLen, Index_Info* indexInfo, vector<Jump_Code>& extendedHeadJumpCodeVec,
		vector<int>& extendedHeadMismatchPosVecInRead, vector<char>& extendedHeadMismatchCharVec) // FIX ME: 04 / 15 / 2015
	{	
		cout << "extendHeadWithAlignInferInfo_backward( .. starts ...." << endl;
		cout << "try to search and return AlignInferInfoVecIndex ..." << endl;
		int foundIndex = this->searchAndReturnAlignInferInfoVecIndex(tmpChrInt, tmpDonerEndPos, tmpAcceptorStartPos);
		cout << "foundIndex: " << foundIndex << endl;
		bool extend_bool = alignInferInfoVec[foundIndex].extendHeadWithAlignInferInfo_backward(
			readSeq_inProcess, unfixedHeadSeqLen, indexInfo, extendedHeadJumpCodeVec,
			extendedHeadMismatchPosVecInRead, extendedHeadMismatchCharVec);
		return extend_bool;
	}*/

	void initiateAlignInferJunctionHashInfo(int chrNum)
	{
		//cout << " alignInferJunctionMapVec.size(): " << alignInferJunctionMapVec.size() << endl;
		for(int tmp = 0; tmp < chrNum; tmp++)
		{
			AlignInferJunctionMap newAlignJunctionMap;
			alignInferJunctionMapVec.push_back(newAlignJunctionMap);
		}
		currentAlignInferInfoVecIndex = 0;
		//cout << " alignInferJunctionMapVec.size(): " << alignInferJunctionMapVec.size() << endl;	
	}

	void initiateAlignInferJunctionInfo(int chrNum)
	{
		//cout << " alignInferJunctionMapVec.size(): " << alignInferJunctionMapVec.size() << endl;
		for(int tmp = 0; tmp < chrNum; tmp++)
		{
			AlignInferJunctionMap newAlignJunctionMap;
			alignInferJunctionMapVec.push_back(newAlignJunctionMap);
		}
		currentAlignInferInfoVecIndex = 0;
		//cout << " alignInferJunctionMapVec.size(): " << alignInferJunctionMapVec.size() << endl;	
	}

	void insertJuncFromAlignmentFileVec(vector<string>& alignmentFileVec, Index_Info* indexInfo,
		int maxReadBaseNumInPathStructure)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_maxReadBaseNumInPathStructure(
				alignmentFileVec[tmp], indexInfo, maxReadBaseNumInPathStructure);
		}
	}

	void insertJuncFromAlignmentFileVec_chrNamePosOnly(
		vector<string>& alignmentFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePosOnly(
				alignmentFileVec[tmp], indexInfo);
		}
	}	

	void insertJuncFromAlignmentFileVec_chrNamePos_supportNum_anchorSize_XM(
		vector<string>& alignmentFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePos_suportNum_anchorSize_XM(
				alignmentFileVec[tmp], indexInfo);
		}
	}

	void insertJuncFromAlignmentFileVec_chrNamePos_supportNum_anchorSize(
		vector<string>& alignmentFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePos_suportNum_anchorSize(
				alignmentFileVec[tmp], indexInfo);
		}
	}

	void insertJuncFromAlignmentFileVec_chrNamePos_supportNum(
		vector<string>& alignmentFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePos_suportNum(
				alignmentFileVec[tmp], indexInfo);
		}
	}		

	void insertJuncFromAlignmentFile_chrNamePosOnly(
		string& alignmentFile, Index_Info* indexInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			//if(sam_ifs.eof())
			//	break;
			//cout << "tmpAlignStr: " << tmpAlignStr << endl;
			if(tmpAlignStr == "")
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			this->getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePosOnly(
				tmpAlignStr, indexInfo);				
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_chrNamePos_suportNum_anchorSize_XM(
		string& alignmentFile, Index_Info* indexInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			//if(sam_ifs.eof())
			//	break;
			//cout << "tmpAlignStr: " << tmpAlignStr << endl;
			if(tmpAlignStr == "")
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			this->getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize_XM(
				tmpAlignStr, indexInfo);				
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_chrNamePos_suportNum_anchorSize(
		string& alignmentFile, Index_Info* indexInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			//if(sam_ifs.eof())
			//	break;
			//cout << "tmpAlignStr: " << tmpAlignStr << endl;
			if(tmpAlignStr == "")
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			this->getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize(
				tmpAlignStr, indexInfo);				
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_chrNamePos_suportNum(
		string& alignmentFile, Index_Info* indexInfo)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			//if(sam_ifs.eof())
			//	break;
			//cout << "tmpAlignStr: " << tmpAlignStr << endl;
			if(tmpAlignStr == "")
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			this->getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum(
				tmpAlignStr, indexInfo);				
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_maxReadBaseNumInPathStructure(
		string& alignmentFile, Index_Info* indexInfo, int maxReadBaseNumInPathStructure)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			//if(sam_ifs.eof())
			//	break;
			//cout << "tmpAlignStr: " << tmpAlignStr << endl;
			if(tmpAlignStr == "")
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			this->getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_maxReadBaseNumInPathStruct(
				tmpAlignStr, indexInfo, maxReadBaseNumInPathStructure);				
		}
		sam_ifs.close();
	}

	// void insertJuncFromAlignmentFileVec_storePeReadAlignInfoWithLowSupSJ(vector<string>& alignmentFileVec, Index_Info* indexInfo,
	// 	int maxReadBaseNumInPathStructure)
	// {
	// 	vector<PeReadAlign2Refine_Info*> peReadAlign2RefineInfoVec;
	// 	int peReadAlign2RefineInfoVec_currentIndex = 0;
	// 	for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
	// 	{
	// 		cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
	// 		insertJuncFromAlignmentFile_storePeReadAlignInfoWithLowSupSJ(
	// 			alignmentFileVec[tmp], indexInfo, maxReadBaseNumInPathStructure,
	// 			peReadAlign2RefineInfoVec, peReadAlign2RefineInfoVec_currentIndex);
	// 	}
	// }	

	// void insertJuncFromAlignmentFile_storePeReadAlignInfoWithLowSupSJ(
	// 	string& alignmentFile, Index_Info* indexInfo, int maxReadBaseNumInPathStructure,
	// 	vector<PeReadAlign2Refine_Info*>& peReadAlign2RefineInfoVec, 
	// 	int& peReadAlign2RefineInfoVec_currentIndex)
	// {
	// 	ifstream sam_ifs(alignmentFile.c_str());
	// 	string firstNonHeaderLineStr;
	// 	while(1)
	// 	{
	// 		string tmpHeaderLineStr;
	// 		getline(sam_ifs, tmpHeaderLineStr);
	// 		if(tmpHeaderLineStr.at(0) == '@')
	// 			continue;
	// 		else
	// 		{	
	// 			firstNonHeaderLineStr = tmpHeaderLineStr;
	// 			break;
	// 		}
	// 	}
	// 	string firstReadName = this->getReadNameFromAlignStr(firstNonHeaderLineStr);
	// 	vector<string> peReadAlignmentStrVec;
	// 	string readName_old;
	// 	while(1)
	// 	{
	// 		if(sam_ifs.eof())
	// 		{
	// 			PeReadAlign2Refine_Info* tmpPeReadAlign2RefineInfo
	// 				= new tmpPeReadAlign2RefineInfo();
	// 			tmpPeReadAlign2RefineInfo->initiateWithAlignStrVec(peReadAlignmentStrVec);
	// 			tmpPeReadAlign2RefineInfo->checkSJconfidence();
	// 			bool readCrossingLowConfiSJ_insert2alignInferJuncHash_bool
	// 				= this->insertPeReadAlign2RefineInfoIfCrossingLowConfidenceSJ(tmpPeReadAlign2RefineInfo);
	// 			if(!readCrossingLowConfiSJ_insert2alignInferJuncHash_bool)
	// 				delete tmpPeReadAlign2RefineInfo;
	// 			break;
	// 		}
	// 		string tmpAlignStr;
	// 		getline(sam_ifs, tmpAlignStr);
	// 		if(tmpAlignStr == "")
	// 		{
	// 			PeReadAlign2Refine_Info* tmpPeReadAlign2RefineInfo
	// 				= new tmpPeReadAlign2RefineInfo();
	// 			tmpPeReadAlign2RefineInfo->initiateWithAlignStrVec(peReadAlignmentStrVec);
	// 			tmpPeReadAlign2RefineInfo->checkSJconfidence();
	// 			bool readCrossingLowConfiSJ_insert2alignInferJuncHash_bool
	// 				= this->insertPeReadAlign2RefineInfoIfCrossingLowConfidenceSJ(tmpPeReadAlign2RefineInfo);
	// 			if(!readCrossingLowConfiSJ_insert2alignInferJuncHash_bool)
	// 				delete tmpPeReadAlign2RefineInfo;
	// 			break;
	// 		}

	// 		string tmpReadName = this->getReadNameFromAlignStr(tmpAlignStr);
	// 		if(tmpReadName == readName_old)
	// 		{
	// 			peReadAlignmentStrVec.push_back(tmpAlignStr);
	// 		}
	// 		else
	// 		{
	// 			PeReadAlign2Refine_Info* tmpPeReadAlign2RefineInfo
	// 				= new tmpPeReadAlign2RefineInfo();
	// 			tmpPeReadAlign2RefineInfo->initiateWithAlignStrVec(peReadAlignmentStrVec);
	// 			tmpPeReadAlign2RefineInfo->checkSJconfidence();
	// 			bool readCrossingLowConfiSJ_insert2alignInferJuncHash_bool
	// 				= this->insertPeReadAlign2RefineInfoIfCrossingLowConfidenceSJ(tmpPeReadAlign2RefineInfo);
	// 			if(!readCrossingLowConfiSJ_insert2alignInferJuncHash_bool)
	// 				delete tmpPeReadAlign2RefineInfo;
	// 			peReadAlignmentStrVec.clear();
	// 			peReadAlignmentStrVec.push_back(tmpAlignStr);
	// 			readName_old = tmpReadName;			
	// 		}
	// 	}
	// 	sam_ifs.close();
	// }

	void getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePosOnly(
		const string& samStr, Index_Info* indexInfo)
	{
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];

		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);

		vector< pair<int,int> > tmpSJposPairVec;
		vector< int > tmpSJindexVec_cigarStringJumpCodeVec;
		this->generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, 
			tmpSJposPairVec, tmpSJindexVec_cigarStringJumpCodeVec);
		this->insertSJvec_chrNamePosOnly(mapChrNameInt, tmpSJposPairVec, indexInfo);
	}

	void getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum(
		const string& samStr, Index_Info* indexInfo)
	{
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];

		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);

		vector< pair<int,int> > tmpSJposPairVec;
		vector< int > tmpSJindexVec_cigarStringJumpCodeVec;
		this->generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, 
			tmpSJposPairVec, tmpSJindexVec_cigarStringJumpCodeVec);
		this->insertSJvec_chrNamePos_supportNum(mapChrNameInt, tmpSJposPairVec, indexInfo);
	}

	void getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize(
		const string& samStr, Index_Info* indexInfo)
	{
		// get SJposVec from SAM string
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];
		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);		
		vector< pair<int,int> > tmpSJposPairVec;
		vector< int > tmpSJindexVec_cigarStringJumpCodeVec;
		this->generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, 
			tmpSJposPairVec, tmpSJindexVec_cigarStringJumpCodeVec);
		vector< pair<int,int> > tmpSJanchorSizePairVec;
		this->generateSJanchorSizeVec(tmpSJindexVec_cigarStringJumpCodeVec,
			tmpSJanchorSizePairVec, cigarStringJumpCodeVec);
		this->insertSJvec_chrNamePos_supportNum_anchorSize(mapChrNameInt, tmpSJposPairVec,
			tmpSJanchorSizePairVec, indexInfo);
	}

	void getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize_XM(
		const string& samStr, Index_Info* indexInfo)
	{
		// get SJposVec from SAM string
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];
		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);		
		vector< pair<int,int> > tmpSJposPairVec;
		vector< int > tmpSJindexVec_cigarStringJumpCodeVec;
		this->generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, 
			tmpSJposPairVec, tmpSJindexVec_cigarStringJumpCodeVec);
		vector< pair<int,int> > tmpSJanchorSizePairVec;
		this->generateSJanchorSizeVec(tmpSJindexVec_cigarStringJumpCodeVec,
			tmpSJanchorSizePairVec, cigarStringJumpCodeVec);
		
		int tmpXMint;
		if(tmpSJanchorSizePairVec.size() > 0)
		{
			int XM_loc = samStr.find("XM:i:",0);
			if(XM_loc == string::npos)
				tmpXMint = 0;
			else
			{	
				int XM_nextTab_loc = samStr.find("\t", XM_loc+1);
				string XMfield = samStr.substr(XM_loc, XM_nextTab_loc - 1 - XM_loc + 1);
				string XMintStr = XMfield.substr(5);
				tmpXMint = atoi(XMintStr.c_str());
			}
		}		
		this->insertSJvec_chrNamePos_supportNum_anchorSize_XM(mapChrNameInt, tmpSJposPairVec,
			tmpSJanchorSizePairVec, tmpXMint, indexInfo);
	}

	void getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_maxReadBaseNumInPathStruct(
		const string& samStr, Index_Info* indexInfo, int maxReadBaseNum)
	{
		// get SJposVec from SAM string
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];

		// cout << endl <<  "mapChrNameStr: " << mapChrNameStr << endl;
		// cout << "mapChrPos: " << mapChrPos << endl;
		// cout << "cigarString: " << cigarString << endl;

		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		// if(cigarStringJumpCodeVec[0].type == "I")
		// 	cout << "Alignment: " << endl << samStr << endl;
		// for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		// {
		// 	cout << "jumpCode: " << cigarStringJumpCodeVec[tmp].toString() << endl;
		// }

		vector< pair<int,int> > tmpSJposPairVec;
		vector< int > tmpSJindexVec_cigarStringJumpCodeVec;
		this->generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, 
			tmpSJposPairVec, tmpSJindexVec_cigarStringJumpCodeVec);
		vector< pair<int,int> > tmpSJanchorSizePairVec;
		this->generateSJanchorSizeVec(tmpSJindexVec_cigarStringJumpCodeVec,
			tmpSJanchorSizePairVec, cigarStringJumpCodeVec);
		// for(int tmp = 0; tmp < tmpSJposPairVec.size(); tmp++)
		// {
		// 	cout << "SJ: " << tmpSJposPairVec[tmp].first << " ~ " << tmpSJposPairVec[tmp].second << endl;
		// } 

		vector< vector<Jump_Code> > tmpSJjumpCodeVecVec_backward;
		vector< vector<Jump_Code> > tmpSJjumpCodeVecVec_forward;
		this->generateSJjumpCodeVec_maxReadBaseNumInPathStruct(
			cigarStringJumpCodeVec, tmpSJindexVec_cigarStringJumpCodeVec,
			tmpSJjumpCodeVecVec_backward, tmpSJjumpCodeVecVec_forward,
			maxReadBaseNum);

		this->insertSJvecWithJumpCodeVec_maxReadBaseNumInPathStruct_maxAnchorSize(mapChrNameInt, tmpSJposPairVec,
			tmpSJjumpCodeVecVec_backward, tmpSJjumpCodeVecVec_forward, indexInfo, maxReadBaseNum, tmpSJanchorSizePairVec);
	}

	void generateSJanchorSizeVec(vector<int>& SJindexVec_cigarStringJumpCodeVec,
		vector< pair<int,int> >& SJanchorSizePairVec, vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int SJindexVecSize = SJindexVec_cigarStringJumpCodeVec.size();
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		if(SJindexVecSize == 0)
		{}
		else if(SJindexVecSize == 1)
		{
			int tmpDonerAnchorFirstJumpCodeIndex = 0;
			int tmpDonerAnchorLastJumpCodeIndex = SJindexVec_cigarStringJumpCodeVec[0]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex = SJindexVec_cigarStringJumpCodeVec[0]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex = cigarStringJumpCodeVecSize - 1;
			int tmpDonerAnchorSize = this->getEndLocInReadOfSpecificJumpCode(
				cigarStringJumpCodeVec, tmpDonerAnchorLastJumpCodeIndex);
			int tmpAcceptorAnchorSize
				= this->getEndLocInReadOfSpecificJumpCode(
					cigarStringJumpCodeVec, tmpAcceptorAnchorLastJumpCodeIndex)
					- this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpAcceptorAnchorFirstJumpCodeIndex-1);
			// cout << "tmpDonerAnchorSize: " << tmpDonerAnchorSize << endl;
			// cout << "tmpAcceptorAnchorSize: " << tmpAcceptorAnchorSize << endl;
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize, tmpAcceptorAnchorSize));
		}
		else if(SJindexVecSize == 2)
		{
			int tmpDonerAnchorFirstJumpCodeIndex_1stSJ = 0;
			int tmpDonerAnchorLastJumpCodeIndex_1stSJ = SJindexVec_cigarStringJumpCodeVec[0]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ = SJindexVec_cigarStringJumpCodeVec[0]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_1stSJ = SJindexVec_cigarStringJumpCodeVec[1]-1;
			int tmpDonerAnchorSize_1stSJ = this->getEndLocInReadOfSpecificJumpCode(
				cigarStringJumpCodeVec, tmpDonerAnchorLastJumpCodeIndex_1stSJ);
			int tmpAcceptorAnchorSize_1stSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					cigarStringJumpCodeVec, tmpAcceptorAnchorLastJumpCodeIndex_1stSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ-1);
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_1stSJ, tmpAcceptorAnchorSize_1stSJ));
			int tmpAcceptorAnchorFirstJumpCodeIndex_2ndSJ = SJindexVec_cigarStringJumpCodeVec[1]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_2ndSJ = cigarStringJumpCodeVecSize - 1;
			int tmpDonerAnchorSize_2ndSJ = tmpAcceptorAnchorSize_1stSJ;
			int tmpAcceptorAnchorSize_2ndSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					cigarStringJumpCodeVec, tmpAcceptorAnchorLastJumpCodeIndex_2ndSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpAcceptorAnchorFirstJumpCodeIndex_2ndSJ-1);
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_2ndSJ, tmpAcceptorAnchorSize_2ndSJ));
		}
		else if(SJindexVecSize > 2)
		{
			// 1st SJ
			int tmpDonerAnchorFirstJumpCodeIndex_1stSJ = 0;
			int tmpDonerAnchorLastJumpCodeIndex_1stSJ = SJindexVec_cigarStringJumpCodeVec[0]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ = SJindexVec_cigarStringJumpCodeVec[0]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_1stSJ = SJindexVec_cigarStringJumpCodeVec[1]-1;
			int tmpDonerAnchorSize_1stSJ = this->getEndLocInReadOfSpecificJumpCode(
				cigarStringJumpCodeVec, tmpDonerAnchorLastJumpCodeIndex_1stSJ);
			int tmpAcceptorAnchorSize_1stSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					cigarStringJumpCodeVec, tmpAcceptorAnchorLastJumpCodeIndex_1stSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ-1);
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_1stSJ, tmpAcceptorAnchorSize_1stSJ));
			// other SJs
			for(int tmpSJ = 1; tmpSJ <= SJindexVecSize-2; tmpSJ++)
			{
				int tmpDonerAnchorFirstJumpCodeIndex_tmpSJ = SJindexVec_cigarStringJumpCodeVec[tmpSJ-1]+1;
				int tmpDonerAnchorLastJumpCodeIndex_tmpSJ = SJindexVec_cigarStringJumpCodeVec[tmpSJ]-1;
				int tmpAcceptorAnchorFirstJumpCodeIndex_tmpSJ = SJindexVec_cigarStringJumpCodeVec[tmpSJ]+1;
				int tmpAcceptorAnchorLastJumpCodeIndex_tmpSJ = SJindexVec_cigarStringJumpCodeVec[tmpSJ+1]-1;
				int tmpDonerAnchorSize_tmpSJ 
					= this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpDonerAnchorLastJumpCodeIndex_tmpSJ);
						- this->getEndLocInReadOfSpecificJumpCode(
							cigarStringJumpCodeVec, tmpDonerAnchorFirstJumpCodeIndex_tmpSJ-1);
				int tmpAcceptorAnchorSize_tmpSJ 
					= this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpAcceptorAnchorLastJumpCodeIndex_tmpSJ)
						- this->getEndLocInReadOfSpecificJumpCode(
							cigarStringJumpCodeVec, tmpAcceptorAnchorFirstJumpCodeIndex_tmpSJ-1);
				SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_tmpSJ, tmpAcceptorAnchorSize_tmpSJ));							
			}
			// last SJ
			int tmpDonerAnchorFirstJumpCodeIndex_lastSJ = SJindexVec_cigarStringJumpCodeVec[SJindexVecSize-2]+1;
			int tmpDonerAnchorLastJumpCodeIndex_lastSJ = SJindexVec_cigarStringJumpCodeVec[SJindexVecSize-1]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex_lastSJ = SJindexVec_cigarStringJumpCodeVec[SJindexVecSize-1]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_lastSJ = cigarStringJumpCodeVecSize - 1;
			int tmpDonerAnchorSize_lastSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					cigarStringJumpCodeVec, tmpDonerAnchorLastJumpCodeIndex_lastSJ);
					- this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpDonerAnchorFirstJumpCodeIndex_lastSJ-1);
			int tmpAcceptorAnchorSize_lastSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					cigarStringJumpCodeVec, tmpAcceptorAnchorLastJumpCodeIndex_lastSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						cigarStringJumpCodeVec, tmpAcceptorAnchorFirstJumpCodeIndex_lastSJ-1); 					
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_lastSJ, tmpAcceptorAnchorSize_lastSJ));
		}
		else
		{}
	}

	void insertSJvec_chrNamePosOnly(int tmpSJchrNameInt,
		vector< pair<int,int> >& tmpSJposPairVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < tmpSJposPairVec.size(); tmp++)
		{
			this->insertSJ_chrNamePosOnly(
				tmpSJchrNameInt, tmpSJposPairVec[tmp].first,
				tmpSJposPairVec[tmp].second, indexInfo);
		}		
	}

	void insertSJvec_chrNamePos_supportNum(int tmpSJchrNameInt,
		vector< pair<int,int> >& tmpSJposPairVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < tmpSJposPairVec.size(); tmp++)
		{
			this->insertSJ_chrNamePos_supportNum(
				tmpSJchrNameInt, tmpSJposPairVec[tmp].first,
				tmpSJposPairVec[tmp].second, indexInfo);
		}		
	}

	void insertSJvec_chrNamePos_supportNum_anchorSize(int tmpSJchrNameInt,
		vector< pair<int,int> >& tmpSJposPairVec, 
		vector< pair<int,int> >& tmpSJanchorSizePairVec,
		Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < tmpSJposPairVec.size(); tmp++)
		{
			this->insertSJ_chrNamePos_supportNum_anchorSize(
				tmpSJchrNameInt, 
				tmpSJposPairVec[tmp].first,
				tmpSJposPairVec[tmp].second, 
				tmpSJanchorSizePairVec[tmp].first,
				tmpSJanchorSizePairVec[tmp].second,
				indexInfo);
		}		
	}

	void insertSJvec_chrNamePos_supportNum_anchorSize_XM(
		int tmpSJchrNameInt,
		vector< pair<int,int> >& tmpSJposPairVec, 
		vector< pair<int,int> >& tmpSJanchorSizePairVec,
		int tmpXMint, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < tmpSJposPairVec.size(); tmp++)
		{
			this->insertSJ_chrNamePos_supportNum_anchorSize_XM(
				tmpSJchrNameInt, 
				tmpSJposPairVec[tmp].first,
				tmpSJposPairVec[tmp].second, 
				tmpSJanchorSizePairVec[tmp].first,
				tmpSJanchorSizePairVec[tmp].second,
				tmpXMint, indexInfo);
		}		
	}	

	void insertSJvecWithJumpCodeVec_maxReadBaseNumInPathStruct_maxAnchorSize(int mapChrNameInt, 
			vector< pair<int,int> >& tmpSJposPairVec,
			vector< vector<Jump_Code> >& tmpSJjumpCodeVecVec_backward, 
			vector< vector<Jump_Code> >& tmpSJjumpCodeVecVec_forward, 
			Index_Info* indexInfo, int maxReadBaseNum, vector< pair<int,int> >& tmpSJanchorSizePairVec)
	{
		for(int tmp = 0; tmp < tmpSJposPairVec.size(); tmp++)
		{
			//cout << "tmpSj: " << tmp+ 1 << " start: " << tmpSJposPairVec[tmp].first << " end: " << tmpSJposPairVec[tmp].second << endl;
			this->insertSJwithJumpCodeVec_maxReadBaseNumInPathStruct_maxAnchorSize(
				mapChrNameInt, tmpSJposPairVec[tmp].first,
				tmpSJposPairVec[tmp].second, tmpSJjumpCodeVecVec_backward[tmp], 
				tmpSJjumpCodeVecVec_forward[tmp], indexInfo, maxReadBaseNum,
				tmpSJanchorSizePairVec[tmp].first, tmpSJanchorSizePairVec[tmp].second);
		}
	}

	void insertSJ_chrNamePosOnly(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart,
		Index_Info* indexInfo)
	{
		AlignInferJunctionMap::iterator alignInferJuncMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
		{
			//cout << "not found" << endl;
			AlignInfer_Info tmpAlignInferInfo;
			tmpAlignInferInfo.initiateAlignInferInfo_chrNamePosOnly(
				mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, indexInfo);
			//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
			alignInferInfoVec.push_back(tmpAlignInferInfo);
			currentAlignInferInfoVecIndex ++;
		}
		else // old donerEndPos found, search for acceptorStartPos
		{
			AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
				= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
			//cout << "found" << endl;
			if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
			{
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePosOnly(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, indexInfo);				
				(alignInferJuncMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else
			{	
				// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
				// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
				// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
				// //alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
			}
		}
	}

	void insertSJ_chrNamePos_supportNum(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart,
		Index_Info* indexInfo)
	{
		AlignInferJunctionMap::iterator alignInferJuncMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
		{
			//cout << "not found" << endl;
			AlignInfer_Info tmpAlignInferInfo;
			tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum(
				mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, indexInfo);
			//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
			alignInferInfoVec.push_back(tmpAlignInferInfo);
			currentAlignInferInfoVecIndex ++;
		}
		else // old donerEndPos found, search for acceptorStartPos
		{
			AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
				= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
			//cout << "found" << endl;
			if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
			{
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, indexInfo);				
				(alignInferJuncMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else
			{	
				// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
				// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
				// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
				alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
			}
		}
	}

	void insertSJ_chrNamePos_strand(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart,
		string& tmpStrand, Index_Info* indexInfo)
	{
		AlignInferJunctionMap::iterator alignInferJuncMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
		{
			//cout << "not found" << endl;
			AlignInfer_Info tmpAlignInferInfo;
			tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum(
				mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, tmpStrand, indexInfo);
			//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
			alignInferInfoVec.push_back(tmpAlignInferInfo);
			currentAlignInferInfoVecIndex ++;
		}
		else // old donerEndPos found, search for acceptorStartPos
		{
			AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
				= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
			//cout << "found" << endl;
			if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
			{
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, tmpStrand, indexInfo);				
				(alignInferJuncMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else
			{	
				// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
				// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
				// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
				alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
			}
		}
	}

	void insertSJ_chrNamePos_supportNum_anchorSize(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart,
		int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart,
		Index_Info* indexInfo)
	{
		AlignInferJunctionMap::iterator alignInferJuncMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
		{
			//cout << "not found" << endl;
			AlignInfer_Info tmpAlignInferInfo;
			tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize(
				mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
				tmpSJanchorSizeDonerEnd, tmpSJanchorSizeAcceptorStart,
				indexInfo);
			//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
			alignInferInfoVec.push_back(tmpAlignInferInfo);
			currentAlignInferInfoVecIndex ++;
		}
		else // old donerEndPos found, search for acceptorStartPos
		{
			AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
				= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
			//cout << "found" << endl;
			if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
			{
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
					tmpSJanchorSizeDonerEnd, tmpSJanchorSizeAcceptorStart, indexInfo);				
				(alignInferJuncMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else
			{	
				// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
				// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
				// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
				alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
				alignInferInfoVec[acceptorStartPosMapIter->second].updateAnchorSize(
					tmpSJanchorSizeDonerEnd, tmpSJanchorSizeAcceptorStart);
			}
		}
	}

	void insertSJ_chrNamePos_supportNum_anchorSize_XM(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart,
		int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart,
		int tmpXMint, Index_Info* indexInfo)
	{
		AlignInferJunctionMap::iterator alignInferJuncMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
		{
			//cout << "not found" << endl;
			AlignInfer_Info tmpAlignInferInfo;
			tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize_XM(
				mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
				tmpSJanchorSizeDonerEnd, tmpSJanchorSizeAcceptorStart,
				tmpXMint, indexInfo);
			//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
			alignInferInfoVec.push_back(tmpAlignInferInfo);
			currentAlignInferInfoVecIndex ++;
		}
		else // old donerEndPos found, search for acceptorStartPos
		{
			AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
				= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
			//cout << "found" << endl;
			if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
			{
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize_XM(
					mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
					tmpSJanchorSizeDonerEnd, tmpSJanchorSizeAcceptorStart, 
					tmpXMint, indexInfo);				
				(alignInferJuncMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else
			{	
				// alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
				// 	tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
				// 	tmpDonerAnchorSize, tmpAcceptorAnchorSize);
				alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
				alignInferInfoVec[acceptorStartPosMapIter->second].updateAnchorSize(
					tmpSJanchorSizeDonerEnd, tmpSJanchorSizeAcceptorStart);
				alignInferInfoVec[acceptorStartPosMapIter->second].updateXM(tmpXMint);
			}
		}
	}

	void insertSJwithJumpCodeVec_maxReadBaseNumInPathStruct_maxAnchorSize(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		vector<Jump_Code>& tmpSJjumpCodeVec_backward, 
		vector<Jump_Code>& tmpSJjumpCodeVec_forward, 
		Index_Info* indexInfo, int maxReadBaseNum, 
		int tmpDonerAnchorSize, int tmpAcceptorAnchorSize)
	{
		AlignInferJunctionMap::iterator alignInferJuncMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if( alignInferJuncMapIter == (alignInferJunctionMapVec[mapChrNameInt]).end() ) // new donerEndPos to insert
		{
			//cout << "not found" << endl;
			AlignInfer_Info tmpAlignInferInfo;
			tmpAlignInferInfo.initiateAlignInferInfo_maxAnchorSize(mapChrNameInt, tmpSJposDonerEnd, 
				tmpSJposAcceptorStart, tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward,
				tmpDonerAnchorSize, tmpAcceptorAnchorSize, indexInfo);
			//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
			alignInferInfoVec.push_back(tmpAlignInferInfo);
			currentAlignInferInfoVecIndex ++;
		}
		else // old donerEndPos found, search for acceptorStartPos
		{
			AcceptorStartPos2AlignInferInfoMap::iterator acceptorStartPosMapIter 
				= (alignInferJuncMapIter->second).find(tmpSJposAcceptorStart);
			//cout << "found" << endl;
			if(acceptorStartPosMapIter == (alignInferJuncMapIter->second).end()) // new acceptor pos found, insert it
			{
				AlignInfer_Info tmpAlignInferInfo;
				tmpAlignInferInfo.initiateAlignInferInfo_maxAnchorSize(mapChrNameInt, tmpSJposDonerEnd, 
					tmpSJposAcceptorStart, tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward,
					tmpDonerAnchorSize, tmpAcceptorAnchorSize, indexInfo);				
				(alignInferJuncMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
				alignInferInfoVec.push_back(tmpAlignInferInfo);
				currentAlignInferInfoVecIndex ++;
			}
			else
			{	
				alignInferInfoVec[acceptorStartPosMapIter->second].updateOldAlignInfo_maxAnchorSize(
					tmpSJjumpCodeVec_backward, tmpSJjumpCodeVec_forward, maxReadBaseNum,
					tmpDonerAnchorSize, tmpAcceptorAnchorSize);
				//alignInferInfoVec[acceptorStartPosMapIter->second].updateSupportNum();
			}
		}
	}

	void generateSJjumpCodeVec_maxReadBaseNumInPathStruct(
			vector<Jump_Code>& cigarStringJumpCodeVec,
			vector< int >& tmpSJindexVec_cigarStringJumpCodeVec,
			vector< vector<Jump_Code> >& SJjumpCodeVecVec_backward,
			vector< vector<Jump_Code> >& SJjumpCodeVecVec_forward, int maxReadBaseNum)
	{
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		//cout << "cigarStringJumpCodeVecSize: " << cigarStringJumpCodeVecSize << endl;
		for(int tmp = 0; tmp < tmpSJindexVec_cigarStringJumpCodeVec.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int tmpSJindex_cigarStringJumpCodeVec = tmpSJindexVec_cigarStringJumpCodeVec[tmp];
			//cout << "tmpSJindex_cigarStringJumpCodeVec " << tmpSJindex_cigarStringJumpCodeVec << endl;
			vector<Jump_Code> tmpSJjumpCodeVec_backward;
			vector<Jump_Code> tmpSJjumpCodeVec_forward;

			this->generateSJjumpCodeVec_backward_maxReadBaseNumInPathStruct(
				cigarStringJumpCodeVec, tmpSJindex_cigarStringJumpCodeVec,
				tmpSJjumpCodeVec_backward, maxReadBaseNum);
			this->generateSJjumpCodeVec_forward_maxReadBaseNumInPathStruct(
				cigarStringJumpCodeVec, tmpSJindex_cigarStringJumpCodeVec,
				tmpSJjumpCodeVec_forward, maxReadBaseNum);

			SJjumpCodeVecVec_backward.push_back(tmpSJjumpCodeVec_backward);
			SJjumpCodeVecVec_forward.push_back(tmpSJjumpCodeVec_forward);
		}
	}	

	void generateSJjumpCodeVec_backward_maxReadBaseNumInPathStruct(
		vector<Jump_Code>& cigarStringJumpCodeVec, 
		int tmpSJindex_cigarStringJumpCodeVec,
		vector<Jump_Code>& tmpSJjumpCodeVec_backward, 
		int maxReadBaseNum)
	{
		int tmpReadBaseNum = 0;
		for(int tmpJumpCodeIndex = tmpSJindex_cigarStringJumpCodeVec-1; tmpJumpCodeIndex >= 0; tmpJumpCodeIndex--)
		{
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpJumpCodeIndex].len;
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpJumpCodeIndex].type;
			if(tmpJumpCodeType == "S")
				return;
			if((tmpJumpCodeType == "I")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "S"))
				tmpReadBaseNum = tmpReadBaseNum + tmpJumpCodeLength;
			if(tmpReadBaseNum < maxReadBaseNum)
			{
				tmpSJjumpCodeVec_backward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
			}
			else if(tmpReadBaseNum == maxReadBaseNum)
			{
				tmpSJjumpCodeVec_backward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
				return;
			}
			else
			{
				int overJumpCodeLength = tmpReadBaseNum - maxReadBaseNum;
				int lastJumpCodeLength = tmpJumpCodeLength - overJumpCodeLength;
				Jump_Code lastJumpCode(lastJumpCodeLength, tmpJumpCodeType);
				tmpSJjumpCodeVec_backward.push_back(lastJumpCode);
				return;
			}
		}
	}

	void generateSJjumpCodeVec_forward_maxReadBaseNumInPathStruct(
		vector<Jump_Code>& cigarStringJumpCodeVec, 
		int tmpSJindex_cigarStringJumpCodeVec,
		vector<Jump_Code>& tmpSJjumpCodeVec_forward,
		int maxReadBaseNum)
	{
		int tmpReadBaseNum = 0;
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		for(int tmpJumpCodeIndex = tmpSJindex_cigarStringJumpCodeVec+1; 
			tmpJumpCodeIndex < cigarStringJumpCodeVecSize; tmpJumpCodeIndex++)
		{
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpJumpCodeIndex].len;
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpJumpCodeIndex].type;
			if(tmpJumpCodeType == "S")
				return;
			if((tmpJumpCodeType == "I")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "S"))
				tmpReadBaseNum = tmpReadBaseNum + tmpJumpCodeLength;
			if(tmpReadBaseNum < maxReadBaseNum)
			{
				tmpSJjumpCodeVec_forward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
			}
			else if(tmpReadBaseNum == maxReadBaseNum)
			{
				tmpSJjumpCodeVec_forward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
				return;
			}
			else
			{
				int overJumpCodeLength = tmpReadBaseNum - maxReadBaseNum;
				int lastJumpCodeLength = tmpJumpCodeLength - overJumpCodeLength;
				Jump_Code lastJumpCode(lastJumpCodeLength, tmpJumpCodeType);
				tmpSJjumpCodeVec_forward.push_back(lastJumpCode);
				return;
			}
		}
	}

	/*
	void insertJuncFromJuncFileVec(vector<string>& juncFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < juncFileVec.size(); tmp++)
		{
			//cout << "tmpAlignmentFile: " << juncFileVec[tmp] << endl;
			insertJuncFromJuncFile(juncFileVec[tmp], indexInfo);
		}		
	}*/

	void insertJuncFromJuncFile(string& juncFile, Index_Info* indexInfo)
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			if(tmpJuncStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash(tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}	

	void insertJuncFromJuncFileVec_chrNamePosOnly(vector<string>& juncFileVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < juncFileVec.size(); tmp++)
		{
			string juncFile = juncFileVec[tmp];	
			ifstream junc_ifs(juncFile.c_str());
			while(1)
			{
				if(junc_ifs.eof())
					break;
				string tmpJuncStr;
				getline(junc_ifs, tmpJuncStr);
				if(junc_ifs.eof())
					break;
				if(tmpJuncStr.substr(0,3) != "chr")
					continue;
				//cout << "start to input SAM string: " << tmpJuncStr << endl;
				this->getAlignInferInfoFromJuncFile_InsertIntoAlignInferJunctionHash_chrNamePosOnly(
					tmpJuncStr, indexInfo);				
			}
			junc_ifs.close();
		}		
	}

	void insertJuncFromJuncFile_chrNamePosOnly(string& juncFile, Index_Info* indexInfo)
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			if(tmpJuncStr.substr(0,3) != "chr")
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJuncFile_InsertIntoAlignInferJunctionHash_chrNamePosOnly(
				tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();		
	}

	void insertJuncFromJuncFile_onlyChrNamePos(string& juncFile, Index_Info* indexInfo)
	// there are only chrName and pos (and junc name, supportNum in some cases) in the file, 
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			if(tmpJuncStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_onlyChrNamePos(tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}

	void insertJuncFromJuncFile_chrNamePos_supportNum_backSpliceOnly(string& juncFile, Index_Info* indexInfo)
	// there are only chrName and pos (and junc name, supportNum in some cases) in the file, 
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			if(tmpJuncStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_backSpliceOnly(
				tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}		

	void insertJuncFromJuncFile_chrNamePos_supportNum_includeBackSplice(string& juncFile, Index_Info* indexInfo)
	// there are only chrName and pos (and junc name, supportNum in some cases) in the file, 
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			if(tmpJuncStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_includeBackSplice(
				tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}		

	void insertJuncFromJuncFile_chrNamePos_supportNum(string& juncFile, Index_Info* indexInfo)
	// there are only chrName and pos (and junc name, supportNum in some cases) in the file, 
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			if(tmpJuncStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum(tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}		

	void insertJuncFromJuncFile_chrNamePos_supportNum_flankStringChange(string& juncFile, Index_Info* indexInfo)
	// there are only chrName and pos (and junc name, supportNum in some cases) in the file, 
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			if(tmpJuncStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_flankStringChange(tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}

	void insertJuncFromJuncFileVec_chrNamePos_supportNum(vector<string>& juncFileVec, Index_Info* indexInfo)
	// there are only chrName and pos (and junc name, supportNum in some cases) in the file, 
	{
		for(int tmp = 0; tmp < juncFileVec.size(); tmp++)
		{	
			string juncFile = juncFileVec[tmp];
			ifstream junc_ifs(juncFile.c_str());
			while(1)
			{
				if(junc_ifs.eof())
					break;
				string tmpJuncStr;
				getline(junc_ifs, tmpJuncStr);
				if(junc_ifs.eof())
					break;
				if(tmpJuncStr.at(0) == '@')
					continue;
				//cout << "start to input SAM string: " << tmpJuncStr << endl;
				this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum(tmpJuncStr, indexInfo);				
			}
			junc_ifs.close();
		}
	}		

	void insertJuncFromJuncFile_chrNamePos_supportNum_anchorSize(string& juncFile, Index_Info* indexInfo)
	// there are only chrName and pos (and junc name, supportNum, anchorSize in some cases) in the file, 
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof())
				break;
			//if(tmpJuncStr.at(0) == '@')
			//	continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize(tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}		

	void insertJuncFromJuncFile_withAlterSpliceSiteAnchorSimilarity(string& juncFile, Index_Info* indexInfo)
	{
		ifstream junc_ifs(juncFile.c_str());
		while(1)
		{
			if(junc_ifs.eof())
				break;
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if(junc_ifs.eof()||(tmpJuncStr == ""))
				break;
			if(tmpJuncStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpJuncStr << endl;
			this->getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_withAlterSpliceSiteAnchorSimilarity(
				tmpJuncStr, indexInfo);				
		}
		junc_ifs.close();
	}	

	void getAlignInferInfoFromJuncFile_InsertIntoAlignInferJunctionHash_chrNamePosOnly(
		const string& juncStr, Index_Info* indexInfo)
	{
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				cout << "error in getAlignInferInfoFromJuncFile_InsertIntoAlignInferJunctionHash_chrNamePosOnly" << endl;
				exit(1);
			}
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		int nextTabLoc = juncStr.find("\t", startLoc);
		if(nextTabLoc == string::npos)
			juncFieldVec.push_back(juncStr.substr(startLoc));
		else
			juncFieldVec.push_back(juncStr.substr(startLoc, nextTabLoc-startLoc));

		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		this->insertSJ_chrNamePosOnly(juncChrNameInt, 
			juncDonerEndPos, juncAcceptorStartPos, indexInfo);
	}


	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash(
		const string& juncStr, Index_Info* indexInfo)
	{
		//cout << "getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash starts ..." << endl;
		//cout << "juncStr: " << juncStr << endl;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 16; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				cout << "error in getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash" << endl;
				exit(1);
			}
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));

		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		string juncSupportNumStr = juncFieldVec[3];
		int juncSupportNum = atoi(juncSupportNumStr.c_str());

		string flankString = juncFieldVec[4];
		string flankStringCaseStr = juncFieldVec[5];
		int flankStringCase = atoi(flankStringCaseStr.c_str());

		string juncNameStr = juncFieldVec[6];

		string juncDonerAnchorSizeMaxStr = juncFieldVec[7];
		string juncAcceptorAnchorSizeMaxStr = juncFieldVec[8];
		int juncDonerAnchorSizeMax = atoi(juncDonerAnchorSizeMaxStr.c_str());
		int juncAcceptorAnchorSizeMax = atoi(juncAcceptorAnchorSizeMaxStr.c_str());

		string juncDonerJumpCodeVecStr = juncFieldVec[9];
		string juncAcceptorJumpCodeVecStr = juncFieldVec[10];
		string donerPathVecStr = juncFieldVec[11];
		string acceptorPathVecStr = juncFieldVec[12]; 
		int donerPathVecNum = atoi(donerPathVecStr.c_str());
		int acceptorPathVecNum = atoi(acceptorPathVecStr.c_str());
		string donerPathSupportNumVecStr = juncFieldVec[13];
		string acceptorPathSupportNumVecStr = juncFieldVec[14]; 
		string donerPathMaxAnchorSizeVecStr = juncFieldVec[15];
		//cout << "donerPathMaxAnchorSizeVecStr: " << donerPathMaxAnchorSizeVecStr << endl;
		string acceptorPathMaxAnchorSizeVecStr = juncStr.substr(startLoc);//juncFieldVec[14];
		//cout << "acceptorPathMaxAnchorSizeVecStr: " << acceptorPathMaxAnchorSizeVecStr << endl;

		vector< vector<Jump_Code> > donerPathVec;
		vector< vector<Jump_Code> > acceptorPathVec;
		startLoc = 0;
		for(int tmp = 0; tmp < donerPathVecNum; tmp++)
		{	
			int commaLoc = juncDonerJumpCodeVecStr.find(",", startLoc);
			string tmpJumpCodeVecField = juncDonerJumpCodeVecStr.substr(startLoc, commaLoc-startLoc);
			vector<Jump_Code> tmpCigarStringJumpCodeVec_backward;
			this->cigarString2jumpCodeVec(tmpJumpCodeVecField, tmpCigarStringJumpCodeVec_backward);
			donerPathVec.push_back(tmpCigarStringJumpCodeVec_backward);
			startLoc = commaLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < acceptorPathVecNum; tmp++)
		{	
			int commaLoc = juncAcceptorJumpCodeVecStr.find(",", startLoc);
			string tmpJumpCodeVecField = juncAcceptorJumpCodeVecStr.substr(startLoc, commaLoc-startLoc);
			vector<Jump_Code> tmpCigarStringJumpCodeVec_forward;
			this->cigarString2jumpCodeVec(tmpJumpCodeVecField, tmpCigarStringJumpCodeVec_forward);
			acceptorPathVec.push_back(tmpCigarStringJumpCodeVec_forward);
			startLoc = commaLoc + 1;
		}

		vector<int> pathSupportNumVec_doner;
		vector<int> pathSupportNumVec_acceptor;		
		startLoc = 0;
		for(int tmp = 0; tmp < donerPathVecNum; tmp++)
		{	
			int commaLoc = donerPathSupportNumVecStr.find(",", startLoc);
			string tmpPathSupportNumStr = donerPathSupportNumVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathSupportNumInt = atoi(tmpPathSupportNumStr.c_str());
			pathSupportNumVec_doner.push_back(tmpPathSupportNumInt);
			startLoc = commaLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < acceptorPathVecNum; tmp++)
		{	
			int commaLoc = acceptorPathSupportNumVecStr.find(",", startLoc);
			string tmpPathSupportNumStr = acceptorPathSupportNumVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathSupportNumInt = atoi(tmpPathSupportNumStr.c_str());
			pathSupportNumVec_acceptor.push_back(tmpPathSupportNumInt);
			startLoc = commaLoc + 1;
		}

		vector<int> donerPathMaxAnchorSizeVec;
		vector<int> acceptorPathMaxAnchorSizeVec;
		startLoc = 0;
		for(int tmp = 0; tmp < donerPathVecNum; tmp++)
		{	
			int commaLoc = donerPathMaxAnchorSizeVecStr.find(",", startLoc);
			string tmpPathMaxAnchorSizeStr = donerPathMaxAnchorSizeVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathMaxAnchorSizeInt = atoi(tmpPathMaxAnchorSizeStr.c_str());
			donerPathMaxAnchorSizeVec.push_back(tmpPathMaxAnchorSizeInt);
			startLoc = commaLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < acceptorPathVecNum; tmp++)
		{	
			int commaLoc = acceptorPathMaxAnchorSizeVecStr.find(",", startLoc);
			string tmpPathMaxAnchorSizeStr = acceptorPathMaxAnchorSizeVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathMaxAnchorSizeInt = atoi(tmpPathMaxAnchorSizeStr.c_str());
			acceptorPathMaxAnchorSizeVec.push_back(tmpPathMaxAnchorSizeInt);
			startLoc = commaLoc + 1;
		}

		//cout << "start to insertAndInitiateNewSJwithMultiPathVecFromJuncFile ..." << endl;

		this->insertAndInitiateNewSJwithMultiPathVecFromJuncFile(
			juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, juncSupportNum,
			donerPathVec, acceptorPathVec, donerPathVecNum, acceptorPathVecNum, 
			pathSupportNumVec_doner, pathSupportNumVec_acceptor, indexInfo,
			donerPathMaxAnchorSizeVec, acceptorPathMaxAnchorSizeVec,
			juncDonerAnchorSizeMax, juncAcceptorAnchorSizeMax,
			flankString, flankStringCase);
	}

	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize(
		const string& juncStr, Index_Info* indexInfo)
	{
		//cout << "getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash starts ..." << endl;
		//cout << "juncStr: " << juncStr << endl;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				cout << "error in getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize" << endl;
				exit(1);
			}
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		int nextTabLoc = juncStr.find("\t", startLoc);
		if(nextTabLoc == string::npos)
			juncFieldVec.push_back(juncStr.substr(startLoc));
		else
			juncFieldVec.push_back(juncStr.substr(startLoc, nextTabLoc - startLoc));

		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		string juncNameStr = juncFieldVec[3];
		string juncSupportNumStr = juncFieldVec[4];
		int juncSupportNum = atoi(juncSupportNumStr.c_str());
		string juncDonerAnchorSizeStr = juncFieldVec[5];
		int juncDonerAnchorSize = atoi(juncDonerAnchorSizeStr.c_str());
		string juncAcceptorAnchorSizeStr = juncFieldVec[6];
		int juncAcceptorAnchorSize = atoi(juncAcceptorAnchorSizeStr.c_str());
		int tmpJuncSize = juncAcceptorStartPos - juncDonerEndPos - 1;
		this->insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum_anchorSize(
			juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, 
			indexInfo, juncSupportNum, juncDonerAnchorSize, juncAcceptorAnchorSize);
	}

	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_backSpliceOnly(
		const string& juncStr, Index_Info* indexInfo)
	{
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				cout << "error in getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_backSpliceOnly" << endl;
				exit(1);
			}
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));

		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		string juncNameStr = juncFieldVec[3];
		string juncSupportNumStr = juncFieldVec[4];
		int juncSupportNum = atoi(juncSupportNumStr.c_str());
		//int juncDonerAnchorSizeMax = 30;
		//int juncAcceptorAnchorSizeMax = 30;
		int tmpJuncSize = juncAcceptorStartPos - juncDonerEndPos - 1;
		if(tmpJuncSize < 0)
		{
			this->insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum(
				juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, indexInfo, juncSupportNum);					
		}		
	}

	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_includeBackSplice(
		const string& juncStr, Index_Info* indexInfo)
	{
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				cout << "error in getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_backSpliceOnly" << endl;
				exit(1);
			}
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		int nextTabLoc = juncStr.find("\t", startLoc);
		if(nextTabLoc == string::npos)
			juncFieldVec.push_back(juncStr.substr(startLoc));
		else
			juncFieldVec.push_back(juncStr.substr(startLoc, nextTabLoc-startLoc));

		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		string juncNameStr = juncFieldVec[3];
		string juncSupportNumStr = juncFieldVec[4];
		int juncSupportNum = atoi(juncSupportNumStr.c_str());
		//int juncDonerAnchorSizeMax = 30;
		//int juncAcceptorAnchorSizeMax = 30;
		int tmpJuncSize = juncAcceptorStartPos - juncDonerEndPos - 1;
		//if(tmpJuncSize < 0)
		//{
			this->insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum(
				juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, indexInfo, juncSupportNum);					
		//}		
	}

	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum(
		const string& juncStr, Index_Info* indexInfo)
	{
		//cout << "getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash starts ..." << endl;
		//cout << "juncStr: " << juncStr << endl;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				cout << "error in getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum" << endl;
				exit(1);
			}
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		int nextTabLoc = juncStr.find("\t", startLoc);
		if(nextTabLoc == string::npos)
			juncFieldVec.push_back(juncStr.substr(startLoc));
		else
			juncFieldVec.push_back(juncStr.substr(startLoc, nextTabLoc-startLoc));
		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		string juncNameStr = juncFieldVec[3];
		string juncSupportNumStr = juncFieldVec[4];
		int juncSupportNum = atoi(juncSupportNumStr.c_str());
		//int juncDonerAnchorSizeMax = 30;
		//int juncAcceptorAnchorSizeMax = 30;
		int tmpJuncSize = juncAcceptorStartPos - juncDonerEndPos - 1;
		#ifdef DETECT_CIRCULAR_RNA
		this->insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum(
			juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, 
			indexInfo, juncSupportNum 
			//juncDonerAnchorSizeMax, juncAcceptorAnchorSizeMax
			);		
		#else
		if(tmpJuncSize >= 50)
		{
			this->insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum(
					juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, 
					indexInfo, juncSupportNum);
		}
		#endif
	}

	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_flankStringChange(
		const string& juncStr, Index_Info* indexInfo)
	{
		//cout << "getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash starts ..." << endl;
		//cout << "juncStr: " << juncStr << endl;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				cout << "error in getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_flankStringChange" << endl;
				exit(1);
			}
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		int nextTabLoc = juncStr.find("\t", startLoc);
		if(nextTabLoc == string::npos)
			juncFieldVec.push_back(juncStr.substr(startLoc));
		else
			juncFieldVec.push_back(juncStr.substr(startLoc, nextTabLoc-startLoc));
		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		//string juncNameStr = juncFieldVec[3];
		string juncFlankStringChange = juncFieldVec[3];
		string juncSupportNumStr = juncFieldVec[4];
		int juncSupportNum = atoi(juncSupportNumStr.c_str());
		//int juncDonerAnchorSizeMax = 30;
		//int juncAcceptorAnchorSizeMax = 30;
		int tmpJuncSize = juncAcceptorStartPos - juncDonerEndPos - 1;
		// #ifdef DETECT_CIRCULAR_RNA
		// this->insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum(
		// 	juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, 
		// 	indexInfo, juncSupportNum 
		// 	//juncDonerAnchorSizeMax, juncAcceptorAnchorSizeMax
		// 	);		
		// #else
		if(tmpJuncSize >= 50)
		{
			this->insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum_flankStringChange(
					juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, 
					indexInfo, juncSupportNum, juncFlankStringChange);
		}
		//#endif
	}

	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_onlyChrNamePos(
		const string& juncStr, Index_Info* indexInfo)
	{
		//cout << "getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash starts ..." << endl;
		//cout << "juncStr: " << juncStr << endl;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));

		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		string juncNameStr = juncFieldVec[3];
		int juncDonerAnchorSizeMax = 30;
		int juncAcceptorAnchorSizeMax = 30;
		int tmpJuncSize = juncAcceptorStartPos - juncDonerEndPos - 1;
		#ifdef DETECT_CIRCULAR_RNA
		this->insertAndInitiateNewSJfromJuncFile_onlyChrNamePos(
			juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, indexInfo,
			juncDonerAnchorSizeMax, juncAcceptorAnchorSizeMax);		
		#else
		if(tmpJuncSize >= 50)
		{
			this->insertAndInitiateNewSJfromJuncFile_onlyChrNamePos(
				juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, indexInfo,
				juncDonerAnchorSizeMax, juncAcceptorAnchorSizeMax);
		}
		#endif
	}

	void insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, Index_Info* indexInfo,
		int tmpJuncSupportNum)
	{
		AlignInfer_Info tmpAlignInferInfo;
		tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum(
			mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
			tmpJuncSupportNum, indexInfo);

		AlignInferJunctionMap::iterator tmpAlignInferJunctionMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if(tmpAlignInferJunctionMapIter == alignInferJunctionMapVec[mapChrNameInt].end()) // tmpSJposDonerPos not found
		{	
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
		}
		else // tmpSJposDonerPos found
		{
			AcceptorStartPos2AlignInferInfoMap::iterator tmpAcceptorStartPos2AlignInferInfoMapIter
				= (tmpAlignInferJunctionMapIter->second).find(tmpSJposAcceptorStart);
			if(tmpAcceptorStartPos2AlignInferInfoMapIter == (tmpAlignInferJunctionMapIter->second).end()) // tmpSJposAcceptorStart not found
			{
				(tmpAlignInferJunctionMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			}
			else
			{
				cout << "both tmpSJposDonerEnd and tmpSJposAcceptorStart found in hash: chrNameInt: "
					<< mapChrNameInt << "  tmpDonerEndPos: " << tmpSJposDonerEnd 
					<< "  tmpAcceptorStartPos: " << tmpSJposAcceptorStart  << endl;
			}
		}

		alignInferInfoVec.push_back(tmpAlignInferInfo);
		currentAlignInferInfoVecIndex ++;
	}

	void insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum_flankStringChange(
		int mapChrNameInt, int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		Index_Info* indexInfo, int tmpJuncSupportNum, string& tmpJuncFlankStringChange)
	{
		AlignInfer_Info tmpAlignInferInfo;
		tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_flankStringChange(
			mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
			tmpJuncSupportNum, indexInfo, tmpJuncFlankStringChange);

		AlignInferJunctionMap::iterator tmpAlignInferJunctionMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if(tmpAlignInferJunctionMapIter == alignInferJunctionMapVec[mapChrNameInt].end()) // tmpSJposDonerPos not found
		{	
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
		}
		else // tmpSJposDonerPos found
		{
			AcceptorStartPos2AlignInferInfoMap::iterator tmpAcceptorStartPos2AlignInferInfoMapIter
				= (tmpAlignInferJunctionMapIter->second).find(tmpSJposAcceptorStart);
			if(tmpAcceptorStartPos2AlignInferInfoMapIter == (tmpAlignInferJunctionMapIter->second).end()) // tmpSJposAcceptorStart not found
			{
				(tmpAlignInferJunctionMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			}
			else
			{
				cout << "both tmpSJposDonerEnd and tmpSJposAcceptorStart found in hash: chrNameInt: "
					<< mapChrNameInt << "  tmpDonerEndPos: " << tmpSJposDonerEnd 
					<< "  tmpAcceptorStartPos: " << tmpSJposAcceptorStart  << endl;
			}
		}

		alignInferInfoVec.push_back(tmpAlignInferInfo);
		currentAlignInferInfoVecIndex ++;
	}	

	void updateEncompassingNum(int index)
	{
		alignInferInfoVec[index].updateEncompassingNum();
	}

	void insertAndInitiateNewSJfromJuncFile_chrNamePos_supportNum_anchorSize(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, Index_Info* indexInfo,
		int tmpJuncSupportNum, int tmpJuncDonerAnchorSizeMax, int tmpJuncAcceptorAnchorSizeMax)
	{
		AlignInfer_Info tmpAlignInferInfo;
		tmpAlignInferInfo.initiateAlignInferInfo_chrNamePos_supportNum_anchorSize(
			mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
			tmpJuncSupportNum, tmpJuncDonerAnchorSizeMax, tmpJuncAcceptorAnchorSizeMax, indexInfo);

		AlignInferJunctionMap::iterator tmpAlignInferJunctionMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if(tmpAlignInferJunctionMapIter == alignInferJunctionMapVec[mapChrNameInt].end()) // tmpSJposDonerPos not found
		{	
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
		}
		else // tmpSJposDonerPos found
		{
			AcceptorStartPos2AlignInferInfoMap::iterator tmpAcceptorStartPos2AlignInferInfoMapIter
				= (tmpAlignInferJunctionMapIter->second).find(tmpSJposAcceptorStart);
			if(tmpAcceptorStartPos2AlignInferInfoMapIter == (tmpAlignInferJunctionMapIter->second).end()) // tmpSJposAcceptorStart not found
			{
				(tmpAlignInferJunctionMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			}
			else
			{
				// cout << "both tmpSJposDonerEnd and tmpSJposAcceptorStart found in hash: chrNameInt: "
				// 	<< mapChrNameInt << "  tmpDonerEndPos: " << tmpSJposDonerEnd 
				// 	<< "  tmpAcceptorStartPos: " << tmpSJposAcceptorStart  << endl;
			}
		}

		alignInferInfoVec.push_back(tmpAlignInferInfo);
		currentAlignInferInfoVecIndex ++;
	}

	void insertAndInitiateNewSJfromJuncFile_onlyChrNamePos(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, Index_Info* indexInfo,
		int tmpJuncDonerAnchorSizeMax, int tmpJuncAcceptorAnchorSizeMax)
	{
		AlignInfer_Info tmpAlignInferInfo;
		tmpAlignInferInfo.initiateAlignInferInfo_onlyChrNamePos(
			mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, indexInfo,
			tmpJuncDonerAnchorSizeMax, tmpJuncAcceptorAnchorSizeMax);

		AlignInferJunctionMap::iterator tmpAlignInferJunctionMapIter 
			= alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if(tmpAlignInferJunctionMapIter == alignInferJunctionMapVec[mapChrNameInt].end()) // tmpSJposDonerPos not found
		{	
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
		}
		else // tmpSJposDonerPos found
		{
			AcceptorStartPos2AlignInferInfoMap::iterator tmpAcceptorStartPos2AlignInferInfoMapIter
				= (tmpAlignInferJunctionMapIter->second).find(tmpSJposAcceptorStart);
			if(tmpAcceptorStartPos2AlignInferInfoMapIter == (tmpAlignInferJunctionMapIter->second).end()) // tmpSJposAcceptorStart not found
			{
				(tmpAlignInferJunctionMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			}
			else
			{
				cout << "both tmpSJposDonerEnd and tmpSJposAcceptorStart found in hash: chrNameInt: "
					<< mapChrNameInt << "  tmpDonerEndPos: " << tmpSJposDonerEnd 
					<< "  tmpAcceptorStartPos: " << tmpSJposAcceptorStart  << endl;
			}
		}

		alignInferInfoVec.push_back(tmpAlignInferInfo);
		currentAlignInferInfoVecIndex ++;
	}

	void getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash_withAlterSpliceSiteAnchorSimilarity(
		const string& juncStr, Index_Info* indexInfo)
	{
		//cout << "getAlignInferInfoFromJunc_InsertIntoAlignInferJunctionHash starts ..." << endl;
		//cout << "juncStr: " << juncStr << endl;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 20; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			//cout << "tmpJuncField: " << tmpJuncField << endl;
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));

		string juncChrNameStr = juncFieldVec[0];
		int juncChrNameInt = indexInfo->convertStringToInt(juncChrNameStr);
		string juncDonerEndPosStr = juncFieldVec[1];
		int juncDonerEndPos = atoi(juncDonerEndPosStr.c_str());
		string juncAcceptorStartPosStr = juncFieldVec[2];
		int juncAcceptorStartPos = atoi(juncAcceptorStartPosStr.c_str());
		string juncSupportNumStr = juncFieldVec[3];
		int juncSupportNum = atoi(juncSupportNumStr.c_str());
		//cout << "juncDonerEndPos: " << juncDonerEndPos << endl;
		//cout << "juncAcceptorStartPos: " << juncAcceptorStartPos << endl;
		string flankString = juncFieldVec[4];
		string flankStringCaseStr = juncFieldVec[5];
		int flankStringCase = atoi(flankStringCaseStr.c_str());

		string juncNameStr = juncFieldVec[6];

		string juncDonerAnchorSizeMaxStr = juncFieldVec[7];
		string juncAcceptorAnchorSizeMaxStr = juncFieldVec[8];
		int juncDonerAnchorSizeMax = atoi(juncDonerAnchorSizeMaxStr.c_str());
		int juncAcceptorAnchorSizeMax = atoi(juncAcceptorAnchorSizeMaxStr.c_str());

		string juncDonerJumpCodeVecStr = juncFieldVec[9];
		string juncAcceptorJumpCodeVecStr = juncFieldVec[10];
		string donerPathVecStr = juncFieldVec[11];
		string acceptorPathVecStr = juncFieldVec[12]; 
		int donerPathVecNum = atoi(donerPathVecStr.c_str());
		int acceptorPathVecNum = atoi(acceptorPathVecStr.c_str());
		string donerPathSupportNumVecStr = juncFieldVec[13];
		string acceptorPathSupportNumVecStr = juncFieldVec[14]; 
		string donerPathMaxAnchorSizeVecStr = juncFieldVec[15];
		//cout << "donerPathMaxAnchorSizeVecStr: " << donerPathMaxAnchorSizeVecStr << endl;
		string acceptorPathMaxAnchorSizeVecStr = juncFieldVec[16];//juncFieldVec[16];
		//cout << "acceptorPathMaxAnchorSizeVecStr: " << acceptorPathMaxAnchorSizeVecStr << endl;
		string extensionSimilarityStr_doner = juncFieldVec[17];
		string extensionSimilarityStr_acceptor = juncFieldVec[18];
		string alterSpliceAnchorSimilarityStr_doner = juncFieldVec[19];
		string alterSpliceAnchorSimilarityStr_acceptor = juncStr.substr(startLoc);//juncFieldVec[14];

		vector< vector<Jump_Code> > donerPathVec;
		vector< vector<Jump_Code> > acceptorPathVec;
		startLoc = 0;
		for(int tmp = 0; tmp < donerPathVecNum; tmp++)
		{	
			int commaLoc = juncDonerJumpCodeVecStr.find(",", startLoc);
			string tmpJumpCodeVecField = juncDonerJumpCodeVecStr.substr(startLoc, commaLoc-startLoc);
			vector<Jump_Code> tmpCigarStringJumpCodeVec_backward;
			this->cigarString2jumpCodeVec(tmpJumpCodeVecField, tmpCigarStringJumpCodeVec_backward);
			donerPathVec.push_back(tmpCigarStringJumpCodeVec_backward);
			startLoc = commaLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < acceptorPathVecNum; tmp++)
		{	
			int commaLoc = juncAcceptorJumpCodeVecStr.find(",", startLoc);
			string tmpJumpCodeVecField = juncAcceptorJumpCodeVecStr.substr(startLoc, commaLoc-startLoc);
			vector<Jump_Code> tmpCigarStringJumpCodeVec_forward;
			this->cigarString2jumpCodeVec(tmpJumpCodeVecField, tmpCigarStringJumpCodeVec_forward);
			acceptorPathVec.push_back(tmpCigarStringJumpCodeVec_forward);
			startLoc = commaLoc + 1;
		}

		vector<int> pathSupportNumVec_doner;
		vector<int> pathSupportNumVec_acceptor;		
		startLoc = 0;
		for(int tmp = 0; tmp < donerPathVecNum; tmp++)
		{	
			int commaLoc = donerPathSupportNumVecStr.find(",", startLoc);
			string tmpPathSupportNumStr = donerPathSupportNumVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathSupportNumInt = atoi(tmpPathSupportNumStr.c_str());
			pathSupportNumVec_doner.push_back(tmpPathSupportNumInt);
			startLoc = commaLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < acceptorPathVecNum; tmp++)
		{	
			int commaLoc = acceptorPathSupportNumVecStr.find(",", startLoc);
			string tmpPathSupportNumStr = acceptorPathSupportNumVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathSupportNumInt = atoi(tmpPathSupportNumStr.c_str());
			pathSupportNumVec_acceptor.push_back(tmpPathSupportNumInt);
			startLoc = commaLoc + 1;
		}

		vector<int> donerPathMaxAnchorSizeVec;
		vector<int> acceptorPathMaxAnchorSizeVec;
		startLoc = 0;
		for(int tmp = 0; tmp < donerPathVecNum; tmp++)
		{	
			int commaLoc = donerPathMaxAnchorSizeVecStr.find(",", startLoc);
			string tmpPathMaxAnchorSizeStr = donerPathMaxAnchorSizeVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathMaxAnchorSizeInt = atoi(tmpPathMaxAnchorSizeStr.c_str());
			donerPathMaxAnchorSizeVec.push_back(tmpPathMaxAnchorSizeInt);
			startLoc = commaLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < acceptorPathVecNum; tmp++)
		{	
			int commaLoc = acceptorPathMaxAnchorSizeVecStr.find(",", startLoc);
			string tmpPathMaxAnchorSizeStr = acceptorPathMaxAnchorSizeVecStr.substr(startLoc, commaLoc-startLoc);
			int tmpPathMaxAnchorSizeInt = atoi(tmpPathMaxAnchorSizeStr.c_str());
			acceptorPathMaxAnchorSizeVec.push_back(tmpPathMaxAnchorSizeInt);
			startLoc = commaLoc + 1;
		}

		vector<Jump_Code> tmpExtensionJumpCodeVec_doner;
		int tmpExtensionPenalty_doner;
		vector<Jump_Code> tmpExtensionJumpCodeVec_acceptor;
		int tmpExtensionPenalty_acceptor;
		vector< pair<int,int> > tmpAlterSpliceSitePairVec_doner;
		vector< vector<Jump_Code> > tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner;
		vector< int > tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner;
		vector< pair<int,int> > tmpAlterSpliceSitePairVec_acceptor;
		vector< vector<Jump_Code> > tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor;
		vector< int > tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor;

		this->generate_extension_alterSpliceSite_anchorSimilarity(
				extensionSimilarityStr_doner, extensionSimilarityStr_acceptor,
				alterSpliceAnchorSimilarityStr_doner, alterSpliceAnchorSimilarityStr_acceptor,
				tmpExtensionJumpCodeVec_doner, tmpExtensionPenalty_doner,
				tmpExtensionJumpCodeVec_acceptor, tmpExtensionPenalty_acceptor,
				tmpAlterSpliceSitePairVec_doner,
				tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner,
				tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner,
				tmpAlterSpliceSitePairVec_acceptor,
				tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor,
				tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor);
		// cout << "tmpExtensionJumpCodeVec_doner.size(): " << tmpExtensionJumpCodeVec_doner.size() << endl;
		// cout << "tmpExtensionJumpCodeVec_acceptor.size(): " << tmpExtensionJumpCodeVec_acceptor.size() << endl;
		// cout << "tmpExtensionPenalty_doner: " << tmpExtensionPenalty_doner << endl;
		// cout << "tmpExtensionPenalty_acceptor: " << tmpExtensionPenalty_acceptor << endl;
		//cout << "start to insertAndInitiateNewSJwithMultiPathVecFromJuncFile ..." << endl;
		this->insertAndInitiateNewSJwithMultiPathVecFromJuncFile_withAlterSpliceSiteAnchorSimilarity(
				juncChrNameInt, juncDonerEndPos, juncAcceptorStartPos, juncSupportNum,
				donerPathVec, acceptorPathVec, donerPathVecNum, acceptorPathVecNum, 
				pathSupportNumVec_doner, pathSupportNumVec_acceptor, indexInfo,
				donerPathMaxAnchorSizeVec, acceptorPathMaxAnchorSizeVec,
				juncDonerAnchorSizeMax, juncAcceptorAnchorSizeMax,
				flankString, flankStringCase,
				tmpExtensionJumpCodeVec_doner, tmpExtensionPenalty_doner,
				tmpExtensionJumpCodeVec_acceptor, tmpExtensionPenalty_acceptor,
				tmpAlterSpliceSitePairVec_doner,
				tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner,
				tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner,
				tmpAlterSpliceSitePairVec_acceptor,
				tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor,
				tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor);
	}

	void generate_extension_alterSpliceSite_anchorSimilarity(
		string& extensionSimilarityStr_doner,
		string& extensionSimilarityStr_acceptor,
		string& alterSpliceAnchorSimilarityStr_doner,
		string& alterSpliceAnchorSimilarityStr_acceptor,
		vector<Jump_Code>& tmpExtensionJumpCodeVec_doner,
		int& tmpExtensionPenalty_doner,
		vector<Jump_Code>& tmpExtensionJumpCodeVec_acceptor,
		int& tmpExtensionPenalty_acceptor,
		vector< pair<int,int> >& tmpAlterSpliceSitePairVec_doner,
		vector< vector<Jump_Code> >& tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner,
		vector< int >& tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner,
		vector< pair<int,int> >& tmpAlterSpliceSitePairVec_acceptor,
		vector< vector<Jump_Code> >& tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor,
		vector< int >& tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor)
	{
		//cout << endl << "generate_extension_alterSpliceSite_anchorSimilarity starts ...." << endl;
		// extension-doner
		//cout << "extensionSimilarityStr_doner: " << extensionSimilarityStr_doner << endl;
		int tmpColonLoc = extensionSimilarityStr_doner.find(":");
		string tmpJumpCodeVecStr = extensionSimilarityStr_doner.substr(0, tmpColonLoc);
		//cout << "tmpJumpCodeVecStr_doner: " << tmpJumpCodeVecStr << endl;
		this->cigarString2jumpCodeVec(tmpJumpCodeVecStr, tmpExtensionJumpCodeVec_doner);
		string tmpPenaltyStr = extensionSimilarityStr_doner.substr(tmpColonLoc+1);
		tmpExtensionPenalty_doner = atoi(tmpPenaltyStr.c_str());
		//cout << "tmpExtensionPenalty_doner: " << tmpExtensionPenalty_doner << endl;
		// extension-acceptor
		//cout << "extensionSimilarityStr_acceptor: " << extensionSimilarityStr_acceptor << endl;
		tmpColonLoc = extensionSimilarityStr_acceptor.find(":");
		tmpJumpCodeVecStr = extensionSimilarityStr_acceptor.substr(0, tmpColonLoc);
		//cout << "tmpJumpCodeVecStr_acceptor: " << tmpJumpCodeVecStr << endl;
		this->cigarString2jumpCodeVec(tmpJumpCodeVecStr, tmpExtensionJumpCodeVec_acceptor);
		tmpPenaltyStr = extensionSimilarityStr_acceptor.substr(tmpColonLoc+1);
		tmpExtensionPenalty_acceptor = atoi(tmpPenaltyStr.c_str());
		//cout << "tmpExtensionPenalty_acceptor: " << tmpExtensionPenalty_acceptor << endl;
		//exit(1);
		// alterSpliceSite-doner
		int startLoc = 0;
		while(1)
		{
			int tmpCommaLoc = alterSpliceAnchorSimilarityStr_doner.find(",", startLoc);
			if(tmpCommaLoc == string::npos)
				break;
			string tmpAlterSpliceSiteFieldStr = alterSpliceAnchorSimilarityStr_doner.substr(
				startLoc, tmpCommaLoc-startLoc);
			int tmpCommaLoc_1 = tmpAlterSpliceSiteFieldStr.find(":", 0);
			int tmpCommaLoc_2 = tmpAlterSpliceSiteFieldStr.find(":", tmpCommaLoc_1+1);
			int tmpCommaLoc_3 = tmpAlterSpliceSiteFieldStr.find(":", tmpCommaLoc_2+1);
			string tmpSJdonerEndPosStr = tmpAlterSpliceSiteFieldStr.substr(0,tmpCommaLoc_1);
			int tmpSJdonerEndPosInt = atoi(tmpSJdonerEndPosStr.c_str());
			string tmpSJacceptorStartPosStr = tmpAlterSpliceSiteFieldStr.substr(tmpCommaLoc_1+1, tmpCommaLoc_2-tmpCommaLoc_1-1);
			int tmpSJacceptorStartPosInt = atoi(tmpSJacceptorStartPosStr.c_str());
			string tmpSJanchorSimilarityJumpCodeVecStr = tmpAlterSpliceSiteFieldStr.substr(tmpCommaLoc_2+1, tmpCommaLoc_3-tmpCommaLoc_2-1);
			vector<Jump_Code> tmpSJanchorSimilarityJumpCodeVec;
			this->cigarString2jumpCodeVec(tmpSJanchorSimilarityJumpCodeVecStr, tmpSJanchorSimilarityJumpCodeVec);
			string tmpSJanchorSimilarityPenaltyStr = tmpAlterSpliceSiteFieldStr.substr(tmpCommaLoc_3+1);
			int tmpSJanchorSimilarityPenaltyInt = atoi(tmpSJanchorSimilarityPenaltyStr.c_str());
			startLoc = tmpCommaLoc + 1;

			tmpAlterSpliceSitePairVec_doner.push_back(pair<int,int>(tmpSJdonerEndPosInt, tmpSJacceptorStartPosInt));
			tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner.push_back(tmpSJanchorSimilarityJumpCodeVec);
			tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner.push_back(tmpSJanchorSimilarityPenaltyInt);
		}

		// alterSpliceSite-acceptor
		startLoc = 0;
		while(1)
		{
			int tmpCommaLoc = alterSpliceAnchorSimilarityStr_acceptor.find(",", startLoc);
			if(tmpCommaLoc == string::npos)
				break;
			string tmpAlterSpliceSiteFieldStr = alterSpliceAnchorSimilarityStr_acceptor.substr(
				startLoc, tmpCommaLoc-startLoc);
			int tmpCommaLoc_1 = tmpAlterSpliceSiteFieldStr.find(":", 0);
			int tmpCommaLoc_2 = tmpAlterSpliceSiteFieldStr.find(":", tmpCommaLoc_1+1);
			int tmpCommaLoc_3 = tmpAlterSpliceSiteFieldStr.find(":", tmpCommaLoc_2+1);
			string tmpSJdonerEndPosStr = tmpAlterSpliceSiteFieldStr.substr(0,tmpCommaLoc_1);
			int tmpSJdonerEndPosInt = atoi(tmpSJdonerEndPosStr.c_str());
			string tmpSJacceptorStartPosStr = tmpAlterSpliceSiteFieldStr.substr(tmpCommaLoc_1+1, tmpCommaLoc_2-tmpCommaLoc_1-1);
			int tmpSJacceptorStartPosInt = atoi(tmpSJacceptorStartPosStr.c_str());
			string tmpSJanchorSimilarityJumpCodeVecStr = tmpAlterSpliceSiteFieldStr.substr(tmpCommaLoc_2+1, tmpCommaLoc_3-tmpCommaLoc_2-1);
			vector<Jump_Code> tmpSJanchorSimilarityJumpCodeVec;
			this->cigarString2jumpCodeVec(tmpSJanchorSimilarityJumpCodeVecStr, tmpSJanchorSimilarityJumpCodeVec);
			string tmpSJanchorSimilarityPenaltyStr = tmpAlterSpliceSiteFieldStr.substr(tmpCommaLoc_3+1);
			int tmpSJanchorSimilarityPenaltyInt = atoi(tmpSJanchorSimilarityPenaltyStr.c_str());
			startLoc = tmpCommaLoc + 1;

			tmpAlterSpliceSitePairVec_acceptor.push_back(pair<int,int>(tmpSJdonerEndPosInt, tmpSJacceptorStartPosInt));
			tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor.push_back(tmpSJanchorSimilarityJumpCodeVec);
			tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor.push_back(tmpSJanchorSimilarityPenaltyInt);
		}
	}

	void insertAndInitiateNewSJwithMultiPathVecFromJuncFile(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, int tmpSJsupportNum,
		vector< vector<Jump_Code> >& tmpSJdonerPathVec, 
		vector< vector<Jump_Code> >& tmpSJacceptorPathVec, 
		int tmpSJdonerPathVecNum, int tmpSJacceptorPathVecNum, 
		vector<int>& tmpSJdonerPathSupportNumVec,
		vector<int>& tmpSJacceptorPathSupportNumVec,
		Index_Info* indexInfo,
		vector<int>& tmpDonerPathMaxAnchorSizeVec,
		vector<int>& tmpAcceptorPathMaxAnchorSizeVec,
		int tmpJuncDonerAnchorSizeMax, int tmpJuncAcceptorAnchorSizeMax,
		string& flankString, int flankStringCase)
	{
		AlignInfer_Info tmpAlignInferInfo;
		tmpAlignInferInfo.initiateAlignInferInfoWithMultiPathVec(
			mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
			tmpSJsupportNum, tmpSJdonerPathVec, tmpSJacceptorPathVec,
			tmpSJdonerPathVecNum, tmpSJacceptorPathVecNum,
			tmpSJdonerPathSupportNumVec, tmpSJacceptorPathSupportNumVec,
			tmpDonerPathMaxAnchorSizeVec, tmpAcceptorPathMaxAnchorSizeVec,
			tmpJuncDonerAnchorSizeMax, tmpJuncAcceptorAnchorSizeMax, indexInfo,
			flankString, flankStringCase);
		//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
		// AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
		// tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
		// 	tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

		// alignInferJunctionMapVec[mapChrNameInt].insert(
		// 	pair< int, AcceptorStartPos2AlignInferInfoMap> (
		// 		tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
		AlignInferJunctionMap::iterator tmpAlignInferJunctionMapIter = alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if(tmpAlignInferJunctionMapIter == alignInferJunctionMapVec[mapChrNameInt].end()) // tmpSJposDonerPos not found
		{	
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
		}
		else // tmpSJposDonerPos found
		{
			AcceptorStartPos2AlignInferInfoMap::iterator tmpAcceptorStartPos2AlignInferInfoMapIter
				= (tmpAlignInferJunctionMapIter->second).find(tmpSJposAcceptorStart);
			if(tmpAcceptorStartPos2AlignInferInfoMapIter == (tmpAlignInferJunctionMapIter->second).end()) // tmpSJposAcceptorStart not found
			{
				(tmpAlignInferJunctionMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			}
			else
			{
				cout << "both tmpSJposDonerEnd and tmpSJposAcceptorStart found in hash: chrNameInt: "
					<< mapChrNameInt << "  tmpDonerEndPos: " << tmpSJposDonerEnd 
					<< "  tmpAcceptorStartPos: " << tmpSJposAcceptorStart  << endl;
			}
		}

		alignInferInfoVec.push_back(tmpAlignInferInfo);
		currentAlignInferInfoVecIndex ++;
	}

	void insertAndInitiateNewSJwithMultiPathVecFromJuncFile_withAlterSpliceSiteAnchorSimilarity(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, int tmpSJsupportNum,
		vector< vector<Jump_Code> >& tmpSJdonerPathVec, 
		vector< vector<Jump_Code> >& tmpSJacceptorPathVec, 
		int tmpSJdonerPathVecNum, int tmpSJacceptorPathVecNum, 
		vector<int>& tmpSJdonerPathSupportNumVec,
		vector<int>& tmpSJacceptorPathSupportNumVec,
		Index_Info* indexInfo,
		vector<int>& tmpDonerPathMaxAnchorSizeVec,
		vector<int>& tmpAcceptorPathMaxAnchorSizeVec,
		int tmpJuncDonerAnchorSizeMax, int tmpJuncAcceptorAnchorSizeMax,
		string& flankString, int flankStringCase,
		vector<Jump_Code>& tmpExtensionJumpCodeVec_doner,
		int& tmpExtensionPenalty_doner,
		vector<Jump_Code>& tmpExtensionJumpCodeVec_acceptor,
		int& tmpExtensionPenalty_acceptor,
		vector< pair<int,int> >& tmpAlterSpliceSitePairVec_doner,
		vector< vector<Jump_Code> >& tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner,
		vector< int >& tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner,
		vector< pair<int,int> >& tmpAlterSpliceSitePairVec_acceptor,
		vector< vector<Jump_Code> >& tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor,
		vector< int >& tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor)
	{
		AlignInfer_Info tmpAlignInferInfo;
		tmpAlignInferInfo.initiateAlignInferInfoWithMultiPathVec_withAlterSpliceSiteAnchorSimilarity(
			mapChrNameInt, tmpSJposDonerEnd, tmpSJposAcceptorStart, 
			tmpSJsupportNum, tmpSJdonerPathVec, tmpSJacceptorPathVec,
			tmpSJdonerPathVecNum, tmpSJacceptorPathVecNum,
			tmpSJdonerPathSupportNumVec, tmpSJacceptorPathSupportNumVec,
			tmpDonerPathMaxAnchorSizeVec, tmpAcceptorPathMaxAnchorSizeVec,
			tmpJuncDonerAnchorSizeMax, tmpJuncAcceptorAnchorSizeMax, indexInfo,
			flankString, flankStringCase,
			tmpExtensionJumpCodeVec_doner, tmpExtensionPenalty_doner,
			tmpExtensionJumpCodeVec_acceptor, tmpExtensionPenalty_acceptor,
			tmpAlterSpliceSitePairVec_doner,
			tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner,
			tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner,
			tmpAlterSpliceSitePairVec_acceptor,
			tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor,
			tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor);
		//cout << "start to insert: " << mapChrNameInt << " " << tmpSJposDonerEnd << " " << tmpSJposAcceptorStart << endl;
		AlignInferJunctionMap::iterator tmpAlignInferJunctionMapIter = alignInferJunctionMapVec[mapChrNameInt].find(tmpSJposDonerEnd);
		if(tmpAlignInferJunctionMapIter == alignInferJunctionMapVec[mapChrNameInt].end()) // tmpSJposDonerPos not found
		{	
			AcceptorStartPos2AlignInferInfoMap tmpAcceptorStartPos2AlignInferInfoMap;
			tmpAcceptorStartPos2AlignInferInfoMap.insert(pair<int, int>(
				tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));

			alignInferJunctionMapVec[mapChrNameInt].insert(
				pair< int, AcceptorStartPos2AlignInferInfoMap> (
					tmpSJposDonerEnd, tmpAcceptorStartPos2AlignInferInfoMap));
		}
		else // tmpSJposDonerPos found
		{
			AcceptorStartPos2AlignInferInfoMap::iterator tmpAcceptorStartPos2AlignInferInfoMapIter
				= (tmpAlignInferJunctionMapIter->second).find(tmpSJposAcceptorStart);
			if(tmpAcceptorStartPos2AlignInferInfoMapIter == (tmpAlignInferJunctionMapIter->second).end()) // tmpSJposAcceptorStart not found
			{
				(tmpAlignInferJunctionMapIter->second).insert(pair<int, int>(
					tmpSJposAcceptorStart, currentAlignInferInfoVecIndex));
			}
			else
			{
				cout << "both tmpSJposDonerEnd and tmpSJposAcceptorStart found in hash: chrNameInt: "
					<< mapChrNameInt << "  tmpDonerEndPos: " << tmpSJposDonerEnd 
					<< "  tmpAcceptorStartPos: " << tmpSJposAcceptorStart  << endl;
			}
		}

		alignInferInfoVec.push_back(tmpAlignInferInfo);
		currentAlignInferInfoVecIndex ++;
	}


	/*
	void generateSJjumpCodeVec(
			vector<Jump_Code>& cigarStringJumpCodeVec,
			vector< int >& tmpSJindexVec_cigarStringJumpCodeVec,
			vector< vector<Jump_Code> >& SJjumpCodeVecVec_backward,
			vector< vector<Jump_Code> >& SJjumpCodeVecVec_forward)
	{
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		//cout << "cigarStringJumpCodeVecSize: " << cigarStringJumpCodeVecSize << endl;
		for(int tmp = 0; tmp < tmpSJindexVec_cigarStringJumpCodeVec.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int tmpSJindex_cigarStringJumpCodeVec = tmpSJindexVec_cigarStringJumpCodeVec[tmp];
			//cout << "tmpSJindex_cigarStringJumpCodeVec " << tmpSJindex_cigarStringJumpCodeVec << endl;
			vector<Jump_Code> tmpSJjumpCodeVec_backward;
			vector<Jump_Code> tmpSJjumpCodeVec_forward;
			for(int tmpJumpCodeIndex = tmpSJindex_cigarStringJumpCodeVec-1; tmpJumpCodeIndex >= 0; tmpJumpCodeIndex--)
			{
				//cout << "getJumpCode: " << cigarStringJumpCodeVec[tmpJumpCodeIndex].toString() << endl;
				tmpSJjumpCodeVec_backward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
			}
			for(int tmpJumpCodeIndex = tmpSJindex_cigarStringJumpCodeVec+1; tmpJumpCodeIndex < cigarStringJumpCodeVecSize; tmpJumpCodeIndex++)
			{
				//cout << "getJumpCode: " << cigarStringJumpCodeVec[tmpJumpCodeIndex].toString() << endl;
				tmpSJjumpCodeVec_forward.push_back(cigarStringJumpCodeVec[tmpJumpCodeIndex]);
			}
			SJjumpCodeVecVec_backward.push_back(tmpSJjumpCodeVec_backward);
			SJjumpCodeVecVec_forward.push_back(tmpSJjumpCodeVec_forward);
		}
	}
	*/

	void generateSJposVecFromJumpCodeVec(int mapChrPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
		vector< pair<int,int> >& SJposPairVec, vector<int>& SJindexVec_cigarStringJumpCodeVec)
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int lastJumpCodeIndex = tmp-1;
				int currentJumpCodeIndex = tmp;
				int tmpDonerEndPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, lastJumpCodeIndex);
				int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, currentJumpCodeIndex) + 1;
				SJposPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
				SJindexVec_cigarStringJumpCodeVec.push_back(tmp);
			}
		}		
	}

	void generateSJindexVecFromJumpCodeVec(int mapChrPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
		vector<int>& SJindexVec_cigarStringJumpCodeVec)
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int lastJumpCodeIndex = tmp-1;
				int currentJumpCodeIndex = tmp;
				//int tmpDonerEndPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, lastJumpCodeIndex);
				//int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, currentJumpCodeIndex) + 1;
				//SJposPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
				SJindexVec_cigarStringJumpCodeVec.push_back(tmp);
			}
		}		
	}		

	int getEndLocInReadOfSpecificJumpCode(
		vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		if(jumpCodeIndex < 0)
			return 0;
		int tmpEndLocInRead = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "D")
			{
			}
			else if(tmpJumpCodeType == "N")
			{
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return tmpEndLocInRead;
	}	

	int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		int tmpEndPos = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
	}

	void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
	{

		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
		
		string candidateJumpCodeType = "SMNIDX";
		while(1)
		{
			jumpCodeEndPosInCigarStr = 
				jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
			if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
				{break;}
			else
			{
				tmpJumpCodeLength = 
					atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	string returnAlignInferInfoHashInfoStr(Index_Info* indexInfo)
	{
		//cout << "start to gather Hash info:" << endl;
		string tmpStr;
		//int tmpJuncNum = 0;
		
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize(indexInfo, tmp+1);
			tmpStr = tmpStr + tmpAlignInferInfoStr;
		}
		//cout << "end of gathering Hash info:" << endl;
		return tmpStr;	
	}

	void compareExtensionAnchorSimilarity(Index_Info* indexInfo)
	{
		int alignInferInfoVecSize = alignInferInfoVec.size();
		#pragma omp parallel for schedule(dynamic)		
		for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
		{
			alignInferInfoVec[tmp].checkExtension_compareAnchorSimilarity(indexInfo);
		}			
	}

	void getAlterSpliceSites_compareAnchorSimilarity(SJhash_Info* SJhashInfo, int offset, Index_Info* indexInfo)
	{
		int alignInferInfoVecSize = alignInferInfoVec.size();
		#pragma omp parallel for schedule(dynamic)		
		for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
		{
			alignInferInfoVec[tmp].getAlterSpliceSites_compareAnchorSimilarity(SJhashInfo, offset, indexInfo);
		}		
	}

	void getAlterSpliceSites_compareAnchorSimilarity_onlyLowSupportSJ(SJhash_Info* SJhashInfo, int offset, Index_Info* indexInfo)
	{
		int alignInferInfoVecSize = alignInferInfoVec.size();
		#pragma omp parallel for schedule(dynamic)
		for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
		{
			alignInferInfoVec[tmp].getAlterSpliceSites_compareAnchorSimilarity_onlyLowSupportSJ(SJhashInfo, offset, indexInfo);
		}		
	}

	void getAlterSpliceSites_compareAnchorSimilarity_allSJ_withDefaultAnchorSize(
		SJhash_Info* SJhashInfo, int offset, int defaultAnchorSize, Index_Info* indexInfo)
	{
		int alignInferInfoVecSize = alignInferInfoVec.size();
		#pragma omp parallel for schedule(dynamic)
		for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
		{
			alignInferInfoVec[tmp].getAlterSpliceSites_compareAnchorSimilarity_allSJ_withDefaultAnchorSize(
				SJhashInfo, offset, defaultAnchorSize, indexInfo);
		}		
	}

	void convert2SJhashInfo(SJhash_Info* SJhashInfo, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			int tmpChrNameInt = alignInferInfoVec[tmp].returnChrNameInt();
			int tmpDonerEndPos = alignInferInfoVec[tmp].returnDonerEndPos();
			int tmpAcceptorStartPos = alignInferInfoVec[tmp].returnAcceptorStartPos();
			//cout << "tmpChrNameInt: " << tmpChrNameInt << endl;
			//cout << "tmpDonerEndPos: " << tmpDonerEndPos << endl;
			//cout << "tmpAcceptorStartPos: " << tmpAcceptorStartPos << endl;
			SJhashInfo->insert2AreaAndStringHash(tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos, indexInfo);
		}		
	}

	void outputExtensionSpliceSiteAnchorSimilarity(const string& tmpSJseqSimilarityStats)
	{
		ofstream outputSJseqSimilarityStats_ofs(tmpSJseqSimilarityStats.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			int tmpDonerAnchorSize = alignInferInfoVec[tmp].returnAnchorSizeMax_doner();
			int tmpAcceptorAnchorSize = alignInferInfoVec[tmp].returnAnchorSizeMax_acceptor();
			int tmpDonerExtensionEditDistance = alignInferInfoVec[tmp].returnExtension_penalty_doner();
			int tmpAcceptorExtensionEditDistance = alignInferInfoVec[tmp].returnExtension_penalty_acceptor();
			double relativeDonerExtensionEditDistanceRatio 
				= (double)tmpDonerExtensionEditDistance / (double)tmpDonerAnchorSize;
			double relativeAcceptorExtensionEditDistanceRatio
				= (double)tmpAcceptorExtensionEditDistance / (double)tmpAcceptorAnchorSize;
			double extensionEditDistanceRatio_max;
			double extensionEditDistanceRatio_min;
			if(relativeDonerExtensionEditDistanceRatio > relativeAcceptorExtensionEditDistanceRatio)
			{
				extensionEditDistanceRatio_max = relativeDonerExtensionEditDistanceRatio;
				extensionEditDistanceRatio_min = relativeAcceptorExtensionEditDistanceRatio;
			}
			else
			{
				extensionEditDistanceRatio_max = relativeAcceptorExtensionEditDistanceRatio;
				extensionEditDistanceRatio_min = relativeDonerExtensionEditDistanceRatio;
			}
			int extensionEditDistance_max = (int)(extensionEditDistanceRatio_max * 30);
			int extensionEditDistance_min = (int)(extensionEditDistanceRatio_min * 30);
			outputSJseqSimilarityStats_ofs << relativeDonerExtensionEditDistanceRatio << "\t"
				<< relativeAcceptorExtensionEditDistanceRatio << "\t"
				<< extensionEditDistance_max << "\t"
				<< extensionEditDistance_min << "\t " << endl;
		}
		outputSJseqSimilarityStats_ofs.close();
	}

	void outputAlignInferInfoHashInfo(Index_Info* indexInfo, 
		const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		//int tmpJuncNum = 0;
		
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr 
				= alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize(
					indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_chrNamePosOnly(Index_Info* indexInfo, 
		const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr 
				= alignInferInfoVec[tmp].returnAlignInferInfoStr_chrNamePosOnly(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_chrNamePos_strand(Index_Info* indexInfo, 
		const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr 
				= alignInferInfoVec[tmp].returnAlignInferInfoStr_chrNamePos_strand(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_chrNamePos_supportNum(Index_Info* indexInfo, 
		const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr 
				= alignInferInfoVec[tmp].returnAlignInferInfoStr_chrNamePos_supportNum(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_chrNamePos_supportNum_encompassingNum(Index_Info* indexInfo, 
		const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr 
				= alignInferInfoVec[tmp].returnAlignInferInfoStr_chrNamePos_supportNum_encompassingNum(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_chrNamePos_supportNum_anchorSize(Index_Info* indexInfo, 
		const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr 
				= alignInferInfoVec[tmp].returnAlignInferInfoStr_chrNamePos_supportNum_anchorSize(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_chrNamePos_supportNum_anchorSize_XM(Index_Info* indexInfo, 
		const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr 
				= alignInferInfoVec[tmp].returnAlignInferInfoStr_chrNamePos_supportNum_anchorSize_XM(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_onlyExtensionAnchorSimilarity(Index_Info* indexInfo, const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		//int tmpJuncNum = 0;
		
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_onlyExtensionAnchorSimilarity(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(
		Index_Info* indexInfo, const string& outputAlignInferHashFile)
	{
		ofstream outputAlignInferHash_ofs(outputAlignInferHashFile.c_str());
		//int tmpJuncNum = 0;
		
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
			outputAlignInferHash_ofs << tmpAlignInferInfoStr << endl;
		}
		outputAlignInferHash_ofs.close();
	}

	void outputAlignInferInfoHashInfo_classified_forRemapping(Index_Info* indexInfo,
		const string& outputSJpath_all_classified, const string& outputSJpath_kept, 
		const string& outputSJpath_filterOut, const string& outputSJpath_filterOut_canBeExtendedDirectly, 
		const string& outputSJpath_filterOut_canMove2alterSpliceSite,
		ofstream& log_ofs)
	{
		int totalSJnum = alignInferInfoVec.size();
		int keptSJnum = 0;
		int filterOutSJnum = 0;
		int filterOutSJnum_canBeExtended = 0;
		int filterOutSJnum_canMove2alterSpliceSite = 0;

		ofstream allClassifiedSJ_ofs(outputSJpath_all_classified.c_str());
		ofstream keptSJ_ofs(outputSJpath_kept.c_str());
		ofstream filterOutSJ_ofs(outputSJpath_filterOut.c_str());
		ofstream filterOutSJ_canBeExtendedDirectly_ofs(outputSJpath_filterOut_canBeExtendedDirectly.c_str());
		ofstream filterOutSJ_canMove2alterSpliceSite_ofs(outputSJpath_filterOut_canMove2alterSpliceSite.c_str());		
		
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			int tmpAlignInferInfo_chrNameInt = alignInferInfoVec[tmp].returnChrNameInt();
			int tmpAlignInferInfo_donerEndPosInt = alignInferInfoVec[tmp].returnDonerEndPos();
			int tmpAlignInferInfo_acceptorStartPosInt = alignInferInfoVec[tmp].returnAcceptorStartPos();
			int tmpAlignInferInfo_supportNumInt = alignInferInfoVec[tmp].returnSupportNum();
			string tmpAlignInferInfo_chrNameStr = indexInfo->returnChrNameStr(tmpAlignInferInfo_chrNameInt);
			string tmpAlignInferInfo_donerPosStr = int_to_str(tmpAlignInferInfo_donerEndPosInt);
			string tmpAlignInferInfo_acceptorPosStr = int_to_str(tmpAlignInferInfo_acceptorStartPosInt);
			string tmpAlignInferInfo_supportNumStr = int_to_str(tmpAlignInferInfo_supportNumInt);
			string tmpAlignInferInfoStr = tmpAlignInferInfo_chrNameStr + "\t"
				+ tmpAlignInferInfo_donerPosStr + "\t" + tmpAlignInferInfo_acceptorPosStr + "\tJUNC_" 
				+ int_to_str(tmp+1) + "\t" + tmpAlignInferInfo_chrNameStr;
			bool filterOut_extended_bool = alignInferInfoVec[tmp].canBeExtended_bool();
			bool filterOut_move2alterSpliceSite_bool = alignInferInfoVec[tmp].exactTheSame2someAlterSpliceSite();
			if(filterOut_extended_bool)
			{
				allClassifiedSJ_ofs << tmpAlignInferInfoStr << "\tF" << endl;
				filterOutSJ_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJ_canBeExtendedDirectly_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJnum ++;
				filterOutSJnum_canBeExtended ++;
			}
			else if(filterOut_move2alterSpliceSite_bool)
			{
				allClassifiedSJ_ofs << tmpAlignInferInfoStr << "\tF" << endl;
				filterOutSJ_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJ_canMove2alterSpliceSite_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJnum ++;
				filterOutSJnum_canMove2alterSpliceSite ++;
			}
			else
			{
				allClassifiedSJ_ofs << tmpAlignInferInfoStr << "\tT" << endl;
				keptSJ_ofs << tmpAlignInferInfoStr << endl;
				keptSJnum ++;
			}
		}

		log_ofs << "totalSJnum: " << totalSJnum << endl;
		log_ofs << "keptSJnum: " << keptSJnum << endl;
		log_ofs << "filterOutSJnum: " << filterOutSJnum << endl;
		log_ofs << "\tfilterOutSJnum_canBeExtended: " << filterOutSJnum_canBeExtended << endl;
		log_ofs << "\tfilterOutSJnum_canMove2alterSpliceSite: " << filterOutSJnum_canMove2alterSpliceSite << endl; 

		allClassifiedSJ_ofs.close();
		keptSJ_ofs.close();
		filterOutSJ_ofs.close();
		filterOutSJ_canBeExtendedDirectly_ofs.close();
		filterOutSJ_canMove2alterSpliceSite_ofs.close();
	}

	void outputAlignInferInfoHashInfo_classified_finalClassification(Index_Info* indexInfo,
		const string& outputSJpath_all_classified, const string& outputSJpath_kept, 
		const string& outputSJpath_filterOut, const string& outputSJpath_filterOut_canBeExtendedDirectly, 
		const string& outputSJpath_filterOut_canMove2alterSpliceSite, 
		const string& outputSJpath_filterOut_lowSupNumNoncanonical,
		int lowSupportThreshold, ofstream& log_ofs)
	{
		int totalSJnum = alignInferInfoVec.size();
		int keptSJnum = 0;
		int filterOutSJnum = 0;
		int filterOutSJnum_canBeExtended = 0;
		int filterOutSJnum_canMove2alterSpliceSite = 0;
		int filterOutSJnum_lowSupNumNoncanonical = 0;

		ofstream allClassifiedSJ_ofs(outputSJpath_all_classified.c_str());
		ofstream keptSJ_ofs(outputSJpath_kept.c_str());
		ofstream filterOutSJ_ofs(outputSJpath_filterOut.c_str());
		ofstream filterOutSJ_canBeExtendedDirectly_ofs(outputSJpath_filterOut_canBeExtendedDirectly.c_str());
		ofstream filterOutSJ_canMove2alterSpliceSite_ofs(outputSJpath_filterOut_canMove2alterSpliceSite.c_str());
		ofstream filterOutSJ_lowSupNumNoncanonical_ofs(outputSJpath_filterOut_lowSupNumNoncanonical.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			int tmpAlignInferInfo_chrNameInt = alignInferInfoVec[tmp].returnChrNameInt();
			int tmpAlignInferInfo_donerEndPosInt = alignInferInfoVec[tmp].returnDonerEndPos();
			int tmpAlignInferInfo_acceptorStartPosInt = alignInferInfoVec[tmp].returnAcceptorStartPos();
			int tmpAlignInferInfo_supportNumInt = alignInferInfoVec[tmp].returnSupportNum();
			string tmpAlignInferInfo_chrNameStr = indexInfo->returnChrNameStr(tmpAlignInferInfo_chrNameInt);
			string tmpAlignInferInfo_donerPosStr = int_to_str(tmpAlignInferInfo_donerEndPosInt);
			string tmpAlignInferInfo_acceptorPosStr = int_to_str(tmpAlignInferInfo_acceptorStartPosInt);
			string tmpAlignInferInfo_supportNumStr = int_to_str(tmpAlignInferInfo_supportNumInt);
			string tmpAlignInferInfoStr = tmpAlignInferInfo_chrNameStr + "\t"
				+ tmpAlignInferInfo_donerPosStr + "\t" + tmpAlignInferInfo_acceptorPosStr + "\tJUNC_" 
				+ int_to_str(tmp+1) + "\t" + tmpAlignInferInfo_chrNameStr;
			bool filterOut_extended_bool 
				= alignInferInfoVec[tmp].canBeExtended_bool();
			bool filterOut_move2alterSpliceSite_bool 
				= alignInferInfoVec[tmp].exactTheSame2someAlterSpliceSite();
			bool filterOut_lowSupNumNonCanonical_bool 
				= alignInferInfoVec[tmp].lowSupportNumNonCanonical_bool(lowSupportThreshold);
			if(filterOut_extended_bool)
			{
				allClassifiedSJ_ofs << tmpAlignInferInfoStr << "\tF" << endl;
				filterOutSJ_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJ_canBeExtendedDirectly_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJnum ++;
				filterOutSJnum_canBeExtended ++;
			}
			else if(filterOut_move2alterSpliceSite_bool)
			{
				allClassifiedSJ_ofs << tmpAlignInferInfoStr << "\tF" << endl;
				filterOutSJ_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJ_canMove2alterSpliceSite_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJnum ++;
				filterOutSJnum_canMove2alterSpliceSite ++;
			}
			else if(filterOut_lowSupNumNonCanonical_bool)
			{
				allClassifiedSJ_ofs << tmpAlignInferInfoStr << "\tF" << endl;
				filterOutSJ_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJ_lowSupNumNoncanonical_ofs << tmpAlignInferInfoStr << endl;
				filterOutSJnum ++;
				filterOutSJnum_lowSupNumNoncanonical ++;
			}
			else
			{
				allClassifiedSJ_ofs << tmpAlignInferInfoStr << "\tT" << endl;
				keptSJ_ofs << tmpAlignInferInfoStr << endl;
				keptSJnum ++;
			} 
		}

		log_ofs << "totalSJnum: " << totalSJnum << endl;
		log_ofs << "keptSJnum: " << keptSJnum << endl;
		log_ofs << "filterOutSJnum: " << filterOutSJnum << endl;
		log_ofs << "\tfilterOutSJnum_canBeExtended: " << filterOutSJnum_canBeExtended << endl;
		log_ofs << "\tfilterOutSJnum_canMove2alterSpliceSite: " << filterOutSJnum_canMove2alterSpliceSite << endl; 
		log_ofs << "\tfilterOutSJnum_lowSupNumNoncanonical: " << filterOutSJnum_lowSupNumNoncanonical << endl;

		allClassifiedSJ_ofs.close();
		keptSJ_ofs.close();
		filterOutSJ_ofs.close();
		filterOutSJ_canBeExtendedDirectly_ofs.close();
		filterOutSJ_canMove2alterSpliceSite_ofs.close();
		filterOutSJ_lowSupNumNoncanonical_ofs.close();
	}	

	void outputSJ_extended(const string& outputSJpath_extended, Index_Info* indexInfo)
	{
		ofstream outputSJpath_extended_ofs(outputSJpath_extended.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if(alignInferInfoVec[tmp].canBeExtended_bool())
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_extended_ofs << tmpAlignInferInfoStr << endl;			
			}
		}
		outputSJpath_extended_ofs.close();
	}

	void outputSJ_nonExtended_nonAlterSpliceSite(const string& outputSJpath_nonExtended_nonAlterSpliceSite, Index_Info* indexInfo)
	{
		ofstream outputSJpath_nonExtended_nonAlterSpliceSite_ofs(outputSJpath_nonExtended_nonAlterSpliceSite.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if((!(alignInferInfoVec[tmp].canBeExtended_bool()))&&(!(alignInferInfoVec[tmp].similar2someAlterSpliceSite())))
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_nonExtended_nonAlterSpliceSite_ofs << tmpAlignInferInfoStr << endl;			
			}
		}	
		outputSJpath_nonExtended_nonAlterSpliceSite_ofs.close();
	}

	void outputSJ_nonExtended_withAlterSpliceSite(const string& outputSJpath_nonExtended_withAlterSpliceSite, Index_Info* indexInfo)
	{
		ofstream outputSJpath_nonExtended_withAlterSpliceSite_ofs(outputSJpath_nonExtended_withAlterSpliceSite.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if((!(alignInferInfoVec[tmp].canBeExtended_bool()))&&(alignInferInfoVec[tmp].similar2someAlterSpliceSite()))
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_nonExtended_withAlterSpliceSite_ofs << tmpAlignInferInfoStr << endl;			
			}
		}
		outputSJpath_nonExtended_withAlterSpliceSite_ofs.close();
	}

	void outputSJ_nonExtended_withAlterSpliceSite_valid(
		const string& outputSJpath_nonExtended_withAlterSpliceSite_valid,
		Index_Info* indexInfo)
	{
		ofstream outputSJpath_nonExtended_withAlterSpliceSite_valid_ofs(outputSJpath_nonExtended_withAlterSpliceSite_valid.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if((!(alignInferInfoVec[tmp].canBeExtended_bool()))&&(alignInferInfoVec[tmp].similar2someAlterSpliceSite())
				&&(alignInferInfoVec[tmp].validSJ()))
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_nonExtended_withAlterSpliceSite_valid_ofs << tmpAlignInferInfoStr << endl;			
			}
		}
		outputSJpath_nonExtended_withAlterSpliceSite_valid_ofs.close();
	}

	void outputSJ_nonExtended_withAlterSpliceSite_invalid(
		const string& outputSJpath_nonExtended_withAlterSpliceSite_invalid,
		Index_Info* indexInfo)
	{
		ofstream outputSJpath_nonExtended_withAlterSpliceSite_invalid_ofs(outputSJpath_nonExtended_withAlterSpliceSite_invalid.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if((!(alignInferInfoVec[tmp].canBeExtended_bool()))&&(alignInferInfoVec[tmp].similar2someAlterSpliceSite())
				&&(!alignInferInfoVec[tmp].validSJ()))
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_nonExtended_withAlterSpliceSite_invalid_ofs << tmpAlignInferInfoStr << endl;			
			}
		}
		outputSJpath_nonExtended_withAlterSpliceSite_invalid_ofs.close();
	}

	void outputSJ_nonExtended_multiAnchorExactlyTheSame(
		const string& outputSJpath_nonExtended_multiAnchorExactlyTheSame, Index_Info* indexInfo)
	{
		ofstream outputSJpath_nonExtended_multiAnchorExactlyTheSame_ofs(outputSJpath_nonExtended_multiAnchorExactlyTheSame.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if((!(alignInferInfoVec[tmp].canBeExtended_bool()))&&(alignInferInfoVec[tmp].exactTheSame2someAlterSpliceSite()))
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_nonExtended_multiAnchorExactlyTheSame_ofs << tmpAlignInferInfoStr << endl;			
			}
		}
		outputSJpath_nonExtended_multiAnchorExactlyTheSame_ofs.close();
	}

	void outputSJtypeNum(const string& outputSJtypeNumPath)
	{
		ofstream outputSJtypeNumPath_ofs(outputSJtypeNumPath.c_str());
		outputSJtypeNumPath_ofs << "Note:" << endl;
		outputSJtypeNumPath_ofs << "LONG_ANCHOR_LEN: " << LONG_ANCHOR_LEN << endl;
		outputSJtypeNumPath_ofs << "SHORT_ANCHOR_LEN: " << SHORT_ANCHOR_LEN << endl;
		outputSJtypeNumPath_ofs << "longAnchorSJ: both anchors longer or equal to LONG_ANCHOR_LEN" << endl;
		outputSJtypeNumPath_ofs << "midAnchorSJ: short anchor shorter than LONG_ANCHOR_LEN, but longer or equal to LONG_ANCHOR_LEN" << endl; 
		outputSJtypeNumPath_ofs << "shortAnchorSJ: short anchor shorter than SHORT_ANCHOR_LEN" << endl;
		outputSJtypeNumPath_ofs << "SHORT_SJ_DISTANCE_LEN: " << SHORT_SJ_DISTANCE_LEN << endl;
		outputSJtypeNumPath_ofs << "MID_SJ_DISTANCE_LEN: " << MID_SJ_DISTANCE_LEN << endl;
		outputSJtypeNumPath_ofs << "LONG_SJ_DISTANCE_LEN: " << LONG_SJ_DISTANCE_LEN << endl;
		outputSJtypeNumPath_ofs << "short DistanceSJ: SJ size less than SHORT_SJ_DISTANCE_LEN" << endl;
		outputSJtypeNumPath_ofs << "mid DistanceSJ: SJ size less than MID_SJ_DISTANCE_LEN but longer or equal to SHORT_SJ_DISTANCE_LEN" << endl; 
		outputSJtypeNumPath_ofs << "long DistanceSJ: SJ size less than LONG_SJ_DISTANCE_LEN but longer or equal to MID_SJ_DISTANCE_LEN" << endl; 
		outputSJtypeNumPath_ofs << "Super long DistanceSJ: SJ size longer or equal to LONG_SJ_DISTANCE_LEN" << endl;

		int totalSJnum = 0;
		int supNumCount = 5;
		int nonSemiOrCanonicalTypeNum = 3;
		int anchorSizeLevelNum = 3;
		int SJdistanceLevelNum = 4;
		vector< vector< vector<int> > > typeNumVec_anchorSize;
		vector< vector< vector<int> > > typeNumVec_SJdistance;
		for(int tmp = 0; tmp <= supNumCount; tmp++)
		{
			vector< vector<int> > nonSemiOrCanonicalTypeNumVec_anchorSize;
			vector< vector<int> > nonSemiOrCanonicalTypeNumVec_SJdistance;
			for(int tmp2 = 0; tmp2 < nonSemiOrCanonicalTypeNum; tmp2++)
			{
				vector<int> anchorSizeLevelTypeNumVec;
				for(int tmp3 = 0; tmp3 < anchorSizeLevelNum; tmp3++)
				{
					anchorSizeLevelTypeNumVec.push_back(0);
				}
				vector<int> SJdistanceLevelTypeNumVec;
				for(int tmp4 = 0; tmp4 < SJdistanceLevelNum; tmp4++)
				{
					SJdistanceLevelTypeNumVec.push_back(0);
				}
				nonSemiOrCanonicalTypeNumVec_anchorSize.push_back(anchorSizeLevelTypeNumVec);
				nonSemiOrCanonicalTypeNumVec_SJdistance.push_back(SJdistanceLevelTypeNumVec);
			}
			typeNumVec_anchorSize.push_back(nonSemiOrCanonicalTypeNumVec_anchorSize);
			typeNumVec_SJdistance.push_back(nonSemiOrCanonicalTypeNumVec_SJdistance);
		}

		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			int tmpAlignInferSupportNum = alignInferInfoVec[tmp].returnSupportNum();
			bool nonCanonicalBool = alignInferInfoVec[tmp].nonCanonical_bool();
			bool semiCanonicalBool = alignInferInfoVec[tmp].semiCanonical_bool();
			bool canonicalBool = alignInferInfoVec[tmp].canonical_bool();
			int anchorSizeLevel = alignInferInfoVec[tmp].returnAnchorSizeLevel();
			int SJdistanceLevel = alignInferInfoVec[tmp].returnSJdistanceLevel();

			int tmpAlignInferSupNumIndex = tmpAlignInferSupportNum-1;
			if(tmpAlignInferSupNumIndex > 5)
				tmpAlignInferSupNumIndex = 5;
			if(nonCanonicalBool)
			{
				((typeNumVec_anchorSize[tmpAlignInferSupNumIndex])[0])[anchorSizeLevel] ++;
				((typeNumVec_SJdistance[tmpAlignInferSupNumIndex])[0])[SJdistanceLevel] ++;
			}
			else if(semiCanonicalBool)
			{
				((typeNumVec_anchorSize[tmpAlignInferSupNumIndex])[1])[anchorSizeLevel] ++;
				((typeNumVec_SJdistance[tmpAlignInferSupNumIndex])[1])[SJdistanceLevel] ++;
			}
			else
			{
				((typeNumVec_anchorSize[tmpAlignInferSupNumIndex])[2])[anchorSizeLevel] ++;
				((typeNumVec_SJdistance[tmpAlignInferSupNumIndex])[2])[SJdistanceLevel] ++;
			}

			totalSJnum ++;
		}
		outputSJtypeNumPath_ofs << "totalSJnum: " << totalSJnum << endl << endl;
		for(int tmp = 0; tmp < supNumCount; tmp++)
		{
			outputSJtypeNumPath_ofs << "SupportNum: " << tmp + 1 << endl;

			outputSJtypeNumPath_ofs << "\tcanonicalSJ: " << (((typeNumVec_anchorSize[tmp])[2])[2] 
				+ ((typeNumVec_anchorSize[tmp])[2])[1] + ((typeNumVec_anchorSize[tmp])[2])[0]) << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tlongAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[2])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tmidAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[2])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tshortAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[2])[0] << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tshort DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[2])[3] << endl;
			outputSJtypeNumPath_ofs << "\t\tmid DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[2])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tlong DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[2])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tSuper long DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[2])[0] << endl << endl;

			outputSJtypeNumPath_ofs << "\tsemiCanonicalSJ: " << (((typeNumVec_anchorSize[tmp])[1])[2] 
				+ ((typeNumVec_anchorSize[tmp])[1])[1] + ((typeNumVec_anchorSize[tmp])[1])[0]) << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tlongAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[1])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tmidAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[1])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tshortAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[1])[0] << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tshort DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[1])[3] << endl;
			outputSJtypeNumPath_ofs << "\t\tmid DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[1])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tlong DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[1])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tSuper long DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[1])[0] << endl << endl;

			outputSJtypeNumPath_ofs << "\tnonCanonicalSJ: " << (((typeNumVec_anchorSize[tmp])[0])[2] 
				+ ((typeNumVec_anchorSize[tmp])[0])[1] + ((typeNumVec_anchorSize[tmp])[0])[0]) << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tlongAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[0])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tmidAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[0])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tshortAnchorSJ: " << ((typeNumVec_anchorSize[tmp])[0])[0] << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tshort DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[0])[3] << endl;
			outputSJtypeNumPath_ofs << "\t\tmid DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[0])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tlong DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[0])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tSuper long DistanceSJ: " << ((typeNumVec_SJdistance[tmp])[0])[0] << endl << endl;
			outputSJtypeNumPath_ofs << endl;
		}
		outputSJtypeNumPath_ofs << "SupportNum > " << supNumCount << endl;

			outputSJtypeNumPath_ofs << "\tcanonicalSJ: " << (((typeNumVec_anchorSize[supNumCount])[2])[2] 
				+ ((typeNumVec_anchorSize[supNumCount])[2])[1] + ((typeNumVec_anchorSize[supNumCount])[2])[0]) << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tlongAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[2])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tmidAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[2])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tshortAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[2])[0] << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tshort DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[2])[3] << endl;
			outputSJtypeNumPath_ofs << "\t\tmid DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[2])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tlong DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[2])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tSuper long DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[2])[0] << endl << endl;

			outputSJtypeNumPath_ofs << "\tsemiCanonicalSJ: " << (((typeNumVec_anchorSize[supNumCount])[1])[2] 
				+ ((typeNumVec_anchorSize[supNumCount])[1])[1] + ((typeNumVec_anchorSize[supNumCount])[1])[0]) << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tlongAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[1])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tmidAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[1])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tshortAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[1])[0] << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tshort DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[1])[3] << endl;
			outputSJtypeNumPath_ofs << "\t\tmid DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[1])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tlong DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[1])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tSuper long DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[1])[0] << endl << endl;
			
			outputSJtypeNumPath_ofs << "\tnonCanonicalSJ: " << (((typeNumVec_anchorSize[supNumCount])[0])[2] 
				+ ((typeNumVec_anchorSize[supNumCount])[0])[1] + ((typeNumVec_anchorSize[supNumCount])[0])[0]) << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tlongAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[0])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tmidAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[0])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tshortAnchorSJ: " << ((typeNumVec_anchorSize[supNumCount])[0])[0] << endl << endl;
			outputSJtypeNumPath_ofs << "\t\tshort DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[0])[3] << endl;
			outputSJtypeNumPath_ofs << "\t\tmid DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[0])[2] << endl;
			outputSJtypeNumPath_ofs << "\t\tlong DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[0])[1] << endl;
			outputSJtypeNumPath_ofs << "\t\tSuper long DistanceSJ: " << ((typeNumVec_SJdistance[supNumCount])[0])[0] << endl << endl;

		outputSJtypeNumPath_ofs.close();
	}

	void outputSJ_kept(const string& outputSJpath_kept, Index_Info* indexInfo)
	{
		ofstream outputSJpath_kept_ofs(outputSJpath_kept.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if(alignInferInfoVec[tmp].validSJ())
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_kept_ofs << tmpAlignInferInfoStr << endl;			
			}
		}
		outputSJpath_kept_ofs.close();
	}

	void outputSJ_filterOut(const string& outputSJpath_filterOut, 
		const string& outputSJpath_filterOut_extension,
		const string& outputSJpath_filterOut_alterSpliceSite,
		const string& outputSJpath_filterOut_lowSupportNonCanonical,
		Index_Info* indexInfo)
	{
		ofstream outputSJpath_filterOut_ofs(
			outputSJpath_filterOut.c_str());
		ofstream outputSJpath_filterOut_extension_ofs(
			outputSJpath_filterOut_extension.c_str());
		ofstream outputSJpath_filterOut_alterSpliceSite_ofs(
			outputSJpath_filterOut_alterSpliceSite.c_str());
		ofstream outputSJpath_filterOut_lowSupportNonCanonical_ofs(
			outputSJpath_filterOut_lowSupportNonCanonical.c_str());		
		int lowSupNum_max = 2;
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if(!(alignInferInfoVec[tmp].validSJ()))
			{
				string tmpAlignInferInfoStr = alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(indexInfo, tmp+1);
				outputSJpath_filterOut_ofs << tmpAlignInferInfoStr << endl;
				bool canBeExtendedBool 
					= alignInferInfoVec[tmp].canBeExtended_bool();
				bool exactTheSame2someAlterSpliceSiteBool 
					= alignInferInfoVec[tmp].exactTheSame2someAlterSpliceSite();
				bool lowSupportNonCanonicalBool
					= alignInferInfoVec[tmp].lowSupportNumNonCanonical_bool(lowSupNum_max);
				if(canBeExtendedBool)
				{
					outputSJpath_filterOut_extension_ofs
						<< tmpAlignInferInfoStr << endl;
				}
				else if(exactTheSame2someAlterSpliceSiteBool)
				{
					outputSJpath_filterOut_alterSpliceSite_ofs
						<< tmpAlignInferInfoStr << endl;
				}
				else if(lowSupportNonCanonicalBool)
				{
					outputSJpath_filterOut_lowSupportNonCanonical_ofs
						<< tmpAlignInferInfoStr << endl;
				}
				else
				{}
			}
		}
		outputSJpath_filterOut_ofs.close();
		outputSJpath_filterOut_extension_ofs.close();
		outputSJpath_filterOut_alterSpliceSite_ofs.close();
		outputSJpath_filterOut_lowSupportNonCanonical_ofs.close();
	}

	bool canBeRefinedByKeptSJspliceSiteWithOffset(int SJindex, int& refinedSJindex)
	{
		int offset = 10;
		int tmpSJchrNameInt = alignInferInfoVec[SJindex].returnChrNameInt();
		int tmpSJdonerEndPos = alignInferInfoVec[SJindex].returnDonerEndPos();
		int tmpSJacceptorStartPos = alignInferInfoVec[SJindex].returnAcceptorStartPos();
		for(int tmpSJdonerEndPos_offset = tmpSJdonerEndPos - offset; 
			tmpSJdonerEndPos_offset <= tmpSJdonerEndPos + offset; 
			tmpSJdonerEndPos_offset ++)
		{
			for(int tmpSJacceptorStartPos_offset = tmpSJacceptorStartPos - offset; 
				tmpSJacceptorStartPos_offset <= tmpSJacceptorStartPos + offset; 
				tmpSJacceptorStartPos_offset ++)
			{
				if((tmpSJdonerEndPos_offset == tmpSJdonerEndPos)&&(tmpSJacceptorStartPos_offset == tmpSJacceptorStartPos))
					continue;
				int tmpFoundSJindex_offset = searchAndReturnAlignInferInfoVecIndex(tmpSJchrNameInt,
					tmpSJdonerEndPos_offset, tmpSJacceptorStartPos_offset);
				if(tmpFoundSJindex_offset < 0)
					continue;
				else
				{
					bool tmpFoundSJvalid_offset_bool = alignInferInfoVec[tmpFoundSJindex_offset].validSJ();
					tmpFoundSJvalid_offset_bool = true;
					if(tmpFoundSJvalid_offset_bool)
					{
						refinedSJindex = tmpFoundSJindex_offset;
						return true;
					}
				}
			}
		}
		return false;
	}

	void outputSJ_filterOut_canBeRefinedWithKeptSJspliceSite(
		string& outputSJpath_filterOut_lowSupportNonCanonical_canBeRefinedWithKeptSJspliceSite,
		string& outputSJpath_filterOut_lowSupportNonCanonical_canNotBeRefinedWithKeptSJspliceSite,
		Index_Info* indexInfo)
	{
		int lowSupNum_max = 2;
		ofstream output_ofs_1(outputSJpath_filterOut_lowSupportNonCanonical_canBeRefinedWithKeptSJspliceSite.c_str());
		ofstream output_ofs_2(outputSJpath_filterOut_lowSupportNonCanonical_canNotBeRefinedWithKeptSJspliceSite.c_str());
		for(int tmp = 0; tmp < alignInferInfoVec.size(); tmp++)
		{
			if(!(alignInferInfoVec[tmp].validSJ()))
			{
				bool canBeExtendedBool 
					= alignInferInfoVec[tmp].canBeExtended_bool();
				bool exactTheSame2someAlterSpliceSiteBool 
					= alignInferInfoVec[tmp].similar2someAlterSpliceSite();
				bool lowSupportNonCanonicalBool
					= alignInferInfoVec[tmp].lowSupportNumNonCanonical_bool(lowSupNum_max);
				if(canBeExtendedBool || exactTheSame2someAlterSpliceSiteBool)
					continue;
				if(lowSupportNonCanonicalBool)
				{
					int tmpRefinedSJindex;
					bool canBeRefinedByKeptSJspliceSiteWithOffset_bool 
						= this->canBeRefinedByKeptSJspliceSiteWithOffset(tmp, tmpRefinedSJindex);
					string tmpAlignInferInfoStr 
						= alignInferInfoVec[tmp].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(
							indexInfo, tmp);
					if(canBeRefinedByKeptSJspliceSiteWithOffset_bool)
					{	
						// string tmpRefinedSJalignInferInfoStr
						// 	= alignInferInfoVec[tmpRefinedSJindex].returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(
						// 		indexInfo, tmpRefinedSJindex);
						output_ofs_1 << tmpAlignInferInfoStr << endl;
						// output_ofs << "Can be refined by: " << endl << tmpRefinedSJalignInferInfoStr << endl;
					}
					else
					{
						output_ofs_2 << tmpAlignInferInfoStr << endl;
					}
				}
				else
				{}
			}
		}
		output_ofs_1.close();
		output_ofs_2.close();
	}

};





#endif