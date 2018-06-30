// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef REFINEFUSIONBREAKPOINT_INFO_H
#define REFINEFUSIONBREAKPOINT_INFO_H

using namespace std;

class RefineFusionBreakPoint_Info
{
private:
	// input original info
	int oriIncompletePairAlignInfo_chrNameInt;
	bool leftReadHeadUnfixed_bool;
	int oriIncompletePairAlignInfo_chrPos_left;
	vector<Jump_Code> oriIncompletePairAlignInfo_jumpCodeVec_left;
	bool rightReadTailUnfixed_bool;
	int oriIncompletePairAlignInfo_chrNamePos_right;
	vector<Jump_Code> oriIncompletePairAlignInfo_jumpCodeVec_right;

	// output fixed fusion info
	bool leftReadHeadUnfixed_fixed_bool;
	bool leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool;
	int leftHeadFusion_chrNameInt;
	int leftHeadFusion_chrPos;
	vector<Jump_Code> leftHeadFusion_jumpCodeVec;
	int keptAlignInfo_chrPos_left;
	vector<Jump_Code> keptAlignInfo_jumpCodeVec_left;

	bool rightReadTailUnfixed_fixed_bool;
	bool rightReadTailUnfixed_fixed_Nor_or_Rcm_bool;
	int rightTailFusion_chrNameInt;
	int rightTailFusion_chrPos;
	vector<Jump_Code> rightTailFusion_jumpCodeVec;
	int keptAlignInfo_chrPos_right;
	vector<Jump_Code> keptAlignInfo_jumpCodeVec_right;
public:	
	RefineFusionBreakPoint_Info()
	{}

	int countFixFusionBuffer(int bufferMax, int matchLen_1, int matchLen_2)
	{
		if(matchLen_1 <= matchLen_2)
		{
			if(matchLen_1 <= bufferMax)
				return matchLen_1;
			else
				return bufferMax;
		}
		else // matchLen_2 < matchLen_1
		{
			if(matchLen_2 <= bufferMax)
				return matchLen_2;
			else
				return bufferMax;
		}
	}

	int compareTwoStringToGetMismatchNum(string& str_1, string& str_2)
	{
		int tmpMismatchNum = 0;
		int tmpStrLen = str_1.length();
		for(int tmp = 0; tmp < tmpStrLen; tmp++)
		{
			string tmpBase_1 = str_1.substr(tmp, 1);
			string tmpBase_2 = str_2.substr(tmp, 1);
			if(tmpBase_1 != tmpBase_2)	
				tmpMismatchNum ++;
		}
		return tmpMismatchNum;
	}

	void fixDoubleAnchor_fusionJunc_noIns(string& toFixSeqInRead, 
		int toFixChrSeq_chrNameInt_left, int toFixChrSeq_chrNameInt_right,
		int toFixChrSeq_chrPos_left, int toFixChrSeq_chrPos_right, 
		vector<int>& mismatchNumVec, vector<string>& flankStringVec,		
		Index_Info* indexInfo)
	{
		//int insertedSeqLenMax = InsertedBaseNumMaxBetweenFusedTranscript;
		int candiFusChrSeq_length = toFixSeqInRead.size() + 2;
		cout << "candiFusChrSeq_length: " << candiFusChrSeq_length << endl;
		string candiFusChrSeq_left = indexInfo->returnChromStrSubstr(
			toFixChrSeq_chrNameInt_left, toFixChrSeq_chrPos_left, candiFusChrSeq_length);
		cout << "candiFusChrSeq_left: " << candiFusChrSeq_left << endl;
		string candiFusChrSeq_right = indexInfo->returnChromStrSubstr(
			toFixChrSeq_chrNameInt_right, toFixChrSeq_chrPos_right - candiFusChrSeq_length + 1, 
			candiFusChrSeq_length);
		cout << "candiFusChrSeq_right: " << candiFusChrSeq_right << endl;
		int toFixSeqInReadLen = toFixSeqInRead.length();

		// fix canonical fusion junc without inserted bases
		for(int tmpToFixSeqAssignedToLeft = 1; tmpToFixSeqAssignedToLeft < toFixSeqInReadLen;
			tmpToFixSeqAssignedToLeft ++)
		{
			int tmpToFixSeqAssignedToRight = toFixSeqInReadLen - tmpToFixSeqAssignedToLeft;
			string toFixFusionChrSeq_left_tmp
				= candiFusChrSeq_left.substr(0, tmpToFixSeqAssignedToLeft);
			cout << "toFixFusionChrSeq_left_tmp: " << toFixFusionChrSeq_left_tmp << endl;
			string toFixFusionChrSeq_right_tmp 
				= candiFusChrSeq_right.substr(candiFusChrSeq_length - tmpToFixSeqAssignedToRight, 
					tmpToFixSeqAssignedToRight);
			cout << "toFixFusionChrSeq_right_tmp: " << toFixFusionChrSeq_right_tmp << endl;
			string toFixFusionChrSeq_tmp = toFixFusionChrSeq_left_tmp + toFixFusionChrSeq_right_tmp;
			string tmpFlankString_left = candiFusChrSeq_left.substr(tmpToFixSeqAssignedToLeft, 2);
			string tmpFlankString_right = candiFusChrSeq_right.substr(
				candiFusChrSeq_length - tmpToFixSeqAssignedToRight - 2, 2);
			string tmpFlankString = tmpFlankString_left + tmpFlankString_right;
			flankStringVec.push_back(tmpFlankString);
			int tmpMismatchNum = this->compareTwoStringToGetMismatchNum(
				toFixSeqInRead, toFixFusionChrSeq_tmp);
			mismatchNumVec.push_back(tmpMismatchNum);
			cout << "tmpToFixSeqAssignedToLeft: " << tmpToFixSeqAssignedToLeft << endl;
			cout << "tmpToFixSeqAssignedToRight: " << tmpToFixSeqAssignedToRight << endl;
			cout << "tmpMismatchNum: " << tmpMismatchNum << endl;
			cout << "tmpFlankString: " << tmpFlankString << endl;
		}
	}

	void fixDoubleAnchor_fusionJunc_withIns(int insLen, string& toFixSeqInRead, 
		int toFixChrSeq_chrNameInt_left, int toFixChrSeq_chrNameInt_right,
		int toFixChrSeq_chrPos_left, int toFixChrSeq_chrPos_right, 
		vector<int>& mismatchNumVec, vector<string>& flankStringVec,		
		Index_Info* indexInfo)
	{
		cout << "tmpIns: " << insLen << endl;
		int candiFusChrSeq_length = toFixSeqInRead.size() + 2 + insLen;
		string candiFusChrSeq_left = indexInfo->returnChromStrSubstr(
			toFixChrSeq_chrNameInt_left, toFixChrSeq_chrPos_left, candiFusChrSeq_length);
		string candiFusChrSeq_right = indexInfo->returnChromStrSubstr(
			toFixChrSeq_chrNameInt_right, toFixChrSeq_chrPos_right - candiFusChrSeq_length + 1, 
			candiFusChrSeq_length);
		cout << "candiFusChrSeq_left: " << candiFusChrSeq_left << endl;
		cout << "candiFusChrSeq_right: " << candiFusChrSeq_right << endl;
		int toFixSeqInReadLen = toFixSeqInRead.length();
		cout << "toFixSeqInReadLen: " << toFixSeqInReadLen << endl;
		cout << "toFixSeqInRead: " << toFixSeqInRead << endl;
		// fix canonical fusion junc without inserted bases
		for(int tmpToFixSeqAssignedToLeft = 1; tmpToFixSeqAssignedToLeft < toFixSeqInReadLen;
			tmpToFixSeqAssignedToLeft ++)
		{
			int tmpToFixSeqAssignedToRight = toFixSeqInReadLen - tmpToFixSeqAssignedToLeft - insLen;
			string tmpToFixSeqInRead_withoutInsertedBase 
				= toFixSeqInRead.substr(0, tmpToFixSeqAssignedToLeft) 
					+ toFixSeqInRead.substr(toFixSeqInReadLen - tmpToFixSeqAssignedToRight, tmpToFixSeqAssignedToRight);
			cout << "tmpToFixSeqInRead_withoutInsertedBase: " << tmpToFixSeqInRead_withoutInsertedBase << endl;
			string toFixFusionChrSeq_left_tmp
				= candiFusChrSeq_left.substr(0, tmpToFixSeqAssignedToLeft);
			string toFixFusionChrSeq_right_tmp 
				= candiFusChrSeq_right.substr(candiFusChrSeq_length - tmpToFixSeqAssignedToRight, 
					tmpToFixSeqAssignedToRight);
			cout << "tmpToFixSeqAssignedToLeft: " << tmpToFixSeqAssignedToLeft << endl;
			cout << "tmpToFixSeqAssignedToRight: " << tmpToFixSeqAssignedToRight << endl;
			cout << "toFixFusionChrSeq_left_tmp: " << toFixFusionChrSeq_left_tmp << endl;
			cout << "toFixFusionChrSeq_right_tmp: " << toFixFusionChrSeq_right_tmp << endl;
			string toFixFusionChrSeq_tmp = toFixFusionChrSeq_left_tmp + toFixFusionChrSeq_right_tmp;
			string tmpFlankString_left = candiFusChrSeq_left.substr(tmpToFixSeqAssignedToLeft, 2);
			string tmpFlankString_right = candiFusChrSeq_right.substr(
				candiFusChrSeq_length - tmpToFixSeqAssignedToRight - 2, 2);
			string tmpFlankString = tmpFlankString_left + tmpFlankString_right;
			flankStringVec.push_back(tmpFlankString);
			int tmpMismatchNum = this->compareTwoStringToGetMismatchNum(
				tmpToFixSeqInRead_withoutInsertedBase, toFixFusionChrSeq_tmp);
			mismatchNumVec.push_back(tmpMismatchNum);
			cout << "tmpMismatchNum: " << tmpMismatchNum << endl;
			cout << "tmpFlankString: " << tmpFlankString << endl;
		}
	}

	void fixDoubleAnchor_fusionJunc_vec(int insLenMax, string& toFixSeqInRead, 
		int toFixChrSeq_chrNameInt_left, int toFixChrSeq_chrNameInt_right,
		int toFixChrSeq_chrPos_left, int toFixChrSeq_chrPos_right, 
		vector< vector<int> >& mismatchNumVecVec, 
		vector< vector<string> >& flankStringVecVec,		
		Index_Info* indexInfo)
	{
		vector<int> tmpMismatchNumVec_noIns;
		vector<string> tmpFlankStringVec_noIns;
		cout << "start to do fixDoubleAnchor_fusionJunc_noIns(" << endl;
		this->fixDoubleAnchor_fusionJunc_noIns(toFixSeqInRead, 
			toFixChrSeq_chrNameInt_left, toFixChrSeq_chrNameInt_right,
			toFixChrSeq_chrPos_left, toFixChrSeq_chrPos_right, 
			tmpMismatchNumVec_noIns, tmpFlankStringVec_noIns, indexInfo);
		mismatchNumVecVec.push_back(tmpMismatchNumVec_noIns);
		flankStringVecVec.push_back(tmpFlankStringVec_noIns);

		for(int tmpIns = 1; tmpIns <= insLenMax; tmpIns ++)
		{
			vector<int> tmpMismatchNumVec;
			vector<string> tmpFlankStringVec;
			cout << "start to do fixDoubleAnchor_fusionJunc_withIns(, tmpIns: " << tmpIns << endl;
			this->fixDoubleAnchor_fusionJunc_withIns(tmpIns, toFixSeqInRead, 
				toFixChrSeq_chrNameInt_left, toFixChrSeq_chrNameInt_right,
				toFixChrSeq_chrPos_left, toFixChrSeq_chrPos_right, 
				tmpMismatchNumVec, tmpFlankStringVec, indexInfo);
			mismatchNumVecVec.push_back(tmpMismatchNumVec);
			flankStringVecVec.push_back(tmpFlankStringVec);
		}
	}	

	void getFixFusionJuncParameter_leftHeadFusion(string& toFixSeqInRead,
		int& toFixChrSeq_chrNameInt_left, int& toFixChrSeq_chrPos_left,
		int& toFixChrSeq_chrNameInt_right, int& toFixChrSeq_chrPos_right,
		IncompleteUniquePairedAlignment2detectFusion_Info& incompleteUniquePairedAlignmentInfo, 
		PE_Read_Info& readInfo_unfixedLeftReadHead, 
		PE_Read_Alignment_Info& peAlignInfo_unfixedLeftReadHead,
		//bool toFixFusionIn_leftHead_or_rightTail_bool, 
		bool toFixFusion_Nor_or_Rcm_bool,
		int bufferMax, Index_Info* indexInfo)
	{
		// if(toFixFusionIn_leftHead_or_rightTail_bool)
		// {
			if(toFixFusion_Nor_or_Rcm_bool) // case -- 2,5
			{
				int firstMatchLen_leftOri
					= incompleteUniquePairedAlignmentInfo.returnFirstMatchLength_left();
				int lastMatchLen_leftHeadFusion
					= peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_lastMatchLen();
				int buffer = this->countFixFusionBuffer(bufferMax, firstMatchLen_leftOri, 
					lastMatchLen_leftHeadFusion);
				cout << "firstMatchLen_leftOri: " << firstMatchLen_leftOri << endl;
				cout << "lastMatchLen_leftHeadFusion: " << lastMatchLen_leftHeadFusion << endl;
				cout << "buffer: " << buffer << endl;
				int tailSoftClipLen_leftHeadFusion
					= peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_tailSoftClipLen();				
				int oriUnfixedLeftHeadLen = incompleteUniquePairedAlignmentInfo.returnUnfixedHeadLength_left();
				int toFixSeqInRead_length = tailSoftClipLen_leftHeadFusion + buffer * 2;
				int toFixSeqInRead_lastLocInRead = oriUnfixedLeftHeadLen + buffer;
				int toFixSeqInRead_startLocInRead = toFixSeqInRead_lastLocInRead - toFixSeqInRead_length + 1;
				toFixSeqInRead = incompleteUniquePairedAlignmentInfo.returnReadSeqSubStr_left(
					toFixSeqInRead_startLocInRead, toFixSeqInRead_lastLocInRead);
				cout << "toFixSeqInRead: " << toFixSeqInRead << endl;
				toFixChrSeq_chrNameInt_left = peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_chrNameInt(indexInfo);
				cout << "toFixChrSeq_chrNameInt_left: " << toFixChrSeq_chrNameInt_left << endl;
				toFixChrSeq_chrNameInt_right = incompleteUniquePairedAlignmentInfo.returnChrNameInt();
				cout << "toFixChrSeq_chrNameInt_right" << toFixChrSeq_chrNameInt_right << endl;
				toFixChrSeq_chrPos_left = peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_endMatchedPosInChr() - buffer + 1;
				cout << "toFixChrSeq_chrPos_left: " << toFixChrSeq_chrPos_left << endl;
				toFixChrSeq_chrPos_right = incompleteUniquePairedAlignmentInfo.returnStartPos_1() + buffer - 1;
				cout << "toFixChrSeq_chrPos_right: " << toFixChrSeq_chrPos_right << endl;
			}
			else // case --  10,11
			{
				int firstMatchLen_leftOri
					= incompleteUniquePairedAlignmentInfo.returnFirstMatchLength_left();
				int firstMatchLen_leftHeadFusion
					= peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_firstMatchLen();					
				int buffer = this->countFixFusionBuffer(bufferMax, firstMatchLen_leftOri, 
					firstMatchLen_leftHeadFusion);

				int headSoftClipLen_leftHeadFusion
					= peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_headSoftClipLen();				
				int oriUnfixedLeftHeadLen = incompleteUniquePairedAlignmentInfo.returnUnfixedHeadLength_left();
				int toFixSeqInRead_length = headSoftClipLen_leftHeadFusion + buffer * 2;
				int toFixSeqInRead_lastLocInRead = oriUnfixedLeftHeadLen + buffer;
				int toFixSeqInRead_startLocInRead = toFixSeqInRead_lastLocInRead - toFixSeqInRead_length + 1;
				toFixSeqInRead = incompleteUniquePairedAlignmentInfo.returnReadSeqSubStr_left(
					toFixSeqInRead_startLocInRead, toFixSeqInRead_lastLocInRead);
				toFixChrSeq_chrNameInt_left = peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_chrNameInt(indexInfo);
				toFixChrSeq_chrNameInt_right = incompleteUniquePairedAlignmentInfo.returnChrNameInt();
				toFixChrSeq_chrPos_left = peAlignInfo_unfixedLeftReadHead.return_onlySeAlignInfo_startMatchedPosInChr() + buffer - 1;
				toFixChrSeq_chrPos_right = incompleteUniquePairedAlignmentInfo.returnStartPos_1() + buffer - 1;
			}
		// }
		// else
		// {
		// 	if(toFixFusion_Nor_or_Rcm_bool)
		// 	{}
		// 	else
		// 	{}
		// }
	}

	bool fixAndChooseBestFusionSite(int allowedMismatchNumMax, 
		vector< vector<int> >& mismatchNumVecVec, 
		vector< vector<string> >& flankStringVecVec,
		int& bestFusionSite_assignedSeqLenToLeft, int& bestFusionSite_insLen)
	{

		return true;
	}

	bool refineFusionBreakPoint_leftHead(int insLenMax, int min_fusion_distance, 
		IncompleteUniquePairedAlignment2detectFusion_Info& incompleteUniquePairedAlignmentInfo, 
		PE_Read_Info& readInfo_unfixedLeftReadHead, 
		PE_Read_Alignment_Info& peAlignInfo_unfixedLeftReadHead, Index_Info* indexInfo)
	{
		int bufferMax = 2;
		cout << "refineFusionBreakPoint_leftHead( starts ......" << endl;
		int tmpSeAlignVec_final_Nor_vecSize 
			= peAlignInfo_unfixedLeftReadHead.return_seAlignVec_final_Nor_size(); //seAlignVec_final_Nor.size();
		int tmpSeAlignVec_final_Rcm_vecSize 
			= peAlignInfo_unfixedLeftReadHead.return_seAlignVec_final_Rcm_size(); // seAlignVec_final_Rcm.size();
		cout << "tmpSeAlignVec_final_Nor_vecSize: " << tmpSeAlignVec_final_Nor_vecSize << endl;
		cout << "tmpSeAlignVec_final_Rcm_vecSize: " << tmpSeAlignVec_final_Rcm_vecSize << endl;
		if(tmpSeAlignVec_final_Nor_vecSize + tmpSeAlignVec_final_Rcm_vecSize != 1)
			return false;
		vector< vector<int> > mismatchNumVecVec;
		vector< vector<string> > flankStringVecVec;
		int bestFusionSite_assignedSeqLenToLeft;
		int bestFusionSite_insLen;
		bool Nor_or_Rcm_bool = (tmpSeAlignVec_final_Nor_vecSize == 1);

		if(Nor_or_Rcm_bool) // fusionCase -- 2,5
		{
			string toFixSeqInRead;
			int toFixChrSeq_chrNameInt_left, toFixChrSeq_chrPos_left;
			int toFixChrSeq_chrNameInt_right, toFixChrSeq_chrPos_right;
			cout << "start to getFixFusionJuncParameter_leftHeadFusion(" << endl;
			this->getFixFusionJuncParameter_leftHeadFusion(toFixSeqInRead,
				toFixChrSeq_chrNameInt_left, toFixChrSeq_chrPos_left,
				toFixChrSeq_chrNameInt_right, toFixChrSeq_chrPos_right,
				incompleteUniquePairedAlignmentInfo, 
				readInfo_unfixedLeftReadHead, peAlignInfo_unfixedLeftReadHead,
				true, bufferMax, indexInfo);
			cout << "start to do fixDoubleAnchor_fusionJunc_vec(" << endl;
			this->fixDoubleAnchor_fusionJunc_vec(insLenMax, toFixSeqInRead, 
				toFixChrSeq_chrNameInt_left, toFixChrSeq_chrNameInt_right,
				toFixChrSeq_chrPos_left, toFixChrSeq_chrPos_right, 
				mismatchNumVecVec, flankStringVecVec, indexInfo);
			int allowedMismatchNumMax 
				= toFixSeqInRead.length() / AllowedMismatchNumMaxPerBaseNum;
			cout << "allowedMismatchNumMax: " << allowedMismatchNumMax << endl;
			cout << "start to do fixAndChooseBestFusionSite( " << endl;
			bool fixAndChooseBestFusionSite_bool = this->fixAndChooseBestFusionSite(
				allowedMismatchNumMax, mismatchNumVecVec, flankStringVecVec,
				bestFusionSite_assignedSeqLenToLeft, bestFusionSite_insLen);
			return fixAndChooseBestFusionSite_bool;
		}
		else // fusionCase -- 10,11
		{

		}
	}

	bool refineFusionBreakPoint_rightTail(int min_fusion_distance, 
		IncompleteUniquePairedAlignment2detectFusion_Info& incompleteUniquePairedAlignmentInfo,
		PE_Read_Info& readInfo_unfixedRightReadTail, 
		PE_Read_Alignment_Info& peAlignInfo_unfixedRightReadTail, Index_Info* indexInfo)
	{
		cout << "refineFusionBreakPoint_rightTail( starts ......" << endl;
		int tmpSeAlignVec_final_Nor_vecSize 
			= peAlignInfo_unfixedRightReadTail.return_seAlignVec_final_Nor_size(); //seAlignVec_final_Nor.size();
		int tmpSeAlignVec_final_Rcm_vecSize 
			= peAlignInfo_unfixedRightReadTail.return_seAlignVec_final_Rcm_size(); //seAlignVec_final_Rcm.size();
		if(tmpSeAlignVec_final_Nor_vecSize + tmpSeAlignVec_final_Rcm_vecSize != 1)
			return false;
		vector< vector<int> > mismatchNumVecVec;
		vector< vector<string> > flankStringVecVec;
		bool Nor_or_Rcm_bool = (tmpSeAlignVec_final_Nor_vecSize == 1);
		if(Nor_or_Rcm_bool) // fusionCase -- 1,4
		{

		}
		else // fusionCase -- 7,8
		{

		}		
	}

	void refineFusionBreakPoint(int insLenMax, int min_fusion_distance, 
		IncompleteUniquePairedAlignment2detectFusion_Info& incompleteUniquePairedAlignmentInfo, 
		bool leftReadHeadUnfixed_bool, PE_Read_Info& readInfo_unfixedLeftReadHead, 
		PE_Read_Alignment_Info& peAlignInfo_unfixedLeftReadHead,
		bool rightReadTailUnfixed_bool, PE_Read_Info& readInfo_unfixedRightReadTail, 
		PE_Read_Alignment_Info& peAlignInfo_unfixedRightReadTail, Index_Info* indexInfo)
	{
		cout << "refineFusionBreakPoint( starts ......" << endl;
		if(leftReadHeadUnfixed_bool)
			leftReadHeadUnfixed_fixed_bool = this->refineFusionBreakPoint_leftHead(
				insLenMax, min_fusion_distance, 
				incompleteUniquePairedAlignmentInfo, readInfo_unfixedLeftReadHead,
				peAlignInfo_unfixedLeftReadHead, indexInfo);
		else
		{
			// rightReadHeadUnfixed_fixed_bool = this->refineFusionBreakPoint_rightTail(
			// 	insLenMax, min_fusion_distance,
			// 	incompleteUniquePairedAlignmentInfo, readInfo_unfixedRightReadTail,
			// 	peAlignInfo_unfixedRightReadTail, indexInfo);
		}
	}

	string returnFusionSamStr()
	{}

	string returnFusionBreakPointStr()
	{}
	
	bool return_leftReadHeadUnfixed_bool()
	{
		return leftReadHeadUnfixed_bool;
	}
	bool return_rightReadHeadUnfixed_bool()
	{
		return rightReadTailUnfixed_bool;
	}	
	bool return_leftReadHeadUnfixed_fixed_bool()
	{
		return leftReadHeadUnfixed_fixed_bool;
	}
	bool return_rightReadTailUnfixed_fixed_bool()
	{
		return rightReadTailUnfixed_fixed_bool;
	}
	bool return_leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool()
	{
		return leftReadHeadUnfixed_fixed_Nor_or_Rcm_bool;
	}
	bool return_rightReadTailUnfixed_fixed_Nor_or_Rcm_bool()
	{
		return rightReadTailUnfixed_fixed_Nor_or_Rcm_bool;
	}
};
#endif