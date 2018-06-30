// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FUSIONBREAKPOINT_INFO_H
#define FUSIONBREAKPOINT_INFO_H

using namespace std;

class FusionBreakPoint_Info
{
private:
	int chrNameInt_1;
	int chrNameInt_2;

	int breakPoint_1;
	int breakPoint_2;

	int supportNum; // support num for fusion-spanning reads
	int supportNum_encompassing; 

	string strand_1;
	string strand_2;

	string flankString;

	int anchorLength_1;
	int anchorLength_2;

	int XM_min;
	int XM_max;

	vector<int> alterFusionIndexVec_1stEndShared;
	vector<int> alterFusionIndexVec_2ndEndShared;

	string fusionCaseStr;
	int spanningReadsMapRangeMax_1;
	int spanningReadsMapRangeMax_2;

public:
	FusionBreakPoint_Info()
	{
		supportNum = 0;
		supportNum_encompassing = 0;
		anchorLength_1 = 0;
		anchorLength_2 = 0;
		flankString = "NULL";
		XM_min = 0;
		XM_max = 0;
	    fusionCaseStr = "NULL";
	    spanningReadsMapRangeMax_1 = 0;
	    spanningReadsMapRangeMax_2 = 0;
	}

	void copyAlterFusionIndexVec_exceptThisFusionIndex(
		int thisFusionIndex,
		vector<int>& tmpAlterFusionIndexVec_1stEndShared,
		vector<int>& tmpAlterFusionIndexVec_2ndEndShared)
	{
		//cout << "tmpAlterFusionIndexVec_1stEndShared.size(): " << tmpAlterFusionIndexVec_1stEndShared.size() << endl;
		for(int tmp = 0; tmp < tmpAlterFusionIndexVec_1stEndShared.size(); tmp++)
		{
			int tmpIndex = tmpAlterFusionIndexVec_1stEndShared[tmp];
			//cout << "tmpIndex: " << tmpIndex << endl;
			if(tmpIndex != thisFusionIndex)
				alterFusionIndexVec_1stEndShared.push_back(tmpIndex);
		}
		//cout << "tmpAlterFusionIndexVec_2ndEndShared.size(): " << tmpAlterFusionIndexVec_2ndEndShared.size() << endl;
		for(int tmp = 0; tmp < tmpAlterFusionIndexVec_2ndEndShared.size(); tmp++)
		{
			int tmpIndex = tmpAlterFusionIndexVec_2ndEndShared[tmp];
			//cout << "tmpIndex: " << tmpIndex << endl;
			if(tmpIndex != thisFusionIndex)
				alterFusionIndexVec_2ndEndShared.push_back(tmpIndex);
		}		
	}

	int returnAlterFusionNum_1stEndShared()
	{
		return alterFusionIndexVec_1stEndShared.size();
	}

	int returnAlterFusionNum_2ndEndShared()
	{
		return alterFusionIndexVec_2ndEndShared.size();
	}

	int returnAlterFusionIndex_1stEndShared(int tmp)
	{
		return alterFusionIndexVec_1stEndShared[tmp];
	}

	int returnAlterFusionIndex_2ndEndShared(int tmp)
	{
		return alterFusionIndexVec_2ndEndShared[tmp];
	}

	string returnFusionInfoStr_forAlterFusionDetection(Index_Info* indexInfo)
	{
		string tmpChrNameStr_1 = indexInfo->returnChrNameStr(chrNameInt_1);
		string tmpChrNameStr_2 = indexInfo->returnChrNameStr(chrNameInt_2);
		string tmpBreakPointStr_1 = int_to_str(breakPoint_1);
		string tmpBreakPointStr_2 = int_to_str(breakPoint_2);
		string tmpFusionInfoStr = tmpChrNameStr_1 + ":" + tmpBreakPointStr_1 + ":"
			+ tmpChrNameStr_2 + ":" + tmpBreakPointStr_2 + ":" + strand_1 + strand_2;
		return tmpFusionInfoStr;
	}

	void updateXM(int tmpXMmin, int tmpXMmax)
	{
		if(tmpXMmin < XM_min)
			XM_min = tmpXMmin;
		if(tmpXMmax > XM_max)
			XM_max = tmpXMmax;
	}

	void updateAnchorSize(int tmpAnchorLength_1, int tmpAnchorLength_2)
	{
		if(tmpAnchorLength_1 > anchorLength_1)
			anchorLength_1 = tmpAnchorLength_1;
		if(tmpAnchorLength_2 > anchorLength_2)
			anchorLength_2 = tmpAnchorLength_2;
	}

	void clearSupportNum()
	{
		supportNum = 0;
	}

	void clearSupportNum_encompassing()
	{
		supportNum_encompassing = 0;
	}

	string returnFusionBreakPointStr(Index_Info* indexInfo)
	{
		string tmpStr;
		string chrNameStr_1 = indexInfo->returnChrNameStr(chrNameInt_1);
		string chrNameStr_2 = indexInfo->returnChrNameStr(chrNameInt_2);
		tmpStr = chrNameStr_1 + "\t" + chrNameStr_2 + "\t"
			+ int_to_str(breakPoint_1) + "\t" + int_to_str(breakPoint_2) + "\t"
			+ strand_1 + "\t" + strand_2 + "\t" + flankString + "\t"
			+ int_to_str(anchorLength_1) + "\t" + int_to_str(anchorLength_2) + "\t"
			+ int_to_str(supportNum) + "\t" + int_to_str(supportNum_encompassing) + "\t"
			+ int_to_str(XM_min) + "\t" + int_to_str(XM_max);
		return tmpStr;
	}

	string returnFusionBreakPointStr_withAlterFusionJunc(Index_Info* indexInfo)
	{
		string tmpStr;
		string chrNameStr_1 = indexInfo->returnChrNameStr(chrNameInt_1);
		string chrNameStr_2 = indexInfo->returnChrNameStr(chrNameInt_2);
		tmpStr = chrNameStr_1 + "\t" + chrNameStr_2 + "\t"
			+ int_to_str(breakPoint_1) + "\t" + int_to_str(breakPoint_2) + "\t"
			+ strand_1 + "\t" + strand_2 + "\t" + flankString + "\t"
			+ int_to_str(anchorLength_1) + "\t" + int_to_str(anchorLength_2) + "\t"
			+ int_to_str(supportNum) + "\t" + int_to_str(supportNum_encompassing) + "\t"
			+ int_to_str(XM_min) + "\t" + int_to_str(XM_max);
		return tmpStr;
	}	

	void addSupportNum(int tmpSupNum)
	{
		supportNum += tmpSupNum;
	}

	void addSupportNum_encompassing(int tmpSupNum)
	{
		supportNum_encompassing += tmpSupNum;
	}

	int returnChrNameInt_1()
	{
		return chrNameInt_1;
	}

	int returnChrNameInt_2()
	{
		return chrNameInt_2;
	}

	int returnBreakPoint_1()
	{
		return breakPoint_1;
	}

	int returnBreakPoint_2()
	{
		return breakPoint_2;
	}

	int returnSupportNum()
	{
		return supportNum;
	}

	int returnSupportNum_encompassing()
	{
		return supportNum_encompassing;
	}

	string returnStrand_1()
	{
		return strand_1;
	}

	string returnStrand_2()
	{
		return strand_2;
	}

	int returnAnchorLength_1()
	{
		return anchorLength_1;
	}

	int returnAnchorLength_2()
	{
		return anchorLength_2;
	}	

	int returnXMmin()
	{
		return XM_min;
	}

	int returnXMmax()
	{
		return XM_max;
	}

	void initiateWithChrNameBreakPoint(
		string& tmpChrNameStr_1, string& tmpChrNameStr_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupportNum, int tmpSupportNum_encompassing,
		string& tmpFlankString, int tmpAnchorLength_1, int tmpAnchorLength_2,
		Index_Info* indexInfo)
	{
		chrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
		chrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
		breakPoint_1 = tmpBreakPoint_1;
		breakPoint_2 = tmpBreakPoint_2;
	    supportNum = tmpSupportNum;
	    supportNum_encompassing = tmpSupportNum_encompassing;
	    strand_1 = tmpStrand_1;
	    strand_2 = tmpStrand_2;
	    flankString = tmpFlankString;
	    anchorLength_1 = tmpAnchorLength_1;
	    anchorLength_2 = tmpAnchorLength_2;
	}

	void initiateWithChrNameBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupportNum, int tmpSupportNum_encompassing,
		string& tmpFlankString, int tmpAnchorLength_1, int tmpAnchorLength_2)
	{
		chrNameInt_1 = tmpChrNameInt_1;
		chrNameInt_2 = tmpChrNameInt_2;
		breakPoint_1 = tmpBreakPoint_1;
		breakPoint_2 = tmpBreakPoint_2;
	    supportNum = tmpSupportNum;
	    supportNum_encompassing = tmpSupportNum_encompassing;
	    strand_1 = tmpStrand_1;
	    strand_2 = tmpStrand_2;
	    flankString = tmpFlankString;
	    anchorLength_1 = tmpAnchorLength_1;
	    anchorLength_2 = tmpAnchorLength_2;
	}

	void initiateWithChrNameBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupportNum, int tmpSupportNum_encompassing,
		string& tmpFlankString, int tmpAnchorLength_1, int tmpAnchorLength_2,
		int tmpXMmin, int tmpXMmax)
	{
		chrNameInt_1 = tmpChrNameInt_1;
		chrNameInt_2 = tmpChrNameInt_2;
		breakPoint_1 = tmpBreakPoint_1;
		breakPoint_2 = tmpBreakPoint_2;
	    supportNum = tmpSupportNum;
	    supportNum_encompassing = tmpSupportNum_encompassing;
	    strand_1 = tmpStrand_1;
	    strand_2 = tmpStrand_2;
	    flankString = tmpFlankString;
	    anchorLength_1 = tmpAnchorLength_1;
	    anchorLength_2 = tmpAnchorLength_2;
	    XM_min = tmpXMmin;
	    XM_max = tmpXMmax;
	}

	void initiateWithChrNameBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupportNum, int tmpSupportNum_encompassing,
		string& tmpFlankString, int tmpAnchorLength_1, int tmpAnchorLength_2,
		int tmpXMmin, int tmpXMmax, string& tmpFusionCaseStr,
		int tmpMapRangeMax_1, int tmpMapRangeMax_2)
	{
		chrNameInt_1 = tmpChrNameInt_1;
		chrNameInt_2 = tmpChrNameInt_2;
		breakPoint_1 = tmpBreakPoint_1;
		breakPoint_2 = tmpBreakPoint_2;
	    supportNum = tmpSupportNum;
	    supportNum_encompassing = tmpSupportNum_encompassing;
	    strand_1 = tmpStrand_1;
	    strand_2 = tmpStrand_2;
	    flankString = tmpFlankString;
	    anchorLength_1 = tmpAnchorLength_1;
	    anchorLength_2 = tmpAnchorLength_2;
	    XM_min = tmpXMmin;
	    XM_max = tmpXMmax;
	    fusionCaseStr = tmpFusionCaseStr;
	    spanningReadsMapRangeMax_1 = tmpMapRangeMax_1;
	    spanningReadsMapRangeMax_2 = tmpMapRangeMax_2;
	}		

	void initiateWithChrNameBreakPoint(
		string& tmpChrNameStr_1, string& tmpChrNameStr_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupportNum, int tmpSupportNum_encompassing,
		Index_Info* indexInfo)
	{
		chrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
		chrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
		breakPoint_1 = tmpBreakPoint_1;
		breakPoint_2 = tmpBreakPoint_2;
	    supportNum = tmpSupportNum;
	    supportNum_encompassing = tmpSupportNum_encompassing;
	    strand_1 = tmpStrand_1;
	    strand_2 = tmpStrand_2;
	}

	void initiateWithChrNameBreakPoint(
		int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPoint_1, int tmpBreakPoint_2, 
		string& tmpStrand_1, string& tmpStrand_2, 
		int tmpSupportNum, int tmpSupportNum_encompassing)
	{
		chrNameInt_1 = tmpChrNameInt_1;
		chrNameInt_2 = tmpChrNameInt_2;
		breakPoint_1 = tmpBreakPoint_1;
		breakPoint_2 = tmpBreakPoint_2;
	    supportNum = tmpSupportNum;
	    supportNum_encompassing = tmpSupportNum_encompassing;
	    strand_1 = tmpStrand_1;
	    strand_2 = tmpStrand_2;
	}
};
#endif