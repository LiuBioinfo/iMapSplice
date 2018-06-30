// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// 0 -- exactly the same
// 1 -- totally different
// 2 -- overlapped, but contradictory in some bases
// 3 -- without contradictory, not totally the same, covered bases are a parent set with another
// 4 -- without contradictory, not totally the same, covered bases are a child set with another

#ifndef ALIGNMENT_INFO_H
#define ALIGNMENT_INFO_H

class Alignment_Info
{
private:
	string readName;
	int chrNameInt;
	int chrPos_start;
	int chrPos_end;
	vector<Jump_Code> cigarJumpCodeVec;

public:
	Alignment_Info()
	{}

	void generateAlignmentInfo_gt(string& tmpGtStr, Index_Info* indexInfo)
	{

	}

	void generateAlignmentInfo_MPS3(string& tmpMPS3str, Index_Info* indexInfo)
	{
		this->generateAlignmentInfo_gt(tmpMPS3str, indexInfo);
	}

	bool cigarJumpCodeVecExactlyTheSameWithAnotherAlignmentInfo(
		Alignment_Info& otherAlignmentInfo)
	{
		int thisJumpCodeVecSize = cigarJumpCodeVec.size();
		int otherJumpCodeVecSize = otherAlignmentInfo.returnCigarJumpCodeVecSize();
		if(thisJumpCodeVecSize != otherJumpCodeVecSize)
			return false;
		for(int tmp = 0; tmp < cigarJumpCodeVec.size(); tmp)
		{
			string tmpThisJumpCodeType = cigarJumpCodeVec[tmp].type;
			int tmpThisJumpCodeLength = cigarJumpCodeVec[tmp].len;
			string tmpOtherJumpCodeType = otherAlignmentInfo.returnJumpCodeType();
			int tmpOtherJumpCodeLength = otherAlignmentInfo.returnJumpCodeLength();
			if((tmpThisJumpCodeType != tmpOtherJumpCodeType)||(tmpThisJumpCodeLength != tmpOtherJumpCodeLength))
				return false;
		}
		return true;
	}

	/*
	int compareTwoExonicMapRegionVec(
		vector< pair<int,int> >& exonicMapRegionVec_1,
		vector< pair<int,int> >& exonicMapRegionVec_2)
	{
		// 0 -- the same
		// 1 -- contradictory 
		// 2 -- not contradictory, exonicMapRegionVec_1 is a parent set of exonicMapRegionVec_2
		// 3 -- not contradictory, exonicMapRegionvec_1 is a child set of exonicMapRegionVec_2
		int exonNum_1 = exonicMapRegionVec_1.size();
		int exonNum_2 = exonicMapRegionVec_2.size();
		bool baseInExonicMapRegionVec_1_not_2 = false;
		bool baseInExonicMapRegionVec_2_not_1 = false;
		if(exonNum_1 == exonNum_2)
		{
			for(int tmp = 0; tmp < exonNum_1; tmp++)
			{
				int tmpExon_startPos_1 = exonicMapRegionVec_1[tmp].first;
				int tmpExon_endPos_1 = exonicMapRegionVec_1[tmp].second;
				int tmpExon_startPos_2 = exonicMapRegionVec_2[tmp].first;
				int tmpExon_endPos_2 = exonicMapRegionVec_2[tmp].second;
				if(( != )&&())
			}
		}
		else if(exonNum_1 > exonNum_2)
		{

		}
		else // exonNum_1 < exonNum_2
		{

		}
	}*/

	int returnNextSJindexInJumpCodeVec_withStartIndex(int startIndex)
	{
		for(int tmp = startIndex; tmp < cigarJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
				return tmp;
		}
		return -1;
	}

	int returnCigarJumpCodeVecSize()
	{
		return cigarJumpCodeVec.size();
	}

	int returnExonNumFromCigarJumpCodeVec()
	{
		int tmpExonNum = 0;
		for(int tmp = 0; tmp < cigarJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
				tmpExonNum ++;
		}
		return tmpExonNum;
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

	void getExonicMapRegionVec(vector< pair<int,int> >& exonicMapRegionVec)
	{
		int tmpExonNum = this->returnExonNumFromCigarJumpCodeVec();
		if(tmpExonNum == 0)
			exonicMapRegionVec.push_back(pair<int,int>(chrPos_start, chrPos_end));
		else if(tmpExonNum == 1)
		{
			int tmpSJindexInJumpCodeVec = this->returnNextSJindexInJumpCodeVec_withStartIndex(0);
			int the1stExonicRegion_endPos = this->getEndPosOfSpecificJumpCode(tmpSJindexInJumpCodeVec - 1);
			int tmp2ndExonicRegion_startPos = this->getEndPosOfSpecificJumpCode(tmpSJindexInJumpCodeVec) + 1;
			exonicMapRegionVec.push_back(pair<int,int>(chrPos_start, the1stExonicRegion_endPos));
			exonicMapRegionVec.push_back(pair<int,int>(tmp2ndExonicRegion_startPos, chrPos_end));
		}
		else // exonNum >= 2
		{
			// push back the 1st exonic region
			int tmpSJindexInJumpCodeVec_1st = this->returnNextSJindexInJumpCodeVec_withStartIndex(0);
			int the1stExonicRegion_endPos = this->getEndPosOfSpecificJumpCode(tmpSJindexInJumpCodeVec_1st - 1);
			exonicMapRegionVec.push_back(pair<int,int>(chrPos_start, the1stExonicRegion_endPos));
			// push back mid exonic regions
			int tmpLastSJindexInJumpCode = tmpSJindexInJumpCodeVec_1st;
			int tmpNextSJindexInJumpCode;
			for(int tmp = 0; tmp < tmpExonNum - 1 - 1; tmp++)
			{
				int tmpExonStartIndexInJumpCodeVec = tmpLastSJindexInJumpCode + 1;
				tmpNextSJindexInJumpCode = this->returnNextSJindexInJumpCodeVec_withStartIndex(
					tmpLastSJindexInJumpCode + 1);
				int tmpExonEndIndexInJumpCodeVec = tmpNextSJindexInJumpCode - 1;
				tmpLastSJindexInJumpCode = tmpNextSJindexInJumpCode;
				int tmpExonicRegion_startPos = this->getEndPosOfSpecificJumpCode(tmpExonStartIndexInJumpCodeVec-1) + 1;
				int tmpExonicRegion_endPos = this->getEndPosOfSpecificJumpCode(tmpExonEndIndexInJumpCodeVec);
				exonicMapRegionVec.push_back(pair<int,int>(tmpExonicRegion_startPos, tmpExonicRegion_endPos));
			}
			// push back the last exonic region
			int tmpLastExonicRegion_startPos = this->getEndPosOfSpecificJumpCode(tmpLastSJindexInJumpCode) + 1;
			exonicMapRegionVec.push_back(pair<int,int>(tmpLastExonicRegion_startPos, chrPos_end));
		}
	}

	int compareWithAnotherAlignmentInfo(
		Alignment_Info& anotherAlignmentInfo)
	{
		// 0 -- exactly the same
		// 1 -- totally different
		// 2 -- overlapped, but contradictory in some bases
		// 3 -- without contradictory, not totally the same, covered bases are a parent set with another
		// 4 -- without contradictory, not totally the same, covered bases are a child set with another
		int chrNameInt_other = anotherAlignmentInfo.returnChrNameInt();
		int chrPosStart_other = anotherAlignmentInfo.returnChrPosStart();
		int chrPosEnd_other = anotherAlignmentInfo.returnChrPosEnd();
		if(chrNameInt != chrNameInt_other)
			return 1;
		if((chrPos_start > chrPosEnd_other)||(chrPos_end < chrPosStart_other))
			return 1;
		if((chrPos_start == chrPosStart_other)&&(chrPos_end == chrPosEnd_other)
			&&(this->cigarJumpCodeVecExactlyTheSameWithAnotherAlignmentInfo(anotherAlignmentInfo)))
			return 0;
		if((chrPos_start <= chrPosStart_other)||(chrPosEnd_other <= chrPos_end)) // 2 or 3
		{
			vector< pair<int,int> > exonicMapRegionVec_this;
			vector< pair<int,int> > exonicMapRegionVec_other;
			this->getExonicMapRegionVec(exonicMapRegionVec_this);
			this->getExonicMapRegionVec(exonicMapRegionVec_other);
			// check if each exonicMapRegion in exonicMapRegionVec_other is a subset of exonicMapRegionVec_this;
			for(int tmpExonIndex_other = 0; tmpExonIndex_other < exonicMapRegionVec_other.size(); tmpExonIndex_other ++)
			{
				int tmpExon_startPos_other = exonicMapRegionVec_other[tmpExonIndex_other].first;
				int tmpExon_endPos_other = exonicMapRegionVec_other[tmpExonIndex_other].second;
				bool tmpOtherExonWithinTmpThisExonBool = false;
				for(int tmpExonIndex_this = 0; tmpExonIndex_this < exonicMapRegionVec_this.size(); tmpExonIndex_this ++)
				{
					int tmpExon_startPos_this = exonicMapRegionVec_this[tmpExonIndex_this].first;
					int tmpExon_endPos_this = exonicMapRegionVec_this[tmpExonIndex_this].second;					
					if((tmpExon_startPos_other >= tmpExon_startPos_this)&&(tmpExon_endPos_other <= tmpExon_endPos_this))
					{
						tmpOtherExonWithinTmpThisExonBool = true;
						break;
					}
				}
				if(tmpOtherExonWithinTmpThisExonBool)
				{}	
				else
					return 2;
			}
			return 3;
		}
		else if((chrPosStart_other <= chrPos_start)||(chrPos_end <= chrPosEnd_other)) // 2 or 4
		{
			vector< pair<int,int> > exonicMapRegionVec_this;
			vector< pair<int,int> > exonicMapRegionVec_other;
			this->getExonicMapRegionVec(exonicMapRegionVec_this);
			this->getExonicMapRegionVec(exonicMapRegionVec_other);
			// check if each exonicMapRegion in exonicMapRegionVec_this is a subset of exonicMapRegionVec_other;
			for(int tmpExonIndex_this = 0; tmpExonIndex_this < exonicMapRegionVec_this.size(); tmpExonIndex_this ++)
			{
				int tmpExon_startPos_this = exonicMapRegionVec_this[tmpExonIndex_this].first;
				int tmpExon_endPos_this = exonicMapRegionVec_this[tmpExonIndex_this].second;
				bool tmpThisExonWithinTmpOtherExonBool =false;
				for(int tmpExonIndex_other = 0; tmpExonIndex_other < exonicMapRegionVec_other.size(); tmpExonIndex_other ++)
				{
					int tmpExon_startPos_other = exonicMapRegionVec_other[tmpExonIndex_other].first;
					int tmpExon_endPos_other = exonicMapRegionVec_other[tmpExonIndex_other].second;					
					if((tmpExon_startPos_this >= tmpExon_startPos_other)&&(tmpExon_endPos_this <= tmpExon_endPos_other))
					{
						tmpThisExonWithinTmpOtherExonBool = true;
						break;
					}
				}
				if(tmpThisExonWithinTmpOtherExonBool)
				{}	
				else
					return 2;
			}
			return 4;
		}
		else
		{
			return 2;
		}
	}

	int returnChrPosEnd()
	{
		return chrPos_end;
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnChrPosStart()
	{
		return chrPos_start;
	}

	string returnReadName()
	{
		return readName;
	}

	int returnJumpCodeLength(int index)
	{
		return cigarJumpCodeVec[index].len;
	}

	string returnJumpCodeType(int index)
	{
		return cigarJumpCodeVec[index].type;
	}
};

int getReadNum_gt(string& gtStr_1, string& gtStr_2)
{
	int tabLoc_dot_1 = gtStr_1.find(".");
	int tabLoc_end1sign = gtStr_1.find("/1");

	int tabLoc_dot_2 = gtStr_2.find(".");
	int tabLoc_end2sign = gtStr_2.find("/2");

	if((tabLoc_dot_1 == string::npos)
		||(tabLoc_end1sign == string::npos)
		||(tabLoc_dot_2 == string::npos)
		||(tabLoc_end2sign == string::npos))
		return -1;

	string tmpReadNumStr_1 = gtStr_1.substr(
		tabLoc_dot_1 + 1, tabLoc_end1sign - tabLoc_dot_1 - 1);
	string tmpReadNumStr_2 = gtStr_2.substr(
		tabLoc_dot_2 + 1, tabLoc_end2sign - tabLoc_dot_2 - 1);

	if(tmpReadNumStr_1 == tmpReadNumStr_2)
		return atoi(tmpReadNumStr_1.c_str());
	else
		return -1;
}

int getReadNum_MPS3(string& samStr_1, string& samStr_2)
{
	int tabLoc_dot_1 = samStr_1.find(".");
	int tabLoc_dot_2 = samStr_2.find(".");
	if((tabLoc_dot_1 == string::npos)
		||(tabLoc_dot_2 == string::npos))
		return -1;
	string tmpReadNumStr_1 = samStr_1.substr(tabLoc_dot_1 + 1);
	string tmpReadNumStr_2 = samStr_2.substr(tabLoc_dot_2 + 1);
	
	if(tmpReadNumStr_1 == tmpReadNumStr_2)
		return atoi(tmpReadNumStr_1.c_str());
	else
		return -1;
}

void generateGroundTruthSamInfoVec(
	vector< Alignment_Info* >& groundTruthSamInfoVec_end1,
	vector< Alignment_Info* >& groundTruthSamInfoVec_end2,
	vector< bool >& groundTruthExistBoolVec_end1,
	vector< bool >& groundTruthExistBoolVec_end2,
	string& inputGroundTruthSamFile, Index_Info* indexInfo)
{
	ifstream gt_ifs(inputGroundTruthSamFile.c_str());
	while(!gt_ifs.eof())
	{
		string gtStr_1;
		getline(gt_ifs, gtStr_1);
		if(gtStr_1 == "")
			break;
		string gtStr_2;
		getline(gt_ifs, gtStr_2);
	
		int tmpReadNum = getReadNum_gt(gtStr_1, gtStr_2);
		if(tmpReadNum <= 0)
			continue;
		groundTruthSamInfoVec_end1[tmpReadNum - 1]->generateAlignmentInfo_gt(gtStr_1, indexInfo);
		groundTruthSamInfoVec_end2[tmpReadNum - 1]->generateAlignmentInfo_gt(gtStr_2, indexInfo);
		groundTruthExistBoolVec_end1[tmpReadNum - 1] = true;
		groundTruthExistBoolVec_end2[tmpReadNum - 1] = true;
	}
	gt_ifs.close();
}

void generateAlignerSamInfoVec(
	vector< Alignment_Info* >& alignerSamInfoVec_end1,
	vector< Alignment_Info* >& alignerSamInfoVec_end2,
	vector< bool >& alignerSamExistBoolVec_end1,
	vector< bool >& alignerSamExistBoolVec_end2,
	string& inputAlignerSamFile, Index_Info* indexInfo)
{
	ifstream alignerSam_ifs(inputAlignerSamFile.c_str());
	while(!alignerSam_ifs.eof())
	{
		string tmpAlignerSamStr_1;
		getline(alignerSam_ifs, tmpAlignerSamStr_1);
		if(tmpAlignerSamStr_1.substr(0,1) == "@")
			continue;
		if(tmpAlignerSamStr_1 == "")
			break;
		string tmpAlignerSamStr_2;
		getline(alignerSam_ifs, tmpAlignerSamStr_2);

		int tmpReadNum = getReadNum_MPS3(tmpAlignerSamStr_1, tmpAlignerSamStr_2);
		if(tmpReadNum <= 0)
			continue;
		alignerSamInfoVec_end1[tmpReadNum - 1]->generateAlignmentInfo_MPS3(tmpAlignerSamStr_1, indexInfo);
		alignerSamInfoVec_end2[tmpReadNum - 1]->generateAlignmentInfo_MPS3(tmpAlignerSamStr_2, indexInfo);		
		alignerSamExistBoolVec_end1[tmpReadNum - 1] = true;
		alignerSamExistBoolVec_end2[tmpReadNum - 1] = true;		
	}
	alignerSam_ifs.close();
}
#endif