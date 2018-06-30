// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef STARCHIMERICJUNC_INFO_H
#define STARCHIMERICJUNC_INFO_H

using namespace std;

typedef map<int, int> JuncEndPos2indexMap;
typedef map<int, JuncEndPos2indexMap > JuncStartPos2endPosMap;

class StarChimericJunc_Info
{
private:
	int chrNameInt_1;
	int chrNameInt_2;
	int breakPointPos_1;
	int breakPointPos_2;
	string strand_1;
	string strand_2;
	int juncType;
	int supportNum_spanning;
	int supportNum_encompassing;

public:
	StarChimericJunc_Info()
	{}

	int returnJuncType()
	{
		return juncType;
	}

	void initiate(int tmpChrNameInt_1, int tmpChrNameInt_2,
		int tmpBreakPointPos_1, int tmpBreakPointPos_2,
		string& tmpStrand_1, string& tmpStrand_2, int tmpJuncType)
	{
		chrNameInt_1 = tmpChrNameInt_1;
		chrNameInt_2 = tmpChrNameInt_2;
		breakPointPos_1 = tmpBreakPointPos_1;
		breakPointPos_2 = tmpBreakPointPos_2;
		strand_1 = tmpStrand_1;
		strand_2 = tmpStrand_2;
		juncType = tmpJuncType;

		supportNum_spanning = 0;
		supportNum_encompassing = 0;
		if(juncType == -1)
			supportNum_encompassing = 1;
		else if(juncType >= 0)
			supportNum_spanning = 1;
		else
		{
			cout << "error in initiate in StarChimericJunc_Info" << endl;
			exit(1);
		}
	}

	void initiate(int tmpChrNameStr_1, int tmpChrNameStr_2,
		int tmpBreakPointPos_1, int tmpBreakPointPos_2,
		string& tmpStrand_1, string& tmpStrand_2, int tmpJuncType,
		Index_Info* indexInfo)
	{
		int tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
		int tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
		this->initiate(tmpChrNameInt_1, tmpChrNameInt_2,
			tmpBreakPointPos_1, tmpBreakPointPos_2,
			tmpStrand_1, tmpStrand_2, tmpJuncType);
	}	

	void supportNumIncrement_spanning()
	{
		supportNum_spanning ++;
	}

	void supportNumIncrement_encompassing()
	{
		supportNum_encompassing ++;
	}
};

class StarChimericJuncHash_Info
{
private:
	vector<StarChimericJunc_Info*> starChimericJuncInfoVec;

	vector< vector<JuncStartPos2endPosMap> > juncStartPos2endPosMapVecVec;
public:
	StarChimericJuncHash_Info()
	{}

	void output_includingEncompassed(string& outputFilePath,
		Index_Info* indexInfo) // STAR chimeric results contain encompassing juncs
	{
		ofstream junc_ofs(outputFilePath.c_str());
		for(int tmp = 0; tmp < starChimericJuncInfoVec.size(); tmp++)
		{
			string tmpJuncStr 
				= starChimericJuncInfoVec[tmp]->returnChimericJuncStr(indexInfo);
			junc_ofs << tmpJuncStr << endl;
		}
		junc_ofs.close();		
	}

	void output_onlySpanned(string& outputFilePath,
		Index_Info* indexInfo) // STAR chimeric results contain encompassing juncs
	{
		ofstream junc_ofs(outputFilePath.c_str());
		for(int tmp = 0; tmp < starChimericJuncInfoVec.size(); tmp++)
		{
			int tmpJuncType = starChimericJuncInfoVec[tmp]->returnJuncType();
			if(tmpJuncType >= 0)
			{	
				string tmpJuncStr 
					= starChimericJuncInfoVec[tmp]->returnChimericJuncStr(indexInfo);
				junc_ofs << tmpJuncStr << endl;
			}
		}
		junc_ofs.close();
	}

	void initiate(int chromNumTotal)
	{
		for(int tmpChr = 0; tmpChr < chromNumTotal; tmpChr++)
		{
			vector<JuncStartPos2endPosMap> tmpJuncStartPos2endPosMapVec;
			for(int tmpChr_2 = 0; tmpChr_2 < chromNumTotal; tmpChr_2++)
			{
				JuncStartPos2endPosMap tmpJuncStartPos2endPosMap;
				tmpJuncStartPos2endPosMapVec.push_back(tmpJuncStartPos2endPosMap);
			}
			juncStartPos2endPosMapVecVec.push_back(tmpJuncStartPos2endPosMapVec);
		}
	}

	int searchAndReturnIndexInChimericJuncInfoVec(
		int tmpChrNameInt_1, int tmpBreakPointPos_1,
		int tmpChrNameInt_2, int tmpBreakPointPos_2)
	{
		JuncStartPos2endPosMapIter::iterator tmpIter_1 
			= ((juncStartPos2endPosMapVecVec[tmpChrNameInt_1])[tmpChrNameInt_2]).find(
				tmpBreakPointPos_1);
		if(tmpIter_1 == ((juncStartPos2endPosMapVecVec[tmpChrNameInt_1])[tmpChrNameInt_2]).end())
			return -1;
		else
		{
			JuncEndPos2indexMap::iterator tmpIter_2 = (tmpIter_1->second).find(tmpBreakPointPos_2);
		}
	}

	void updateStarChimericJuncHash_withNewChimericJunc(
		int tmpChrNameInt_1, int tmpBreakPointPos_1,
		int tmpChrNameInt_2, int tmpBreakPointPos_2,
		string& tmpStrand_1, string& tmpStrand_2,
		int tmpJuncType)
	{

	}

	void generateFromStarChimericJuncOutputFile(
		string& tmpChimericJuncOutputFile, Index_Info* indexInfo)
	{
		ifstream junc_ifs(tmpChimericJuncOutputFile.c_str());
		while(!junc_ifs.eof())
		{
			string tmpJuncStr;
			getline(junc_ifs, tmpJuncStr);
			if((junc_ifs.eof())||(tmpJuncStr == ""))
				break;
			vector<string> tmpJuncFieldVec;
			int startLoc = 0;
			for(int tmp = 0; tmp < 7; tmp++)
			{
				int tabLoc = tmpJuncStr.find("\t", startLoc);
				string tmpJuncFieldStr = tmpJuncStr.substr(startLoc, tabLoc - startLoc);
				tmpJuncFieldVec.push_back(tmpJuncFieldStr);
				startLoc = tabLoc + 1;
			}
			string tmpChrNameStr_1 = tmpJuncFieldVec[0];
			int tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
			string tmpBreakPointPosStr_1 = tmpJuncFieldVec[1];
			int tmpBreakPointPos_1 = atoi(tmpBreakPointPosStr_1.c_str());
			string tmpStrand_1 = tmpJuncFieldVec[2];
			string tmpChrNameStr_2 = tmpJuncFieldVec[3];
			int tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
			string tmpBreakPointPosStr_2 = tmpJuncFieldVec[4];
			int tmpBreakPointPos_2 = atoi(tmpBreakPointPosStr_2.c_str());
			string tmpStrand_2 = tmpJuncFieldVec[5];
			string tmpJuncTypeStr = tmpJuncFieldVec[6];
			int tmpJuncType = atoi(tmpJuncTypeStr.c_str());
			this->updateStarChimericJuncHash_withNewChimericJunc(
					tmpChrNameInt_1, tmpBreakPointPos_1, 
					tmpChrNameInt_2, tmpBreakPointPos_2,
					tmpStrand_1, tmpStrand_2, tmpJuncType);
		}
		junc_ifs.close();
	}
};
#endif