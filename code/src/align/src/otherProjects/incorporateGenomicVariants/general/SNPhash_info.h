// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SNPHASH_INFO_H
#define SNPHASH_INFO_H

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

#include "SNP_info.h"

using namespace std;

typedef map<int, int> SNPinfoIndexMap;

typedef map<int, set<int> > SNPareaMap; //( areaNO = pos/100 ) intermediate hash to directly get all possible SNPs near a certain position
//typedef SNPareaHash::iterator SNPareaMapIter

class SNPhash_Info
{
private:
	vector< SNPinfoIndexMap > SNPinfoIndexMapVec;
	vector< SNPareaMap > SNPareaMapVec;

	int areaSize;
public:
	vector<SNP_Info> SNPinfoVec; 

	SNPhash_Info()
	{
		areaSize = 1000;
	}

	int returnSNP_chrNameInt(int index)
	{
		return SNPinfoVec[index].returnChrNameInt();
	}

	int returnSNP_chrPos(int index)
	{
		return SNPinfoVec[index].returnPos();
	}

	string returnSNP_referBase(int index)
	{
		return SNPinfoVec[index].returnReferBase();
	}

	string returnSNP_alterBase(int index)
	{
		return SNPinfoVec[index].returnAlterBase();
	}

	int returnSNPnum()
	{
		return SNPinfoVec.size();
	}

	void updateChromStrWithSNP(Index_Info* indexInfo)
	{

	}

	bool searchAndReturnSNPbase(int tmpChrNameInt, int tmpSNPpos, string& tmpSNPbase)
	{
		int index = this->searchAndReturnSNPinfoVecIndex(tmpChrNameInt, tmpSNPpos);
		if(index == -1)
			return false;
		else
		{
			tmpSNPbase = this->returnAlterBaseWithIndexInVec(index);
			return true;
		}
	}

	string returnAlterBaseWithIndexInVec(int index)
	{
		return SNPinfoVec[index].returnAlterBase();
	}

	int searchAndReturnSNPinfoVecIndex(int tmpChrNameInt, int tmpSNPpos)
	{
		SNPinfoIndexMap::iterator tmpSNPinfoIndexMapIter = SNPinfoIndexMapVec[tmpChrNameInt].find(tmpSNPpos);
		if(tmpSNPinfoIndexMapIter != SNPinfoIndexMapVec[tmpChrNameInt].end()) // found
			return (tmpSNPinfoIndexMapIter->second);
		else
			return -1;
	}

	bool searchSNPvecWithinRegion(int tmpChrNameInt, int tmpStartChrPos, 
		int tmpEndChrPos)
	{
		vector<int> tmpSNPposVec;
		this->returnSNPposVecWithinRegion(tmpChrNameInt, tmpStartChrPos, 
			tmpEndChrPos, tmpSNPposVec);
		if(tmpSNPposVec.size() > 0)
			return true;
		else
			return false;
	}

	void returnSNPvecWithinRegion(int tmpChrNameInt, int tmpStartChrPos, 
		int tmpEndChrPos, vector<int>& candiSNPposVec, vector<string>& candiSNPbaseVec)
	{
		int areaNOmin = (int)(tmpStartChrPos/areaSize);
		int areaNOmax = (int)(tmpEndChrPos/areaSize);
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SNPareaMap::iterator tmpSNPareaMapIter
				= SNPareaMapVec[tmpChrNameInt].find(tmpArea);
			if(tmpSNPareaMapIter != SNPareaMapVec[tmpChrNameInt].end()) // found
			{
				for(set<int>::iterator tmpIntSetIter = (tmpSNPareaMapIter->second).begin();
					tmpIntSetIter != (tmpSNPareaMapIter->second).end(); tmpIntSetIter ++)
				{
					int tmpCandiPos = (*tmpIntSetIter);
					if((tmpCandiPos >= tmpStartChrPos)&&(tmpCandiPos <= tmpEndChrPos))
					{
						candiSNPposVec.push_back(tmpCandiPos);
						string tmpSNPalterBase;
						bool searchAndReturnSNPbase_bool = this->searchAndReturnSNPbase(
							tmpChrNameInt, tmpCandiPos, tmpSNPalterBase);
						if(!searchAndReturnSNPbase_bool)
						{
							cout << "SNP not found in returnSNPvecWithinRegion ...." << endl;
							exit(1);
						}
						else
							candiSNPbaseVec.push_back(tmpSNPalterBase);
					}
				}
			}
		}
	}

	void returnSNPposVecWithinRegion(int tmpChrNameInt, int tmpStartChrPos, 
		int tmpEndChrPos, vector<int>& candiSNPposVec)
	{
		int areaNOmin = (int)(tmpStartChrPos/areaSize);
		int areaNOmax = (int)(tmpEndChrPos/areaSize);
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SNPareaMap::iterator tmpSNPareaMapIter
				= SNPareaMapVec[tmpChrNameInt].find(tmpArea);
			if(tmpSNPareaMapIter != SNPareaMapVec[tmpChrNameInt].end()) // found
			{
				for(set<int>::iterator tmpIntSetIter = (tmpSNPareaMapIter->second).begin();
					tmpIntSetIter != (tmpSNPareaMapIter->second).end(); tmpIntSetIter ++)
				{
					int tmpCandiPos = (*tmpIntSetIter);
					if((tmpCandiPos >= tmpStartChrPos)&&(tmpCandiPos <= tmpEndChrPos))
						candiSNPposVec.push_back(tmpCandiPos);
				}
			}
		}
	}

	int returnLeftMostSNPpos(int tmpChrNameInt, int tmpStartChrPos, int tmpEndChrPos)
	{
		int areaNOmin = (int)(tmpStartChrPos/areaSize);
		int areaNOmax = (int)(tmpEndChrPos/areaSize);

		vector<int> candiSNPposVec;
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SNPareaMap::iterator tmpSNPareaMapIter
				= SNPareaMapVec[tmpChrNameInt].find(tmpArea);
			if(tmpSNPareaMapIter != SNPareaMapVec[tmpChrNameInt].end()) // found
			{
				for(set<int>::iterator tmpIntSetIter = (tmpSNPareaMapIter->second).begin();
					tmpIntSetIter != (tmpSNPareaMapIter->second).end(); tmpIntSetIter ++)
				{
					int tmpCandiPos = (*tmpIntSetIter);
					if((tmpCandiPos >= tmpStartChrPos)&&(tmpCandiPos <= tmpEndChrPos))
						candiSNPposVec.push_back(tmpCandiPos);
				}
			}
			else
			{}
		}	
		//cout << "tmp tmpLeftMostSNPpos found " << endl;
		if(candiSNPposVec.size() == 0)
			return -1;
		else
		{
			int tmpLeftMostSNPpos = candiSNPposVec[0];
			for(int tmp = 1; tmp < candiSNPposVec.size(); tmp++)
			{
				int tmpSNPpos = candiSNPposVec[tmp];
				if(tmpSNPpos < tmpLeftMostSNPpos)
					tmpLeftMostSNPpos = tmpSNPpos;
			}
			return tmpLeftMostSNPpos;
		}	
	}

	int returnRightMostSNPpos(int tmpChrNameInt, int tmpStartChrPos, int tmpEndChrPos)
	{
		int areaNOmin = (int)(tmpStartChrPos/areaSize);
		int areaNOmax = (int)(tmpEndChrPos/areaSize);
		//cout << "areaNOmin: " << areaNOmin << endl;
		//cout << "areaNOmax: " << areaNOmax << endl;
		vector<int> candiSNPposVec;
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			SNPareaMap::iterator tmpSNPareaMapIter
				= SNPareaMapVec[tmpChrNameInt].find(tmpArea);
			if(tmpSNPareaMapIter != SNPareaMapVec[tmpChrNameInt].end()) // found
			{
				for(set<int>::iterator tmpIntSetIter = (tmpSNPareaMapIter->second).begin();
					tmpIntSetIter != (tmpSNPareaMapIter->second).end(); tmpIntSetIter ++)
				{
					int tmpCandiPos = (*tmpIntSetIter);
					if((tmpCandiPos >= tmpStartChrPos)&&(tmpCandiPos <= tmpEndChrPos))
						candiSNPposVec.push_back(tmpCandiPos);
				}
			}
			else
			{}
		}	
		//cout << "candiSNPposVec.size(): " << candiSNPposVec.size() << endl;
		if(candiSNPposVec.size() == 0)
			return -1;
		else
		{
			//cout << "tmp tmpRightMostSNPpos found " << endl;
			int tmpRightMostSNPpos = candiSNPposVec[0];
			for(int tmp = 1; tmp < candiSNPposVec.size(); tmp++)
			{
				int tmpSNPpos = candiSNPposVec[tmp];
				if(tmpSNPpos > tmpRightMostSNPpos)
					tmpRightMostSNPpos = tmpSNPpos;
			}
			return tmpRightMostSNPpos;
		}	
	}

	void initiate_SNPinfoIndexMapVec_SNPareaMapVec(int chromNum)
	{
		this->initiateSNPinfoIndexMapVec(chromNum);
		this->initiateSNPareaMapVec(chromNum);
	}

	void initiateSNPinfoIndexMapVec(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp ++)
		{
			SNPinfoIndexMap tmpSNPinfoIndexMap;
			SNPinfoIndexMapVec.push_back(tmpSNPinfoIndexMap);
		}
	}

	void initiateSNPareaMapVec(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp ++)
		{
			SNPareaMap tmpSNPareaMap;
			SNPareaMapVec.push_back(tmpSNPareaMap);
		}
	}

	void insertSNPpos2areaHash(int tmpChrNameInt, int tmpPos)
	{
		int tmpSNPareaNO = (int)(tmpPos/areaSize);
		SNPareaMap::iterator tmpSNPareaMapIter;
		tmpSNPareaMapIter = SNPareaMapVec[tmpChrNameInt].find(tmpSNPareaNO);
		if(tmpSNPareaMapIter == SNPareaMapVec[tmpChrNameInt].end())
		{
			set<int> newPosSet;
			newPosSet.insert(tmpPos);
			SNPareaMapVec[tmpChrNameInt].insert(pair<int, set<int> > (tmpSNPareaNO, newPosSet));
		}
		else
		{
			if((tmpSNPareaMapIter->second).find(tmpPos) == (tmpSNPareaMapIter->second).end())
				(tmpSNPareaMapIter->second).insert(tmpPos);
			else
			{}
		}
	}

	void initiateAndAddSNPinfo_fromNonDuplicateSNPlistFile(int tmpChrNameInt, int tmpPos, 
		string& tmpAlterBase, Index_Info* indexInfo) // only used when inserting SNP from SNP list 
	{
		SNPinfoIndexMap::iterator tmpSNPinfoIndexMapIter 
			= SNPinfoIndexMapVec[tmpChrNameInt].find(tmpPos);
		if(tmpSNPinfoIndexMapIter == SNPinfoIndexMapVec[tmpChrNameInt].end()) // not found
		{
			SNP_Info tmpSNPinfo;
			tmpSNPinfo.initiate(tmpChrNameInt, tmpPos, tmpAlterBase, indexInfo);
			SNPinfoVec.push_back(tmpSNPinfo);
			int tmpSNPinfoIndex = SNPinfoVec.size() - 1;
			SNPinfoIndexMapVec[tmpChrNameInt].insert(pair<int,int>(tmpPos, tmpSNPinfoIndex));
		}
		else // found 
		{}
	}

	void initiateAndAddSNPinfo2areaMap_fromNonDuplicateSNPlistFile(int tmpChrNameInt, int tmpPos)
	{
		int tmpAreaID = tmpPos / areaSize;
		SNPareaMap::iterator tmpSNPareaMapIter = SNPareaMapVec[tmpChrNameInt].find(tmpAreaID);
		if(tmpSNPareaMapIter == SNPareaMapVec[tmpChrNameInt].end()) // not found
		{
			set<int> tmpNewSet;
			tmpNewSet.insert(tmpPos);
			SNPareaMapVec[tmpChrNameInt].insert(pair<int, set<int> >(tmpAreaID, tmpNewSet));
		}
		else // found
		{
			set<int>::iterator tmpIntSetIter = (tmpSNPareaMapIter->second).find(tmpPos);
			if(tmpIntSetIter == (tmpSNPareaMapIter->second).end()) // not found
				(tmpSNPareaMapIter->second).insert(tmpPos);
			else // found
			{}
		}
	}

	void output(string& outputFilePath, Index_Info* indexInfo)
	{
		ofstream SNP_ofs(outputFilePath.c_str());
		int SNPinfoVecSize = SNPinfoVec.size();
		for(int tmp = 0; tmp < SNPinfoVecSize; tmp++)
		{
			string tmpSNPinfoStr = SNPinfoVec[tmp].returnSNPinfoStr(indexInfo);
			SNP_ofs << tmpSNPinfoStr << endl;
		}
		SNP_ofs.close();
	}

	void generateSNPhash_BeersSNPfile(string& tmpNonDuplicateSNPlistFile, Index_Info* indexInfo)
	{
		ifstream SNP_ifs(tmpNonDuplicateSNPlistFile.c_str());
		while(!SNP_ifs.eof())
		{
			string tmpSNPstr;
			getline(SNP_ifs, tmpSNPstr);
			if(tmpSNPstr == "")
				break;
			//cout << "tmpSNPstr: " << endl << tmpSNPstr << endl;
			int commaLoc = tmpSNPstr.find(":");
			string tmpChrNameStr = tmpSNPstr.substr(0,commaLoc);
			int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
			int firstTabLoc = tmpSNPstr.find("\t");
			int secondTabLoc = tmpSNPstr.find("\t", firstTabLoc + 1);
			//int thirdTabLoc = tmpSNPstr.find("\t", secondTabLoc + 1);
			string tmpPosStr = tmpSNPstr.substr(firstTabLoc+1, secondTabLoc - firstTabLoc - 1);
			int tmpPos = atoi(tmpPosStr.c_str());
			string tmpBaseChangeStr = tmpSNPstr.substr(secondTabLoc);
			//cout << "tmpBaseChangeStr: " << tmpBaseChangeStr << endl;
			string tmpRefBase = tmpBaseChangeStr.substr(1,1);
			string tmpAlterBase = tmpBaseChangeStr.substr(4,1);
			//cout << "tmpChrNameStr: " << tmpChrNameStr << endl;
			//cout << "tmpPos: " << tmpPos << endl;
			//cout << "tmpRefBase: " << tmpRefBase << endl;
			//cout << "tmpAlterBase: " << tmpAlterBase << endl;
			if(tmpChrNameInt >= 0)
			{
				this->initiateAndAddSNPinfo_fromNonDuplicateSNPlistFile(tmpChrNameInt, tmpPos, tmpAlterBase, indexInfo);
				this->initiateAndAddSNPinfo2areaMap_fromNonDuplicateSNPlistFile(tmpChrNameInt, tmpPos);
			}
			else
			{
				cout << "invalid SNPs, chrom name not in the provided prebuilt index" << endl;
				cout << tmpSNPstr << endl;
			}
		}
		for(int tmp = 0; tmp < SNPareaMapVec.size(); tmp++)
		{
			cout << "tmpChrNameInt: " << tmp << endl;
			cout << "areaSize: " << SNPareaMapVec[tmp].size() << endl;
		}
		SNP_ifs.close();
	}

	void generateSNPhash_formattedSNPfile(string& tmpFormattedSNPlistFile, Index_Info* indexInfo)
	{
		ifstream SNP_ifs(tmpFormattedSNPlistFile.c_str());
		while(!SNP_ifs.eof())
		{
			string tmpSNPstr;
			getline(SNP_ifs, tmpSNPstr);
			if(tmpSNPstr == "")
				break;
			int tabLoc_1 = tmpSNPstr.find("\t");
			int tabLoc_2 = tmpSNPstr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpSNPstr.find("\t", tabLoc_2 + 1);			
			string tmpChrNameStr = tmpSNPstr.substr(0, tabLoc_1);
			int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
			string tmpPosStr = tmpSNPstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			int tmpPos = atoi(tmpPosStr.c_str());			
			string tmpRefBase = tmpSNPstr.substr(tabLoc_2 + 1, 1);
			string tmpAlterBase = tmpSNPstr.substr(tabLoc_3 + 1, 1);
			//cout << "tmpPos: " << tmpPos << endl;
			if(tmpChrNameInt >= 0)
			{
				this->initiateAndAddSNPinfo_fromNonDuplicateSNPlistFile(tmpChrNameInt, tmpPos, tmpAlterBase, indexInfo);
				this->initiateAndAddSNPinfo2areaMap_fromNonDuplicateSNPlistFile(tmpChrNameInt, tmpPos);
			}
			else
			{
				cout << "invalid SNPs, chrom name not in the provided prebuilt index" << endl;
				cout << tmpSNPstr << endl;
			}
		}
		for(int tmp = 0; tmp < SNPareaMapVec.size(); tmp++)
		{
			cout << "tmpChrNameInt: " << tmp << endl;
			cout << "areaSize: " << SNPareaMapVec[tmp].size() << endl;
		}
		SNP_ifs.close();
	}

	void generateSNPhash_formattedSNPfile_2(string& tmpFormattedSNPlistFile, Index_Info* indexInfo)
	{
		ifstream SNP_ifs(tmpFormattedSNPlistFile.c_str());
		while(!SNP_ifs.eof())
		{
			string tmpSNPstr;
			getline(SNP_ifs, tmpSNPstr);
			if(tmpSNPstr == "")
				break;
			int tabLoc_1 = tmpSNPstr.find("\t");
			int tabLoc_2 = tmpSNPstr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpSNPstr.find("\t", tabLoc_2 + 1);			
			string tmpChrNameStr = tmpSNPstr.substr(0, tabLoc_1);
			int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
			string tmpPosStr = tmpSNPstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			int tmpPos = atoi(tmpPosStr.c_str());			
			string tmpRefBase = tmpSNPstr.substr(tabLoc_2 + 1, 1);
			string tmpAlterBase = tmpSNPstr.substr(tabLoc_3 + 1, 1);

			if(tmpChrNameInt >= 0)
			{
				this->initiateAndAddSNPinfo_fromNonDuplicateSNPlistFile(tmpChrNameInt, tmpPos, tmpAlterBase, indexInfo);
				this->initiateAndAddSNPinfo2areaMap_fromNonDuplicateSNPlistFile(tmpChrNameInt, tmpPos);
			}
			else
			{
				cout << "invalid SNPs, chrom name not in the provided prebuilt index" << endl;
				cout << tmpSNPstr << endl;
			}
		}
		for(int tmp = 0; tmp < SNPareaMapVec.size(); tmp++)
		{
			cout << "tmpChrNameInt: " << tmp << endl;
			cout << "areaSize: " << SNPareaMapVec[tmp].size() << endl;
		}
		SNP_ifs.close();
	}

	void outputSimulatedGenomReadSeq(string& outputSimulatedReadFilePath, Index_Info* indexInfo)
	{
		ofstream simuRead_ofs(outputSimulatedReadFilePath.c_str());
		int SNPinfoVecSize = SNPinfoVec.size();
		for(int tmp = 0; tmp < SNPinfoVecSize; tmp++)
		{
			string tmpSimulatedRead = SNPinfoVec[tmp].simulatedRead(indexInfo);
			simuRead_ofs << tmpSimulatedRead << endl;
		}
		simuRead_ofs.close();
	}
};
#endif