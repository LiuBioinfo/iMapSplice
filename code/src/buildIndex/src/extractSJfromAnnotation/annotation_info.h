// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ANNOTATION_INFO_H
#define ANNOTATION_INFO_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include "transcript_set.h"

//typedef map<int, set<int,int> > SpliceEndHash;
typedef map<int, set<int> > PosHash;
typedef map<int, set<int> > SpliceAreaStartHash;
typedef PosHash::iterator PosHashIter;
typedef SpliceAreaStartHash::iterator AreaHashIter;

class Annotation_Info
{
private:
	vector< PosHash > SJposHash_vec;
 	int areaSize;
	vector< SpliceAreaStartHash > SJstartPosAreaHash_vec;
public:
	Annotation_Info()
	{
		areaSize = 1000;
	}

	void outputSJ(ofstream& SJ_ofs, Index_Info* indexInfo)
	{
		int tmpJuncNum = 0;
		for(int tmpPosHashIndex = 0; tmpPosHashIndex < SJposHash_vec.size(); tmpPosHashIndex++)
		{
			string tmpChrNameStr = indexInfo->returnChromName(tmpPosHashIndex);
			for(PosHashIter tmpPosHashIter = SJposHash_vec[tmpPosHashIndex].begin();
					tmpPosHashIter != SJposHash_vec[tmpPosHashIndex].end();
					tmpPosHashIter ++)
			{
				int tmpSJdonerSitePos = tmpPosHashIter->first; 
				for(set<int>::iterator intSetIter = (tmpPosHashIter->second).begin();
					intSetIter != (tmpPosHashIter->second).end(); intSetIter++)
				{
					int tmpSJacceptorSitePos = (*intSetIter);
					SJ_ofs << tmpChrNameStr << "\t" << tmpSJdonerSitePos << "\t" << tmpSJacceptorSitePos
						<< "\t" << "Junc_" << tmpJuncNum << endl;
					tmpJuncNum ++;
				}
			}
		}
	}

	bool SJfoundInAnnotation(int chrNameInt, int chromPos_doner, int chromPos_acceptor)
	{
		PosHashIter tmpSJhashIter = SJposHash_vec[chrNameInt].find(chromPos_doner);
		if(tmpSJhashIter == SJposHash_vec[chrNameInt].end())
		{
			return false;
		}
		else
		{
			set<int>::iterator tmpPosSetIter = (tmpSJhashIter->second).find(chromPos_acceptor);
			if(tmpPosSetIter == (tmpSJhashIter->second).end())
			{
				return false;
			}
			else
				return true;
		}
	}

	void initiateAndReadAnnotationFile(Index_Info* indexInfo, 
		ifstream& annotation_ifs)
	{
		int chromNum = indexInfo->returnChromNum();
		this->initiateHash_vec(chromNum);

		string headLine;
		getline(annotation_ifs, headLine);
		while(!annotation_ifs.eof())
		{
			string tmpAnnotationEntry;
			getline(annotation_ifs, tmpAnnotationEntry);
			//if(tmpAnnotationEntry == "")
			//	break;
			string entryString;
			getline(annotation_ifs, entryString);
			if(annotation_ifs.eof())
				break;
			int tabLocation1 = entryString.find('\t', 0);
			int tabLocation2 = entryString.find('\t', tabLocation1+1);
			int tabLocation3 = entryString.find('\t', tabLocation2+1);
			string chrIntString = entryString.substr(0, tabLocation1);
			string spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
			string spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
			int chrInt = indexInfo->convertStringToInt(chrIntString);
			int spliceStartPos = atoi(spliceStartPosString.c_str());
			int spliceEndPos = atoi(spliceEndPosString.c_str());
			this->insert2AreaAndPosHash(chrInt, spliceStartPos, spliceEndPos);			
		}
	}

	void insert2AreaHash(int chrInt, int spliceStartPos)
	{
		int spliceStartPosAreaNO = (int)(spliceStartPos/areaSize);
		//int spliceEndPosAreaNO = (int)(spliceEndPos/areaSize);
		AreaHashIter foundAreaHashIter;
		foundAreaHashIter = SJstartPosAreaHash_vec[chrInt].find(spliceStartPosAreaNO);
		if(foundAreaHashIter == SJstartPosAreaHash_vec[chrInt].end()) // area not found 
		{
			set<int> newSet;
			newSet.insert(spliceStartPos);
			SJstartPosAreaHash_vec[chrInt].insert(pair<int, set<int> > (spliceStartPosAreaNO, newSet));
		}
		else
		{
			set<int>::iterator intSetIter;
			intSetIter = (foundAreaHashIter->second).find(spliceStartPos);
			if(intSetIter == (foundAreaHashIter->second).end()) // spliceStartPos not found
			{
				(foundAreaHashIter->second).insert(spliceStartPos);
			}
		}
	}

	void insert2PosHash(int chrInt, int spliceStartPos, int spliceEndPos)
	{
		PosHashIter foundPosHashIter;
		foundPosHashIter = SJposHash_vec[chrInt].find(spliceStartPos);
		if(foundPosHashIter == SJposHash_vec[chrInt].end()) // spliceStartPos not found
		{
			set<int> newSet;
			newSet.insert(spliceEndPos);
			SJposHash_vec[chrInt].insert(pair<int, set<int> > (spliceStartPos, newSet));
		}
		else
		{
			set<int>::iterator intSetIter;
			intSetIter = (foundPosHashIter->second).find(spliceEndPos);
			if(intSetIter == (foundPosHashIter->second).end())
			{
				(foundPosHashIter->second).insert(spliceEndPos);
			}			
		}
	}

	void insert2AreaAndPosHash(int chrInt, int spliceStartPos, int spliceEndPos)
	{
		this->insert2AreaHash(chrInt, spliceStartPos);
		this->insert2PosHash(chrInt, spliceStartPos, spliceEndPos);
	}

	void initiateHash_vec(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp ++)
		{
			SpliceAreaStartHash newSJareaStartHash;
			SJstartPosAreaHash_vec.push_back(newSJareaStartHash);
		}
		for(int tmp = 0; tmp < chromNum; tmp ++)
		{
			PosHash newSJposHash;
			SJposHash_vec.push_back(newSJposHash);
		}
	}

	void initiateHash_vec(Index_Info* indexInfo)
	{
		int chromNum = indexInfo->returnChromNum();
		this->initiateHash_vec(chromNum);		
	}

	void extractSJfromTranscriptSet(Transcript_Set* transcriptSet, const string& transcriptType, Index_Info* indexInfo)
	{
		if(transcriptType == "GTF")
		{
			this->extractSJ_GTF(transcriptSet, indexInfo);
		}
		else
		{
			cout << "invalid transcriptType in extractSJfromTranscriptSet of annotation_info.h" << endl;
			exit(1);
		}
	}

	void extractSJ_GTF(Transcript_Set* transcriptSet, Index_Info* indexInfo)
	{
		int transcript_size = transcriptSet->returnTranscriptNum();
		for(int tmp = 0; tmp < transcript_size; tmp++)
		{
			//cout << "tmpTranscript_size: " << tmp << endl;
			string tmpSJchrNameStr = transcriptSet->returnChromName(tmp);
			int tmpChrNameInt = indexInfo->convertStringToInt(tmpSJchrNameStr);
			vector< pair<int,int> > tmpSJvec;
			transcriptSet->getSJsitePairVec(tmp, tmpSJvec);
			for(int tmpSJ = 0; tmpSJ < tmpSJvec.size(); tmpSJ++)
			{
				//cout << "tmpSJ: " << tmpSJ << endl;
				int tmpSJdonerSite = tmpSJvec[tmpSJ].first;
				int tmpSJacceptorSite = tmpSJvec[tmpSJ].second;	
				//cout << "SJ: " << tmpChrNameInt << "\t" << tmpSJdonerSite << "\t" << tmpSJacceptorSite << endl;
				this->insert2AreaAndPosHash(tmpChrNameInt, tmpSJdonerSite, tmpSJacceptorSite);
			}
		}
	}
};

#endif