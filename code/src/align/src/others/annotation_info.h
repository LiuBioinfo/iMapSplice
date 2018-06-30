// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <string>
#include <string.h>
#include <map>
#include <set>

using namespace std;

typedef map<int, set<int> > AnnotatedSJ_map;
typedef AnnotatedSJ_map::iterator AnnotatedSJ_map_iter;
typedef set<int>::iterator IntSetIter;

class AnnotatedSJ_Info
{
private:
	vector< AnnotatedSJ_map* > annotatedSJ_map_vec;
public:

	AnnotatedSJ_Info(Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < (indexInfo->chrNameStr).size(); tmp++)
		{
			AnnotatedSJ_map* newAnnotatedSJ_map = new AnnotatedSJ_map;
			annotatedSJ_map_vec.push_back(newAnnotatedSJ_map);
		}
	}

	void extractSJfromSJfile(ifstream& SJ_ifs, Index_Info* indexInfo)
	{

		return;
	}

	void extractSJfromAnnotation_GAF(ifstream& annotation_ifs, Index_Info* indexInfo)
	{
		string headLine;
		getline(annotation_ifs, headLine);
		while(!annotation_ifs.eof())
		{
			vector<string> SJchrNameVec;
			vector<int> SJdonerEndPosVec;
			vector<int> SJacceptorStartPosVec;
			string tmpAnnotationEntry;
			getline(annotation_ifs, tmpAnnotationEntry);
			//cout << "SJ annotation line: " << endl << tmpAnnotationEntry << endl;
			this->parseGAF2SJ(tmpAnnotationEntry, SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec);
			//this->parseGTF2SJ(tmpAnnotationEntry, SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec);
			//this->parseGFF2SJ(tmpAnnotationEntry, SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec);

			this->insertSJt2Map_vec(SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec, indexInfo);
		}
	}

	void buildSJmap(ifstream& SJ_ifs, Index_Info* indexInfo)
	{
		//cout << "start to buildSJmap ..." << endl;
		while(!SJ_ifs.eof())
		{
			string tmpSJentry;
			getline(SJ_ifs, tmpSJentry);

			vector<string> tmpSJentryFieldStrVec;
			int searchStartLoc = 0;
			int tabLoc = 0;
			for(int tmp = 0; tmp < 3; tmp++)
			{
				tabLoc = tmpSJentry.find('\t', searchStartLoc);
				string tmpFieldStr = tmpSJentry.substr(searchStartLoc, tabLoc - 1 - searchStartLoc + 1);
				tmpSJentryFieldStrVec.push_back(tmpFieldStr);
				searchStartLoc = tabLoc + 1;
			}

			string chrName = tmpSJentryFieldStrVec[0];
			int SJ_donerEndPos = atoi(tmpSJentryFieldStrVec[1].c_str());
			int SJ_acceptorStartPos = atoi(tmpSJentryFieldStrVec[2].c_str());
			//cout << "insertSJ2Map: " << chrName << " " << SJ_donerEndPos << " " << SJ_acceptorStartPos << endl;
			this->insertSJ2Map(chrName, SJ_donerEndPos, SJ_acceptorStartPos, indexInfo);
		}		
		return;
	}

	void getTranscriptFromSJ_GAF(ifstream& annotation_ifs, Index_Info* indexInfo, ofstream& outputTranscript_ofs)
	{
		string headLine;
		getline(annotation_ifs, headLine);
		while(!annotation_ifs.eof())
		{
			vector<string> SJchrNameVec;
			vector<int> SJdonerEndPosVec;
			vector<int> SJacceptorStartPosVec;
			string tmpAnnotationEntry;
			getline(annotation_ifs, tmpAnnotationEntry);
			//cout << "SJ annotation line: " << endl << tmpAnnotationEntry << endl;
			this->parseGAF2SJ(tmpAnnotationEntry, SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec);
			//this->parseGTF2SJ(tmpAnnotationEntry, SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec);
			//this->parseGFF2SJ(tmpAnnotationEntry, SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec);

			vector<string> SJchrNameVec_found;
			vector<int> SJdonerEndPosVec_found;
			vector<int> SJacceptorStartPosVec_found;			

			//this->insertSJt2Map_vec(SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec, indexInfo);
			bool searchSJ_bool = this->searchSJ_vec(SJchrNameVec, SJdonerEndPosVec, SJacceptorStartPosVec, 
				indexInfo, SJchrNameVec_found, SJdonerEndPosVec_found, SJacceptorStartPosVec_found);

			if(searchSJ_bool)
			{
				for(int tmp = 0; tmp < SJchrNameVec_found.size(); tmp++)
				{
					outputTranscript_ofs << endl << "missingSJ: " << endl
						<< SJchrNameVec_found[tmp] << " " 
						<< SJdonerEndPosVec_found[tmp] << " " 
						<< SJacceptorStartPosVec_found[tmp] << endl;
				}
				outputTranscript_ofs << tmpAnnotationEntry << endl;
			}
		}		
	}

	void parseGAF2SJ(const string& GAFentry,
		vector<string>& SJchrNameVec, vector<int>& SJdonerEndPosVec,
		vector<int>& SJacceptorStartPosVec)
	{
		//cout << "start to parse ..." << endl;
		//string nameStrm, 
		string chromStr;//, strandStr;
		//int txStart, txEnd, cdStart, cdEnd, exonCound;
		string exonCountStr;
		string exonStartsStr, exonEndsStr;//, others;

		vector<string> fieldStrVec;
		int searchStartLoc = 0;
		int tabLoc = 0;
		for(int tmp = 0; tmp < 10; tmp++)
		{
			tabLoc = GAFentry.find('\t', searchStartLoc);
			string tmpFieldStr = GAFentry.substr(searchStartLoc, tabLoc - 1 - searchStartLoc + 1);
			fieldStrVec.push_back(tmpFieldStr);
			searchStartLoc = tabLoc + 1;
		}

		chromStr = fieldStrVec[1];
		exonCountStr = fieldStrVec[7];
		exonStartsStr = fieldStrVec[8];
		exonEndsStr = fieldStrVec[9];

		//cout << "chromStr: " << chromStr << endl 
		//	<< "exonCountStr: " << exonCountStr << endl
		//	<< "exonStartsStr: " << exonStartsStr << endl
		//	<< "exonEndsStr: " << exonEndsStr << endl;
		int exonCountInt = atoi(exonCountStr.c_str());

		int searchStartLoc_exonStart = 0;
		int commaLoc_exonStart = 0;

		vector<int> exonStartPosVec;
		for(int tmp = 0; tmp < exonCountInt; tmp++)
		{
			commaLoc_exonStart = exonStartsStr.find(',', searchStartLoc_exonStart);
			
			string tmpExonStartPosStr = exonStartsStr.substr(
				searchStartLoc_exonStart, commaLoc_exonStart - 1 - searchStartLoc_exonStart + 1);
			
			int tmpExonStartPosInt = atoi(tmpExonStartPosStr.c_str());
			exonStartPosVec.push_back(tmpExonStartPosInt);

			searchStartLoc_exonStart = commaLoc_exonStart + 1;
		}

		int searchStartLoc_exonEnd  = 0;
		int commaLoc_exonEnd = 0;

		vector<int> exonEndPosVec;
		for(int tmp = 0; tmp < exonCountInt; tmp++)
		{
			commaLoc_exonEnd = exonEndsStr.find(',', searchStartLoc_exonEnd);

			string tmpExonEndPosStr = exonEndsStr.substr(
				searchStartLoc_exonEnd, commaLoc_exonEnd - 1 - searchStartLoc_exonEnd + 1);

			int tmpExonEndPosInt = atoi(tmpExonEndPosStr.c_str());
			exonEndPosVec.push_back(tmpExonEndPosInt);

			searchStartLoc_exonEnd = commaLoc_exonEnd + 1;
		}

		int spliceJunctionCount = exonCountInt - 1;		
		for(int tmp = 0; tmp < spliceJunctionCount; tmp++)
		{
			int tmpSJdonerEndPos = exonEndPosVec[tmp];
			int tmpSJacceptorStartPos = exonStartPosVec[tmp+1];
			SJchrNameVec.push_back(chromStr);
			SJdonerEndPosVec.push_back(tmpSJdonerEndPos);
			SJacceptorStartPosVec.push_back(tmpSJacceptorStartPos);
		}
		return;
	}

	void insertSJt2Map_vec(vector<string>& SJchrNameVec, vector<int>& SJdonerEndPosVec,
		vector<int>& SJacceptorStartPosVec, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < SJchrNameVec.size(); tmp++)
		{
			this->insertSJ2Map(SJchrNameVec[tmp], SJdonerEndPosVec[tmp], 
				SJacceptorStartPosVec[tmp], indexInfo);
		}
	}

	bool searchSJ_vec(vector<string>& SJchrNameVec, vector<int>& SJdonerEndPosVec,
		vector<int>& SJacceptorStartPosVec, Index_Info* indexInfo,
		vector<string>& SJchrNameVec_found, vector<int>& SJdonerEndPosVec_found,
		vector<int>& SJacceptorStartPosVec_found)
	{
		bool SJfound_bool = false;
		for(int tmp = 0; tmp < SJchrNameVec.size(); tmp++)
		{
			bool tmpSJfound_bool = this->searchSJ(SJchrNameVec[tmp], SJdonerEndPosVec[tmp], 
				SJacceptorStartPosVec[tmp], indexInfo);
			if(tmpSJfound_bool)
			{
				SJchrNameVec_found.push_back(SJchrNameVec[tmp]);
				SJdonerEndPosVec_found.push_back(SJdonerEndPosVec[tmp]);
				SJacceptorStartPosVec_found.push_back(SJacceptorStartPosVec[tmp]);
				return true;
			}
		}		
		return false;
	}

	bool searchSJ(const string& chrName, int SJ_donerEndPos, int SJ_acceptorStartPos, Index_Info* indexInfo)
	{
		//cout << "start to search: " << chrName << " " << SJ_donerEndPos << " " << SJ_acceptorStartPos << endl;
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(chrName);
		int tmpSJ_donerEndPos = SJ_donerEndPos;
		int tmpSJ_acceptorStartPos = SJ_acceptorStartPos;		
		if((tmpSJ_chrNameInt < 0) || (tmpSJ_chrNameInt >= (indexInfo->chrNameStr).size()))
			return false;
		AnnotatedSJ_map_iter tmpAnnoMapIter 
			= annotatedSJ_map_vec[tmpSJ_chrNameInt]->find(tmpSJ_donerEndPos);		
		if(tmpAnnoMapIter != annotatedSJ_map_vec[tmpSJ_chrNameInt]->end()) // donerEndPos exists
		{
			IntSetIter tmpSetIter = (tmpAnnoMapIter->second).find(tmpSJ_acceptorStartPos);
			if(tmpSetIter != (tmpAnnoMapIter->second).end())
			{
				return true;
			}
			else
			{
				return false;	
			}
		}
		else	
		{
			return false;
		}
	}

	void insertSJ2Map(const string& chrName, int SJ_donerEndPos, int SJ_acceptorStartPos, 
		Index_Info* indexInfo)
	{
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(chrName);
		int tmpSJ_donerEndPos = SJ_donerEndPos;
		int tmpSJ_acceptorStartPos = SJ_acceptorStartPos;

		if((tmpSJ_chrNameInt < 0) || (tmpSJ_chrNameInt >= (indexInfo->chrNameStr).size()))
			return;
		AnnotatedSJ_map_iter tmpAnnoMapIter 
			= annotatedSJ_map_vec[tmpSJ_chrNameInt]->find(tmpSJ_donerEndPos);
		if(tmpAnnoMapIter != annotatedSJ_map_vec[tmpSJ_chrNameInt]->end()) // donerEndPos exists
		{
			IntSetIter tmpSetIter = (tmpAnnoMapIter->second).find(tmpSJ_acceptorStartPos);
			if(tmpSetIter != (tmpAnnoMapIter->second).end())
			{
				//cout << "error in inserting " << chrName << " " << SJ_donerEndPos << " " << SJ_acceptorStartPos << endl;
			}
			else
			{
				(tmpAnnoMapIter->second).insert(tmpSJ_acceptorStartPos);
			}
		}
		else	
		{
			set<int> tmpAcceptorStartPosSet;
			tmpAcceptorStartPosSet.insert(tmpSJ_acceptorStartPos);
			annotatedSJ_map_vec[tmpSJ_chrNameInt]->insert(
				pair<int, set<int> >(tmpSJ_donerEndPos, tmpAcceptorStartPosSet));
		}
		return;
	}

	void outputSJ(ofstream& annotationSJ_ofs, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < annotatedSJ_map_vec.size(); tmp++)
		{
			string tmpSJ_chrName = (indexInfo->chrNameStr)[tmp];
			for(AnnotatedSJ_map_iter tmpAnnoMapIter = annotatedSJ_map_vec[tmp]->begin();
				tmpAnnoMapIter != annotatedSJ_map_vec[tmp]->end(); tmpAnnoMapIter ++)
			{
				int tmpSJ_donerEndPos = tmpAnnoMapIter->first;
				for(IntSetIter tmpSetIter = (tmpAnnoMapIter->second).begin();
					tmpSetIter != (tmpAnnoMapIter->second).end(); tmpSetIter ++)
				{
					int tmpSJ_acceptorStartPos = (*tmpSetIter);
					annotationSJ_ofs << tmpSJ_chrName << "\t" 
						<< tmpSJ_donerEndPos << "\t"
						<< tmpSJ_acceptorStartPos << "\tN/A" << endl; 
				}
			}
		}
	}


};