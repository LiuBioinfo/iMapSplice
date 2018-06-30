// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input:  cirRNA counts files from different groups or samples;
// output: read counts for merged cirRNAs (for some cirRNAs, read counts may be 0 for some samples)

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

#include "cirRNA_transcript_info.h"
#include "../../general/splice_info.h"

using namespace std;

typedef map<string, int> CirRNAname2readCountMap;

int index_cirRNA_numCount_2Samples(bool exists_in_sample1_orNot, 
	bool exists_in_sample2_orNot)
{
	if((exists_in_sample1_orNot)&&(!exists_in_sample2_orNot))
		return 2;
	else if( (!exists_in_sample1_orNot) && exists_in_sample2_orNot)
		return 3;
	else if(exists_in_sample1_orNot && exists_in_sample2_orNot)
		return 4;
	else
	{
		cout << "not exists in both samples ! error !" << endl;
		exit(1);
	}
}

int index_cirRNA_numCount_3samples(bool exists_in_sample1_orNot, 
	bool exists_in_sample2_orNot, bool exists_in_sample3_orNot)
{
	if(exists_in_sample1_orNot && exists_in_sample2_orNot && exists_in_sample3_orNot)
		return 3;
	else if(exists_in_sample1_orNot && exists_in_sample2_orNot && (!exists_in_sample3_orNot))
		return 4;
	else if(exists_in_sample1_orNot && (!exists_in_sample2_orNot) && exists_in_sample3_orNot)
		return 5;
	else if((!exists_in_sample1_orNot) && exists_in_sample2_orNot && exists_in_sample3_orNot)
		return 6;
	else if(exists_in_sample1_orNot && (!exists_in_sample2_orNot) && (!exists_in_sample3_orNot))
		return 7;
	else if((!exists_in_sample1_orNot) && exists_in_sample2_orNot && (!exists_in_sample3_orNot))
		return 8;
	else if((!exists_in_sample1_orNot) && (!exists_in_sample2_orNot) && exists_in_sample3_orNot)
		return 9;
	else
	{
		cout << "not exists in all three samples ! error !" << endl;
		exit(1);
	}	
}

void cirRNAcountIncrement(bool exists_in_sample1_orNot, bool exists_in_sample2_orNot,
	bool exists_in_sample3_orNot, vector<int>& cirRNAcountVec, int sampleSize)
{
	if(sampleSize == 2)
	{
		int index = index_cirRNA_numCount_2Samples(exists_in_sample1_orNot,
			exists_in_sample2_orNot);
		cirRNAcountVec[index] ++;
	}
	else if(sampleSize == 3)
	{
		int index  = index_cirRNA_numCount_3samples(exists_in_sample1_orNot,
			exists_in_sample2_orNot, exists_in_sample3_orNot);
		cirRNAcountVec[index]++;
	}
	else
	{
		cout << "error in numOfsamples, for now, only supports 2 or 3 samples " << endl;
		exit(1);
	}
}

int main(int argc, char** argv)
{
	if((argc != 8)&&(argc != 11))
	{
		cout << "for now, only supports 2 samples or 3 samples, if wants to apply it for more samples,";
		cout << " please disable the function of counting cirRNA numbers in different categories (cirRNAcountVec)" << endl;
		cout << "Executable outputFolderPath name_1 cirRNA_readCount_1 totalReadNum_1 ";
		cout << " name_2 cirRNA_readCount_2 totalReadNum_2 (name_3 cirRNA_readCount_3 totalReadNum_3)" << endl;
		exit(1);
	}
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[1];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	int numOfSamples;
	int cirRNAcountVecSize;
	vector<int> cirRNAcountVec; // refer 
	if(argc == 8)
	{
		numOfSamples = 2;
		cirRNAcountVecSize = 5;
	}
	else if(argc == 11)
	{
		numOfSamples = 3;	
		cirRNAcountVecSize = 10;
	}
	else 
	{
		cout << "incorrect sample number, now only supports 2 or 3 samples, if wants to apply it for more samples, ";
		cout <<	" please disable the function of counting cirRNA numbers in different categories (cirRNAcountVec)" << endl;
		exit(1);
	}
	for(int tmp = 0; tmp < cirRNAcountVecSize; tmp++)
		cirRNAcountVec.push_back(0);

	set<string> cirRNAnameSet_merged;
	set<string> cirRNAnameSet_shared;

	vector<string> sampleNameVec;
	vector<string> cirRNAcountFileVec;
	vector<int> totalReadNumVec;
	vector<double> normalizedFactorVec;
	for(int tmpArgc = 2; tmpArgc < argc; tmpArgc+=3)
	{
		string tmpSampleName = argv[tmpArgc];
		string tmpCirRNAcountFile = argv[tmpArgc+1];
		string tmpSampleTotalReadNumStr = argv[tmpArgc+2];
		sampleNameVec.push_back(tmpSampleName);
		cirRNAcountFileVec.push_back(tmpCirRNAcountFile);
		int tmpSampleTotalReadNumInt = atoi(tmpSampleTotalReadNumStr.c_str());
		totalReadNumVec.push_back(tmpSampleTotalReadNumInt);
	}

	int firstTotalReadNum = totalReadNumVec[0];
	normalizedFactorVec.push_back(1.0);
	for(int tmp = 1; tmp < totalReadNumVec.size(); tmp++)
	{
		int tmpTotalReadNum = totalReadNumVec[tmp];
		double tmp_normalized_factor = ((double)firstTotalReadNum / (double)tmpTotalReadNum);
		normalizedFactorVec.push_back(tmp_normalized_factor);
	}

	int sampleSize = sampleNameVec.size();
	vector<CirRNAname2readCountMap> cirRNAname2readCountMapVec;
	for(int tmpSample = 0; tmpSample < sampleSize; tmpSample++)
	{
		CirRNAname2readCountMap tmpSample_cirRNAname2readCountMap;
		string tmpSample_cirRNAreadCountFile = cirRNAcountFileVec[tmpSample];
		ifstream tmpSample_cirRNAreadCount_ifs(tmpSample_cirRNAreadCountFile.c_str());
		while(!tmpSample_cirRNAreadCount_ifs.eof())
		{			
			string tmpCirRNAreadCountStr;
			getline(tmpSample_cirRNAreadCount_ifs, tmpCirRNAreadCountStr);
			if(tmpCirRNAreadCountStr == "")
			 	break;
			int tabLoc = tmpCirRNAreadCountStr.find("\t", 0);
			string tmpCirRNAnameStr = tmpCirRNAreadCountStr.substr(0, tabLoc);
			string tmpReadCountStr = tmpCirRNAreadCountStr.substr(tabLoc+1);
			int tmpReadCountInt = atoi(tmpReadCountStr.c_str());
			cirRNAcountVec[tmpSample] ++;
			tmpSample_cirRNAname2readCountMap.insert(pair<string,int>(
				tmpCirRNAnameStr, tmpReadCountInt));
			cirRNAnameSet_merged.insert(tmpCirRNAnameStr);
		}
		cirRNAname2readCountMapVec.push_back(tmpSample_cirRNAname2readCountMap);
		tmpSample_cirRNAreadCount_ifs.close();
	}

	string logFileName = outputFolderStr + "log.txt";
	ofstream log_ofs(logFileName.c_str());
	string mergedReadCountFileName = outputFolderStr + "merged.readCount";
	ofstream mergedReadCount_ofs(mergedReadCountFileName.c_str());
	string mergedReadCountFileName_normalized = outputFolderStr + "merged.normalized.readCount";
	ofstream mergedReadCount_normalized_ofs(mergedReadCountFileName_normalized.c_str());
	string mergedReadCountFileName_normalized_foldChange = outputFolderStr + "merged.normalized.readCount.foldChange";
	ofstream mergedReadCount_normalized_foldChange_ofs(mergedReadCountFileName_normalized_foldChange.c_str());	

	string mergedReadCountFileName_shared = outputFolderStr + "shared.readCount";
	ofstream mergedReadCount_shared_ofs(mergedReadCountFileName_shared.c_str());
	string mergedReadCountFileName_shared_normalized = outputFolderStr + "shared.normalized.readCount";
	ofstream mergedReadCount_shared_normalized_ofs(mergedReadCountFileName_shared_normalized.c_str());
	string mergedReadCountFileName_shared_normalized_foldChange = outputFolderStr + "shared.normalized.readCount.foldChange";
	ofstream mergedReadCount_shared_normalized_foldChange_ofs(mergedReadCountFileName_shared_normalized_foldChange.c_str());

	mergedReadCount_ofs << "CirRNA_name";
	mergedReadCount_normalized_ofs << "CirRNA_name";
	mergedReadCount_normalized_foldChange_ofs << "CirRNA_name";
	mergedReadCount_shared_ofs << "CirRNA_name";
	mergedReadCount_shared_normalized_ofs << "CirRNA_name";
	mergedReadCount_shared_normalized_foldChange_ofs << "CirRNA_name";
	for(int tmpSample = 0; tmpSample < sampleSize; tmpSample++)
	{
		mergedReadCount_ofs << "\t" << sampleNameVec[tmpSample];
		mergedReadCount_normalized_ofs << "\t" << sampleNameVec[tmpSample];
		mergedReadCount_normalized_foldChange_ofs << "\t" << sampleNameVec[tmpSample];
		mergedReadCount_shared_ofs << "\t" << sampleNameVec[tmpSample];
		mergedReadCount_shared_normalized_ofs << "\t" << sampleNameVec[tmpSample];
		mergedReadCount_shared_normalized_foldChange_ofs << "\t" << sampleNameVec[tmpSample];
	}
	for(int tmpSample = 0; tmpSample < sampleSize; tmpSample ++)
	{
		for(int tmpSample2 = tmpSample + 1; tmpSample2 < sampleSize; tmpSample2 ++)
		{
			mergedReadCount_normalized_foldChange_ofs << "\t" 
				<< sampleNameVec[tmpSample] << "/" << sampleNameVec[tmpSample2];
			mergedReadCount_shared_normalized_foldChange_ofs << "\t" 
				<< sampleNameVec[tmpSample] << "/" << sampleNameVec[tmpSample2];
		}
	}
	mergedReadCount_ofs << endl;
	mergedReadCount_normalized_ofs << endl;
	mergedReadCount_normalized_foldChange_ofs << endl;
	mergedReadCount_shared_ofs << endl;
	mergedReadCount_shared_normalized_ofs << endl;
	mergedReadCount_shared_normalized_foldChange_ofs << endl;
	for(set<string>::iterator tmpSetIter = cirRNAnameSet_merged.begin();
		tmpSetIter != cirRNAnameSet_merged.end(); tmpSetIter++)
	{
		bool exists_in_sample1_orNot;
		bool exists_in_sample2_orNot;
		bool exists_in_sample3_orNot;

		string tmpCirRNAname = (*tmpSetIter);
		mergedReadCount_ofs << tmpCirRNAname;
		mergedReadCount_normalized_ofs << tmpCirRNAname;
		mergedReadCount_normalized_foldChange_ofs << tmpCirRNAname;
		//cout << "tmpCirRNAname: " << tmpCirRNAname << endl;
		vector<int> tmpReadCountVec;
		vector<double> tmpReadCountVec_normalized;
		bool missingInSomeSampleBool = false;
		for(int tmpSample = 0; tmpSample < sampleSize; tmpSample++)
		{	
			CirRNAname2readCountMap::iterator tmpMapIter = cirRNAname2readCountMapVec[tmpSample].find(tmpCirRNAname);
			if(tmpMapIter == cirRNAname2readCountMapVec[tmpSample].end()) // not found
			{
				if(tmpSample == 0)
					exists_in_sample1_orNot = false;
				else if(tmpSample == 1)
					exists_in_sample2_orNot = false;
				else
					exists_in_sample3_orNot = false;

				mergedReadCount_ofs << "\t0";
				mergedReadCount_normalized_ofs << "\t0";
				mergedReadCount_normalized_foldChange_ofs << "\t0";
				tmpReadCountVec.push_back(0);
				tmpReadCountVec_normalized.push_back(0.0);
				missingInSomeSampleBool = true;
			}
			else // found
			{
				if(tmpSample == 0)
					exists_in_sample1_orNot = true;
				else if(tmpSample == 1)
					exists_in_sample2_orNot = true;
				else
					exists_in_sample3_orNot = true;
				int tmpCirRNAreadCount = tmpMapIter->second;
				double tmpNormalizedReadCount = tmpCirRNAreadCount * normalizedFactorVec[tmpSample];
				mergedReadCount_ofs << "\t" << tmpCirRNAreadCount;
				mergedReadCount_normalized_ofs << "\t" << tmpNormalizedReadCount;
				mergedReadCount_normalized_foldChange_ofs << "\t" << tmpNormalizedReadCount;
				tmpReadCountVec.push_back(tmpCirRNAreadCount);
				tmpReadCountVec_normalized.push_back(tmpNormalizedReadCount);
			}
		}

		cirRNAcountIncrement(exists_in_sample1_orNot, exists_in_sample2_orNot,
			exists_in_sample3_orNot, cirRNAcountVec, numOfSamples);
		
		mergedReadCount_ofs << endl;
		mergedReadCount_normalized_ofs << endl;
		if(!missingInSomeSampleBool)
		{
			mergedReadCount_shared_ofs << tmpCirRNAname;
			mergedReadCount_shared_normalized_ofs << tmpCirRNAname;
			mergedReadCount_shared_normalized_foldChange_ofs << tmpCirRNAname;
			cirRNAnameSet_shared.insert(tmpCirRNAname);
			for(int tmp = 0; tmp < tmpReadCountVec.size(); tmp++)
			{
				mergedReadCount_shared_ofs << "\t" << tmpReadCountVec[tmp];
				mergedReadCount_shared_normalized_ofs << "\t" << tmpReadCountVec_normalized[tmp];
				mergedReadCount_shared_normalized_foldChange_ofs << "\t" << tmpReadCountVec_normalized[tmp]; 
			}
			for(int tmpSample_1 = 0; tmpSample_1 < sampleSize; tmpSample_1++)
			{
				for(int tmpSample_2 = tmpSample_1 + 1; tmpSample_2 < sampleSize; tmpSample_2++)
				{
					double tmpFoldChange = (tmpReadCountVec_normalized[tmpSample_1] / tmpReadCountVec_normalized[tmpSample_2]);
					mergedReadCount_normalized_foldChange_ofs << "\t" << tmpFoldChange;
					mergedReadCount_shared_normalized_foldChange_ofs << "\t" << tmpFoldChange;
				}
			}
			mergedReadCount_normalized_foldChange_ofs << endl;		
			mergedReadCount_shared_ofs << endl;
			mergedReadCount_shared_normalized_ofs << endl;
			mergedReadCount_shared_normalized_foldChange_ofs << endl;
		}
		else
		{
			for(int tmpSample_1 = 0; tmpSample_1 < sampleSize; tmpSample_1++)
			{
				for(int tmpSample_2 = tmpSample_1 + 1; tmpSample_2 < sampleSize; tmpSample_2++)
				{
					mergedReadCount_normalized_foldChange_ofs << "\tNA";
				}
			}	
			mergedReadCount_normalized_foldChange_ofs << endl;
		}
	}
	mergedReadCount_ofs.close();
	mergedReadCount_shared_ofs.close();
	mergedReadCount_shared_normalized_ofs.close();
	mergedReadCount_shared_normalized_foldChange_ofs.close();

	for(int tmpSample = 0; tmpSample < sampleSize; tmpSample++)
	{
		string tmpFileName = outputFolderStr + sampleNameVec[tmpSample] + "_merged.readCount";
		string tmpFileName_shared = outputFolderStr + sampleNameVec[tmpSample] + "_shared" + ".readCount";
		ofstream tmpReadCount_ofs(tmpFileName.c_str());
		ofstream tmpReadCount_shared_ofs(tmpFileName_shared.c_str());
		//cout << " cirRNAnameSet_merged.size(): " << cirRNAnameSet_merged.size() << endl;
		for(set<string>::iterator tmpSetIter = cirRNAnameSet_merged.begin();
			tmpSetIter != cirRNAnameSet_merged.end(); tmpSetIter++)
		{
			string tmpCirRNAname = (*tmpSetIter);
			//cout << "tmpCirRNAname: " << tmpCirRNAname << endl;
			CirRNAname2readCountMap::iterator tmpMapIter = cirRNAname2readCountMapVec[tmpSample].find(tmpCirRNAname);
			if(tmpMapIter == cirRNAname2readCountMapVec[tmpSample].end()) // not found
			{
				//cout << "notFound !" << endl;
				tmpReadCount_ofs << tmpCirRNAname << "\t0" << endl;
			}
			else // found
			{
				//cout << "Found !" << endl;
				int tmpCirRNAreadCount = tmpMapIter->second;
				tmpReadCount_ofs << tmpCirRNAname << "\t" << tmpCirRNAreadCount << endl;
			}
		}
		tmpReadCount_ofs.close();

		for(set<string>::iterator tmpSetIter = cirRNAnameSet_shared.begin();
			tmpSetIter != cirRNAnameSet_shared.end(); tmpSetIter++)
		{
			string tmpCirRNAname = (*tmpSetIter);
			//cout << "tmpCirRNAname: " << tmpCirRNAname << endl;
			CirRNAname2readCountMap::iterator tmpMapIter = cirRNAname2readCountMapVec[tmpSample].find(tmpCirRNAname);
			if(tmpMapIter == cirRNAname2readCountMapVec[tmpSample].end()) // not found
			{
				cout << "notFound ! error ! each SJ in shared set should exist in all samples" << endl;
				tmpReadCount_shared_ofs << tmpCirRNAname << "\t0" << endl;
			}
			else // found
			{
				//cout << "Found !" << endl;
				int tmpCirRNAreadCount = tmpMapIter->second;
				tmpReadCount_shared_ofs << tmpCirRNAname << "\t" << tmpCirRNAreadCount << endl;
			}
		}
		tmpReadCount_shared_ofs.close();
	}

	if(numOfSamples == 2)
	{
		cout << "# of total cirRNA in " << sampleNameVec[0] << ": " << cirRNAcountVec[0] << endl;
		cout << "# of total cirRNA in " << sampleNameVec[1] << ": " << cirRNAcountVec[1] << endl;
		cout << endl;
		cout << "# of cirRNA only in " << sampleNameVec[0] << ": " << cirRNAcountVec[3] << endl;
		cout << "# of cirRNA only in " << sampleNameVec[1] << ": " << cirRNAcountVec[4] << endl;
		cout << endl;
		cout << "# of cirRNA in both samples: " << cirRNAcountVec[2] << endl;
	}
	else if(numOfSamples == 3)
	{
		cout << "# of total cirRNA in " << sampleNameVec[0] << ": " << cirRNAcountVec[0] << endl;
		cout << "# of total cirRNA in " << sampleNameVec[1] << ": " << cirRNAcountVec[1] << endl;
		cout << "# of total cirRNA in " << sampleNameVec[2] << ": " << cirRNAcountVec[2] << endl;
		cout << endl;
		cout << "# of cirRNA only in " << sampleNameVec[0] << ": " << cirRNAcountVec[7] << endl;
		cout << "# of cirRNA only in " << sampleNameVec[1] << ": " << cirRNAcountVec[8] << endl;
		cout << "# of cirRNA only in " << sampleNameVec[2] << ": " << cirRNAcountVec[9] << endl;
		cout << endl;
		cout << "# of cirRNA only in " << sampleNameVec[0] << " & " << sampleNameVec[1] << ": "
			<< cirRNAcountVec[4] << endl; 
		cout << "# of cirRNA only in " << sampleNameVec[0] << " & " << sampleNameVec[2] << ": "
			<< cirRNAcountVec[5] << endl;
		cout << "# of cirRNA only in " << sampleNameVec[1] << " & " << sampleNameVec[2] << ": "
			<< cirRNAcountVec[6] << endl;
		cout << endl;			
		cout << "# of cirRNA shared by all 3 samples: " << cirRNAcountVec[3] << endl;
	}
	else
	{
		cout << "incorrect sample number, now only supports 2 or 3 samples, if wants to apply it for more samples, ";
		cout << " please disable the function of counting cirRNA numbers in different categories (cirRNAcountVec)" << endl;
		exit(1);
	}

	return 0;
}