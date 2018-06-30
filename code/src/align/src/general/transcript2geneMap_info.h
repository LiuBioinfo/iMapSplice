// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef TRANSCRIPT2GENEMAP_INFO_H
#define TRANSCRIPT2GENEMAP_INFO_H

//#include "gene_set.h"
#include "gene_count_vec.h"

class Transcript2geneMap_Info
{
private:
	vector< string > geneNameVec;

	vector< int > indexInGeneSet_transcript;
	vector< string > reissuedTranscriptIDinFileVec;
	vector< string > oriTranscriptIDinFileVec;
	vector< string > geneNameInFileVec;
public:
	Transcript2geneMap_Info()
	{

	}

	void outputGeneCount(Gene_Count* tmpGeneCountInfo, string& geneCountFile,
		string& geneCountFile_others)
	{
		ofstream geneCount_ofs(geneCountFile.c_str());
		for(int tmp = 0; tmp < geneNameVec.size(); tmp++)
		{
			string tmpGeneName = geneNameVec[tmp];
			int tmpGeneCount = tmpGeneCountInfo->returnGeneCount(tmp);
			geneCount_ofs << tmpGeneName << "\t" << tmpGeneCount << endl;
		}
		geneCount_ofs.close();

		ofstream geneCount_otherStats_ofs(geneCountFile_others.c_str());
		int readCount_unmapped = tmpGeneCountInfo->returnReadCount_unmapped();
		int readCount_uniqueMapped = tmpGeneCountInfo->returnReadCount_uniqueMapped();
		int readCount_multiMapped = tmpGeneCountInfo->returnReadCount_multiMapped();
		int readCount_total = readCount_unmapped + readCount_uniqueMapped + readCount_multiMapped;
		int readCount_unmapped_perc = ((double)readCount_unmapped / (double)readCount_total) * 100;
		int readCount_uniqueMapped_perc = ((double)readCount_uniqueMapped / (double)readCount_total) * 100;
		int readCount_multiMapped_perc = ((double)readCount_multiMapped / (double)readCount_total) * 100;
		geneCount_otherStats_ofs << "totalMapped:\t" << readCount_uniqueMapped + readCount_multiMapped
			<< " -- " << readCount_uniqueMapped_perc + readCount_multiMapped_perc << "%" << endl;
		geneCount_otherStats_ofs << "\tuniqueMapped:\t" << readCount_uniqueMapped 
			<< " -- " << readCount_uniqueMapped_perc << "%" << endl; 
		geneCount_otherStats_ofs << "\tmultiMapped:\t" << readCount_multiMapped 
			<< " -- " << readCount_multiMapped_perc << "%" << endl; 
		geneCount_otherStats_ofs << endl;
		geneCount_otherStats_ofs << "\tunmapped:\t" << readCount_unmapped
			<< " -- " << readCount_unmapped_perc << "%" << endl;
		geneCount_otherStats_ofs.close();
	}

	int returnGeneNum()
	{
		return geneNameVec.size();
	}

	int getGeneIndexFromTranscriptIndex(int indexInTranscriptVec)
	{
		return indexInGeneSet_transcript[indexInTranscriptVec];
	}

	// void initiate_andGenerateGeneSetInfo(
	// 	Gene_Set* geneSetInfo, ifstream& transcriptome_ifs,
	// 	Transcript_Set* transcriptSetInfo,
	// 	Index_Info* indexInfo, string transcript_type)
	// {
	// 	cout << "to implement" << endl;
	// 	exit(1);
	// }

	void initiate_transcript2geneMap(string& transcript2geneMap_path)
	{
		ifstream transcript2geneMap_ifs(transcript2geneMap_path.c_str());
		while(!transcript2geneMap_ifs.eof())
		{
			string tmpTranscript2geneMapStr;
			getline(transcript2geneMap_ifs, tmpTranscript2geneMapStr);
			//cout << "tmpTranscript2geneMapStr: " << tmpTranscript2geneMapStr << endl;
			if(tmpTranscript2geneMapStr == "")
				break;
			int startLoc = 0;
			int firstTabLoc = tmpTranscript2geneMapStr.find("\t", startLoc+1);
			int secondTabLoc = tmpTranscript2geneMapStr.find("\t", firstTabLoc + 1);
			int geneNameFieldStartLoc = tmpTranscript2geneMapStr.find("gene_name");
			int firstDoubleQuotaLoc = tmpTranscript2geneMapStr.find("\"", geneNameFieldStartLoc + 1);
			int secondDoubleQuotaLoc = tmpTranscript2geneMapStr.find("\"", firstDoubleQuotaLoc + 1);
			string tmpReissuedTranscriptID = tmpTranscript2geneMapStr.substr(0, firstTabLoc);
			//cout << "tmpReissuedTranscriptID: " << tmpReissuedTranscriptID << endl;
			string tmpOriTranscriptID = tmpTranscript2geneMapStr.substr(firstTabLoc + 1, secondTabLoc - 1 - firstTabLoc - 1 + 1);
			//cout << "tmpOriTranscriptID: " << tmpOriTranscriptID << endl;
			string tmpGeneName = tmpTranscript2geneMapStr.substr(firstDoubleQuotaLoc + 1, secondDoubleQuotaLoc - 1 - firstDoubleQuotaLoc - 1 + 1);
			//cout << "tmpGeneName: " << tmpGeneName << endl;
			reissuedTranscriptIDinFileVec.push_back(tmpReissuedTranscriptID);
			oriTranscriptIDinFileVec.push_back(tmpOriTranscriptID);
			geneNameInFileVec.push_back(tmpGeneName);

			bool geneNameAlreadyExistsBool = false;
			for(int tmp = geneNameVec.size() - 1; tmp >= 0; tmp--)
			{
				string tmpExistingGeneName = geneNameVec[tmp];
				if(tmpGeneName == tmpExistingGeneName)
				{
					geneNameAlreadyExistsBool = true;
					indexInGeneSet_transcript.push_back(tmp);
					break;
				}
			}
			if(!geneNameAlreadyExistsBool)
			{
				geneNameVec.push_back(tmpGeneName);
				indexInGeneSet_transcript.push_back(geneNameVec.size()-1);
			}
		}
		transcript2geneMap_ifs.close();
	}
};
#endif