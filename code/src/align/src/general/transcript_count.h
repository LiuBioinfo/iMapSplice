// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef TRANSCRIPT_COUNT_H
#define TRANSCRIPT_COUNT_H

#include "transcript_set.h"

class Transcript_Count
{
private:
	vector<double> readCountVec_inTranscript; // size == normal_transcript_vec.size() in transcript_set
	int unmappedReadCount;
	int uniqueMappedReadCount;
	int multiMappedReadCount;
public:
	Transcript_Count()
	{}

	void output(Transcript_Set* transcriptInfo, ofstream& transcriptCount_ofs,
		ofstream& stats_ofs)
	{
		int mappedReadNum = uniqueMappedReadCount + multiMappedReadCount;
		int totalReadNum = mappedReadNum + unmappedReadCount;
		double unmappedRead_perc = ((double)unmappedReadCount/(double)totalReadNum)*100;
		double uniqueMappedRead_perc = ((double)uniqueMappedReadCount/(double)totalReadNum)*100;
		double multiMappedRead_perc = ((double)multiMappedReadCount/(double)totalReadNum)*100;
		double mappedRead_perc = ((double)mappedReadNum/(double)totalReadNum)*100;
		stats_ofs << "Total read num: " << totalReadNum << endl;
		stats_ofs << "mappedReadNum: " << mappedReadNum << " -- " << mappedRead_perc << "%" << endl;
		stats_ofs << "\tunqiueMappedReadNum: " << uniqueMappedReadCount << " -- " << uniqueMappedRead_perc << "%" << endl;
		stats_ofs << "\tmultiMappedReadNum: " << multiMappedReadCount << " -- " << multiMappedRead_perc << "%" << endl;
		stats_ofs << "unmappedReadNum: " << unmappedReadCount << " -- " << unmappedRead_perc << "%" << endl;

		int transcriptNum = readCountVec_inTranscript.size();
		stats_ofs << "total transcript num: " << transcriptNum << endl;	

		transcriptCount_ofs << "Transcript ID\tcount" << endl;
		for(int tmp = 0; tmp < transcriptNum; tmp++)
		{
			string tmpTranscriptID = transcriptInfo->returnTranscriptID(tmp);
			double tmpTranscriptCount = readCountVec_inTranscript[tmp];
			transcriptCount_ofs << tmpTranscriptID << "\t" << tmpTranscriptCount << endl;
		}
	}

	int returnReadCount_unmapped()
	{
		return unmappedReadCount;
	}

	int returnReadCount_uniqueMapped()
	{
		return uniqueMappedReadCount;
	}

	int returnReadCount_multiMapped()
	{
		return multiMappedReadCount;
	}

	double returnReadCount_inOneTranscript(int indexInTranscriptVec)
	{
		return readCountVec_inTranscript[indexInTranscriptVec];
	}

	int returnReadCountVecSize()
	{
		return readCountVec_inTranscript.size();
	}

	void initiate_withTranscriptSetInfo(Transcript_Set* transcriptSetInfo)
	{
		int transcriptNum = transcriptSetInfo->returnTranscriptNum();
		initiate_withTranscriptNum(transcriptNum);
	}

	void initiate_withTranscriptNum(int transcriptNum)
	{
		for(int tmp = 0; tmp < transcriptNum; tmp++)
		{
			readCountVec_inTranscript.push_back(0.0);
		}
		unmappedReadCount = 0;
		uniqueMappedReadCount = 0;
		multiMappedReadCount = 0;
	}

	void transcriptReadCountIncrement(int indexInTranscriptVec)
	{
		readCountVec_inTranscript[indexInTranscriptVec] = 
			readCountVec_inTranscript[indexInTranscriptVec] + 1.0;
	}

	void unmappedReadCountIncrement()
	{
		unmappedReadCount ++;
	}

	void uniqueMappedReadCountIncrement()
	{
		uniqueMappedReadCount ++;
	}

	void multiMappedReadCountIncrement()
	{
		multiMappedReadCount ++;
	}	

	void addNewReadCount_unmapped(int newReadCount)
	{
		unmappedReadCount += newReadCount;
	}

	void addNewReadCount_uniqueMapped(int newReadCount)
	{
		uniqueMappedReadCount += newReadCount;
	}

	void addNewReadCount_multiMapped(int newReadCount)
	{
		multiMappedReadCount += newReadCount;
	}

	void addNewReadCount(int indexInTranscriptVec, double newReadCount)
	{
		readCountVec_inTranscript[indexInTranscriptVec] = 
			readCountVec_inTranscript[indexInTranscriptVec] + newReadCount;
	}
};

#endif