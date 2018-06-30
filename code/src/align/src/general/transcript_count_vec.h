// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef TRANSCRIPT_COUNT_VEC_H
#define TRANSCRIPT_COUNT_VEC_H

#include "transcript_count.h"

class Transcript_Count_Vec
{
private:
	vector<Transcript_Count*> transcriptCountInfoVec;

public:
	Transcript_Count_Vec()
	{

	}

	void initiate(int countVecNum, Transcript_Set* transcriptSetInfo)// countVecNum = threadsNum
	{
		int transcriptNum = transcriptSetInfo->returnTranscriptNum();
		for(int tmp = 0; tmp < countVecNum; tmp++)
		{
			Transcript_Count* tmpTranscriptCount = new Transcript_Count();
			tmpTranscriptCount->initiate_withTranscriptNum(transcriptNum);
			transcriptCountInfoVec.push_back(tmpTranscriptCount);
		}
	}

	double returnReadCount_inOneTranscript_inOneCountInfo(int indexInCountInfoVec, 
		int indexInTranscriptVec)
	{
		return transcriptCountInfoVec[indexInCountInfoVec]->returnReadCount_inOneTranscript(
			indexInTranscriptVec);
	}

	int returnReadCount_unmapped_inOneCountInfo(int indexInCountInfoVec)
	{
		return transcriptCountInfoVec[indexInCountInfoVec]->returnReadCount_unmapped();
	}

	int returnReadCount_uniqueMapped_inOneCountInfo(int indexInCountInfoVec)
	{
		return transcriptCountInfoVec[indexInCountInfoVec]->returnReadCount_uniqueMapped();
	}

	int returnReadCount_multiMapped_inOneCountInfo(int indexInCountInfoVec)
	{
		return transcriptCountInfoVec[indexInCountInfoVec]->returnReadCount_multiMapped();
	}

	void merge2oneTranscriptCount(Transcript_Count* mergedTranscriptCount)
	{
		int mergedTranscriptCount_transcriptVecSize = mergedTranscriptCount->returnReadCountVecSize();
		// cout << "start to merge" << endl;
		// cout << "transcriptCountInfoVec.size(): " << transcriptCountInfoVec.size() << endl;
		for(int tmp = 0; tmp < transcriptCountInfoVec.size(); tmp++)
		{
			int tmpTranscriptCount_transcriptVecSize = transcriptCountInfoVec[tmp]->returnReadCountVecSize();
			//cout << "tmpTranscriptCount_transcriptVecSize: " << tmpTranscriptCount_transcriptVecSize << endl;
			if(mergedTranscriptCount_transcriptVecSize != tmpTranscriptCount_transcriptVecSize)
			{
				//cout << "inconsistent transcriptVecSize between merged one and the one in countInfoVec: " << tmp << endl;
				exit(1);
			}
			for(int tmpIndex_inTranscriptVec = 0; 
				tmpIndex_inTranscriptVec < tmpTranscriptCount_transcriptVecSize;
				tmpIndex_inTranscriptVec ++)
			{
				double tmpReadCount = this->returnReadCount_inOneTranscript_inOneCountInfo(
					tmp, tmpIndex_inTranscriptVec);
				// if(tmpReadCount != 0)
				// {
				// 	cout << "counts do exist: " << tmpReadCount << endl;
				// }
				mergedTranscriptCount->addNewReadCount(
					tmpIndex_inTranscriptVec, tmpReadCount);
			}
			int tmpUnmappedReadCount = this->returnReadCount_unmapped_inOneCountInfo(tmp);
			int tmpUniqueMappedReadCount = this->returnReadCount_uniqueMapped_inOneCountInfo(tmp);
			int tmpMultiMappedReadCount = this->returnReadCount_multiMapped_inOneCountInfo(tmp);
			mergedTranscriptCount->addNewReadCount_unmapped(tmpUnmappedReadCount);
			mergedTranscriptCount->addNewReadCount_uniqueMapped(tmpUniqueMappedReadCount);
			mergedTranscriptCount->addNewReadCount_multiMapped(tmpMultiMappedReadCount);
		}
	}

	void unmappedReadCountIncrement_inOneCountInfo(int indexInCountInfoVec)
	{
		transcriptCountInfoVec[indexInCountInfoVec]->unmappedReadCountIncrement();
	}

	void uniqueMappedReadCountIncrement_inOneCountInfo(int indexInCountInfoVec)
	{
		transcriptCountInfoVec[indexInCountInfoVec]->uniqueMappedReadCountIncrement();
	}

	void multiMappedReadCountIncrement_inOneCountInfo(int indexInCountInfoVec)
	{
		transcriptCountInfoVec[indexInCountInfoVec]->multiMappedReadCountIncrement();
	}

	void addNewReadCount_inOneCountInfo(
		int indexInTranscriptVec, double newReadCount, int indexInCountInfoVec)
	{
		// cout << "before: transcriptCount[" << indexInCountInfoVec << "][" 
		// 	<< indexInTranscriptVec << "]: " 
		// 	<< transcriptCountInfoVec[indexInCountInfoVec]->returnReadCount_inOneTranscript(
		// 		indexInTranscriptVec) << endl;   
		transcriptCountInfoVec[indexInCountInfoVec]->addNewReadCount(
			indexInTranscriptVec, newReadCount);
		// cout << "after: transcriptCount[" << indexInCountInfoVec << "][" 
		// 	<< indexInTranscriptVec << "]: " 
		// 	<< transcriptCountInfoVec[indexInCountInfoVec]->returnReadCount_inOneTranscript(
		// 		indexInTranscriptVec) << endl;   		
	}


	void memoryFree()
	{
		for(int tmp = 0; tmp < transcriptCountInfoVec.size(); tmp++)
		{
			delete transcriptCountInfoVec[tmp];
		}		
	}
};

#endif