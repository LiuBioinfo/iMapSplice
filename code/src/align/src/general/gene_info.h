// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef GENE_INFO_H
#define GENE_INFO_H

#include "transcript_set.h"

class Gene_Info
{
private:
	vector<int> transcriptInfoIndexVecInTranscriptSet;
	string geneName;
public:
	Gene_Info()
	{}

	void initiateGeneInfo(int tmpTranscriptInfoIndexInTranscriptSet,
		string& tmpGeneName)
	{
		transcriptInfoIndexVecInTranscriptSet.push_back(
			tmpTranscriptInfoIndexInTranscriptSet);
		geneName = tmpGeneName;
	}

	void addNewTranscript(int tmpNewTranscriptIndex)
	{
		transcriptInfoIndexVecInTranscriptSet.push_back(
			tmpNewTranscriptIndex);
	}

	string returnGeneName()
	{
		return geneName;
	}

	string returnGeneID()
	{
		return geneName;	
	}

	string returnTranscriptID(Transcript_Set* transcriptSetInfo, int transcriptIndexInGene)
	{
		int transcriptIndexInTranscriptSet 
			= transcriptInfoIndexVecInTranscriptSet[transcriptIndexInGene];
		return transcriptSetInfo->returnTranscriptID(transcriptIndexInTranscriptSet);
	}

	string returnChromName(Transcript_Set* transcriptSetInfo)
	{
		int firstTranscriptIndexInTranscriptSet
			= transcriptInfoIndexVecInTranscriptSet[0];
		return transcriptSetInfo->returnChromName(firstTranscriptIndexInTranscriptSet);
	}
};
#endif