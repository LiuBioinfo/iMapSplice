// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef GENE_COUNT_VEC_H
#define GENE_COUNT_VEC_H

#include "gene_count.h"

class Gene_Count_Vec
{
private:
	vector< Gene_Count* > geneCountInfoVec;

public:
	Gene_Count_Vec()
	{}

	void initiate(int countVecNum, int geneNum)
	{
		for(int tmp = 0; tmp < countVecNum; tmp++)
		{
			Gene_Count* tmpGeneCount = new Gene_Count();
			tmpGeneCount->initiate_withGeneNum(geneNum);
			geneCountInfoVec.push_back(tmpGeneCount);
		}
	}

	double returnReadCount_inOneGene_inOneCountInfo(int indexInCountInfoVec,
		int indexInGeneVec)
	{
		return geneCountInfoVec[indexInCountInfoVec]->returnReadCount_inOneGene(
			indexInGeneVec);
	}

	int returnReadCount_unmapped_inOneCountInfo(int indexInCountInfoVec)
	{
		return geneCountInfoVec[indexInCountInfoVec]->returnReadCount_unmapped();
	}

	int returnReadCount_uniqueMapped_inOneCountInfo(int indexInCountInfoVec)
	{
		return geneCountInfoVec[indexInCountInfoVec]->returnReadCount_uniqueMapped();
	}

	int returnReadCount_multiMapped_inOneCountInfo(int indexInCountInfoVec)
	{
		return geneCountInfoVec[indexInCountInfoVec]->returnReadCount_multiMapped();
	}

	void merge2oneGeneCount(Gene_Count* mergedGeneCount)
	{
		int mergedGeneCount_geneVecSize = mergedGeneCount->returnReadCountVecSize();
		for(int tmp = 0; tmp < geneCountInfoVec.size(); tmp++)
		{
			int tmpGeneCount_geneVecSize = geneCountInfoVec[tmp]->returnReadCountVecSize();
			if(mergedGeneCount_geneVecSize != tmpGeneCount_geneVecSize)
				exit(1);
			for(int tmpIndex_inGeneVec = 0;
				tmpIndex_inGeneVec < tmpGeneCount_geneVecSize;
				tmpIndex_inGeneVec ++)
			{
				double tmpReadCount = this->returnReadCount_inOneGene_inOneCountInfo(
					tmp, tmpIndex_inGeneVec);
				mergedGeneCount->addNewReadCount(
					tmpIndex_inGeneVec, tmpReadCount);
			}
			int tmpUnmappedReadCount = this->returnReadCount_unmapped_inOneCountInfo(tmp);
			int tmpUniqueMappedReadCount = this->returnReadCount_uniqueMapped_inOneCountInfo(tmp);
			int tmpMultiMappedReadCount = this->returnReadCount_multiMapped_inOneCountInfo(tmp);
			mergedGeneCount->addNewReadCount_unmapped(tmpUnmappedReadCount);
			mergedGeneCount->addNewReadCount_uniqueMapped(tmpUniqueMappedReadCount);
			mergedGeneCount->addNewReadCount_multiMapped(tmpMultiMappedReadCount);			
		}
	}

	void unmappedReadCountIncrement_inOneCountInfo(int indexInCountInfoVec)
	{
		geneCountInfoVec[indexInCountInfoVec]->unmappedReadCountIncrement();
	}

	void uniqueMappedReadCountIncrement_inOneCountInfo(int indexInCountInfoVec)
	{
		geneCountInfoVec[indexInCountInfoVec]->uniqueMappedReadCountIncrement();
	}

	void multiMappedReadCountIncrement_inOneCountInfo(int indexInCountInfoVec)
	{
		geneCountInfoVec[indexInCountInfoVec]->multiMappedReadCountIncrement();
	}

	void addNewReadCount_inOneCountInfo(
		int indexInGeneVec, double newReadCount, int indexInCountInfoVec)
	{
		// cout << "before: geneCount[" << indexInCountInfoVec << "][" 
		// 	<< indexInGeneVec << "]: " 
		// 	<< geneCountInfoVec[indexInCountInfoVec]->returnReadCount_inOneGene(
		// 		indexInGeneVec) << endl;   
		geneCountInfoVec[indexInCountInfoVec]->addNewReadCount(
			indexInGeneVec, newReadCount);
		// cout << "after: geneCount[" << indexInCountInfoVec << "][" 
		// 	<< indexInGeneVec << "]: " 
		// 	<< geneCountInfoVec[indexInCountInfoVec]->returnReadCount_inOneGene(
		// 		indexInGeneVec) << endl;   		
	}


	void memoryFree()
	{
		for(int tmp = 0; tmp < geneCountInfoVec.size(); tmp++)
		{
			delete geneCountInfoVec[tmp];
		}		
	}	
};

#endif