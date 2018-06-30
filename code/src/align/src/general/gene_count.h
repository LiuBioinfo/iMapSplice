// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef GENE_COUNT_H
#define GENE_COUNT_H

//#include "gene_set.h"

class Gene_Count
{
private:
	vector<double> readCountVec_inGene; // size == geneInfoVec.size() in transcript_set
	int unmappedReadCount;
	int uniqueMappedReadCount;
	int multiMappedReadCount;
public:
	Gene_Count()
	{}

	int returnGeneCount(int geneIndex)
	{
		return readCountVec_inGene[geneIndex];
	}
	/*
	void output(Gene_Set* geneSetInfo, ofstream& geneCount_ofs, ofstream& stats_ofs)
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
	
		int geneNum = readCountVec_inGene.size();
		stats_ofs << "total gene num: " << geneNum << endl;	

		geneCount_ofs << "Gene ID\tcount" << endl;
		for(int tmp = 0; tmp < geneNum; tmp++)
		{
			string tmpGeneID = geneSetInfo->returnGeneID(tmp);
			double tmpGeneCount = readCountVec_inGene[tmp];
			geneCount_ofs << tmpGeneID << "\t" << tmpGeneCount << endl;
		}
	}*/

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

	double returnReadCount_inOneGene(int indexInGeneVec)
	{
		return readCountVec_inGene[indexInGeneVec];
	}

	int returnReadCountVecSize()
	{
		return readCountVec_inGene.size();
	}

	/*
	void initiate_withGeneSetInfo(Gene_Set* geneSetInfo)
	{
		int geneNum = geneSetInfo->returnGeneNum();
		initiate_withGeneNum(geneNum);
	}*/	


	void initiate_withGeneNum(int geneNum)
	{
		for(int tmp = 0; tmp < geneNum; tmp++)
		{
			readCountVec_inGene.push_back(0.0);
		}
		unmappedReadCount = 0;
		uniqueMappedReadCount = 0;
		multiMappedReadCount = 0;
	}	

	void geneReadCountIncrement(int indexInGeneVec)
	{
		readCountVec_inGene[indexInGeneVec] = 
			readCountVec_inGene[indexInGeneVec] + 1.0;
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

	void addNewReadCount(int indexInGeneVec, double newReadCount)
	{
		readCountVec_inGene[indexInGeneVec] = 
			readCountVec_inGene[indexInGeneVec] + newReadCount;
	}	
};

#endif