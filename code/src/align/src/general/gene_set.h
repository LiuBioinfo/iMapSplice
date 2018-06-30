// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef GENE_SET_H
#define GENE_SET_H

#include "gene_info.h"

class Gene_set
{
private:
	vector<Gene_Info*> geneInfoVec;
	//map<string, int> geneSetIndexMap;
public:
	Gene_set()
	{}

	int returnGeneNum()
	{
		return geneInfoVec.size();
	}

	int returnGeneInfoVecSize()
	{
		return geneInfoVec.size();
	}

	void addNewGeneInfo(Gene_Info* tmpNewGeneInfo)
	{
		geneInfoVec.push_back(tmpNewGeneInfo);
	}

	void outputGeneSetCountInfo(ofstream& geneSetCount_ofs, 
		Gene_Count* geneCountInfo, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < geneInfoVec.size(); tmp++)
		{
			double tmpGeneCount = geneCountInfo->returnGeneCount(tmp);
			string tmpGeneCountInfoStr = geneInfoVec[tmp]->returnGeneCountInfoStr(
				tmpGeneCount, indexInfo);
			geneSetCount_ofs << tmpGeneCountInfoStr << endl;
		}
	}

	void outputGeneSetCountInfo(string geneSetCount_output_path, 
		Gene_Count* geneCountInfo, Index_Info* indexInfo)
	{
		ofstream geneSetCount_ofs(geneSetCount_output_path.c_str());
		this->outputGeneSetCountInfo(geneSetCount_ofs, geneCountInfo, indexInfo);
		geneSetCount_ofs.close();
	}	

	void outputGeneSetInfo(
		ofstream& geneSet_ofs, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < geneInfoVec.size(); tmp+=)
		{
			string tmpGeneInfoStr = geneInfoVec[tmp]->returnGeneInfoStr(indexInfo);
			geneSet_ofs << tmpGeneInfoStr << endl;
		}
	}

	void outputGeneSetInfo(
		string& geneSet_output_path, Index_Info* indexInfo)
	{
		ofstream geneSet_ofs(geneSet_output_path.c_str());
		this->outputGeneSetInfo(geneSet_ofs, indexInfo)
		geneSet_ofs.close();
	}
};
#endif