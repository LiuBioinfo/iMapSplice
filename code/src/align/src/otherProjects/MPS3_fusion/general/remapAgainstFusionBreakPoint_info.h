// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef REMAPAGAINSTFUSIONBREAKPOINT_INFO_H
#define REMAPAGAINSTFUSIONBREAKPOINT_INFO_H

#include "fusionBreakPointHash_info.h"

using namespace std;

class RemapAgainstFusionBreakPoint_Info
{
private:
	bool unfixedHeadOrTail2remapBool;
	bool oriGene_1_or_2_bool;

	vector<int> candiFusion_startOrEndLocInReadVec_oriGene;
	vector<int> candiFusion_breakPointPosVec_oriGene;
	vector<string> candiFusion_strandVec_oriGene;
	vector<int> candiFusion_chrNameIntVec_theOtherGene;
	vector<int> candiFusion_startPosInChrVec_theOtherGene;
	vector<int> candiFusion_breakPointPosVec_theOtherGene;
	vector<string> candiFusion_strandVec_theOtherGene;
	vector< vector<Jump_Code> > candiFusion_jumpCodeVecVec_theOtherGene;
public:
	RemapAgainstFusionBreakPoint_Info()
	{}

	int return_startOrEndLocInReadVec_oriGene(int index)
	{
		return candiFusion_startOrEndLocInReadVec_oriGene[index];
	}

	int return_breakPointPos_theOtherGene(int index)
	{
		return candiFusion_breakPointPosVec_theOtherGene[index];
	}

	int return_breakPointPos_oriGene(int index)
	{
		return candiFusion_breakPointPosVec_oriGene[index];
	}

	int return_chrNameInt_theOtherGene(int index)
	{
		return candiFusion_chrNameIntVec_theOtherGene[index];
	}

	string return_strand_theOtherGene(int index)
	{
		return candiFusion_strandVec_theOtherGene[index];
	}

	string return_strand_oriGene(int index)
	{
		return candiFusion_strandVec_oriGene[index];
	}

	string returnCandiFusionStrand(int index)
	{
		string strand_gene1;
		string strand_gene2;
		if(oriGene_1_or_2_bool)
		{	
			strand_gene1 = candiFusion_strandVec_oriGene[index];
			strand_gene2 = candiFusion_strandVec_theOtherGene[index];
		}
		else
		{
			strand_gene1 = candiFusion_strandVec_theOtherGene[index];
			strand_gene2 = candiFusion_strandVec_oriGene[index];
		}
		return strand_gene1 + strand_gene2;
	}

	int returnResultsSize()
	{
		return candiFusion_startOrEndLocInReadVec_oriGene.size();
	}

	bool setUnfixedHeadOrTail2remapBool(
		bool tmpUnfixedHeadOrTail2remapBool)
	{
		unfixedHeadOrTail2remapBool = tmpUnfixedHeadOrTail2remapBool;
	}

	bool setOriGene1or2Bool(
		bool tmp_oriGene_1_or_2_bool)
	{
		oriGene_1_or_2_bool = tmp_oriGene_1_or_2_bool;
	}

	void generateCandiFusionResultsVec(
		vector<int>& tmpCandiFusion_startOrEndLocInReadVec_oriGene,
		vector<int>& tmpCandiFusion_breakPointPosVec_oriGene,
		vector<string>& tmpCandiFusion_strandVec_oriGene,
		vector<int>& tmpCandiFusion_chrNameIntVec_theOtherGene,
		vector<int>& tmpCandiFusion_startPosInChrVec_theOtherGene,
		vector<int>& tmpCandiFusion_breakPointPosVec_theOtherGene,
		vector<string>& tmpCandiFusion_strandVec_theOtherGene,
		vector< vector<Jump_Code> >& tmpCandiFusion_jumpCodeVecVec_theOtherGene)
	{
		for(int tmp = 0; tmp < tmpCandiFusion_startOrEndLocInReadVec_oriGene.size(); tmp++)
		{
			candiFusion_startOrEndLocInReadVec_oriGene.push_back(
				tmpCandiFusion_startOrEndLocInReadVec_oriGene[tmp]);
			candiFusion_breakPointPosVec_oriGene.push_back(
				tmpCandiFusion_breakPointPosVec_oriGene[tmp]);
			candiFusion_strandVec_oriGene.push_back(
				tmpCandiFusion_strandVec_oriGene[tmp]);
			candiFusion_chrNameIntVec_theOtherGene.push_back(
				tmpCandiFusion_chrNameIntVec_theOtherGene[tmp]);
			candiFusion_startPosInChrVec_theOtherGene.push_back(
				tmpCandiFusion_startPosInChrVec_theOtherGene[tmp]);
			candiFusion_breakPointPosVec_theOtherGene.push_back(
				tmpCandiFusion_breakPointPosVec_theOtherGene[tmp]);
			candiFusion_strandVec_theOtherGene.push_back(
				tmpCandiFusion_strandVec_theOtherGene[tmp]);
			candiFusion_jumpCodeVecVec_theOtherGene.push_back(
				tmpCandiFusion_jumpCodeVecVec_theOtherGene[tmp]);			
		}
	}

	/*
	void initiateUnfixedHead2remap(
		int tmpChrNameIntInOriGene,
		int tmpLeftMostPosInChr, int tmpRightMostPosInChr,
		int tmpLeftMostLocInRead, int tmpRightMostLocInRead,
		string& tmpReadSeqfusionBreakPointHashInfo)
	{
		unfixedHeadOrTail2remapBool = true;
		chrNameIntInOriGene = tmpChrNameIntInOriGene;
		leftMostPosInChr_breakPointSearchInOriGene = tmpLeftMostPosInChr;
		rightMostPosInChr_breakPointSearchInOriGene = tmpRightMostPosInChr;
		leftMostLocInRead = tmpLeftMostLocInRead;
		rightMostLocInRead = tmpRightMostLocInRead;
		readSeq = tmpReadSeq;
		readSeqLength = tmpReadSeq.length();
	}

	void searchForCandiFusionBreakPointInOriGene(
		//vector<int>& candiBreakPointPosVecInOriGene,
		FusionBreakPointHash_Info* fusionBreakPointHashInfo
		//bool oriGene_1_or_2_bool
		)
	{
		vector<int> candiBreakPointPosVec_oriGeneAsGene1;
		vector< vector<int> > tmpChrNameIntVec_theOtherGene;
			vector<int> tmpBreakPosVec_theOtherGene;
			vector<string> tmpStrandVec_thisGene;
			vector<string> tmpStrandvec_theOtherGene;

		vector<int> candiBreakPointPosVec_oriGeneAsGene2;

		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameIntInOriGene, leftMostPosInChr_breakPointSearchInOriGene,
			rightMostPosInChr_breakPointSearchInOriGene,
			true, candiBreakPointPosVec_oriGeneAsGene1);
		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameIntInOriGene, leftMostPosInChr_breakPointSearchInOriGene,
			rightMostPosInChr_breakPointSearchInOriGene,
			false, candiBreakPointPosVec_oriGeneAsGene1);		

		for(int tmp = 0; tmp < candiBreakPointPosVec_oriGeneAsGene1.size(); tmp++)
		{
			int tmpCandiThisGeneBreakPoint = candiBreakPointPosVec_oriGeneAsGene1[tmp];
			vector<int> tmpChrNameIntVec_theOtherGene;
			vector<int> tmpBreakPosVec_theOtherGene;
			vector<string> tmpStrandVec_thisGene;
			vector<string> tmpStrandvec_theOtherGene;
			fusionBreakPointHashInfo->returnTheOtherFusionGeneBreakPoint(
				chrNameInt, tmpCandiThisGeneBreakPointPos, true,
				tmpStrandVec_thisGene, tmpChrNameIntVec_theOtherGene,
				tmpBreakPosVec_theOtherGene, tmpStrandvec_theOtherGene);			
		}

	}

	void initiateUnfixedTail2remap()
	{
		unfixedHeadOrTail2remapBool = false;
		chrNameIntInOriGene = tmpChrNameIntInOriGene;
		leftMostPosInChr_breakPointSearchInOriGene = tmpLeftMostPosInChr;
		rightMostPosInChr_breakPointSearchInOriGene = tmpRightMostPosInChr;
		leftMostLocInRead = tmpLeftMostLocInRead;
		rightMostLocInRead = tmpRightMostLocInRead;
		readSeq = tmpReadSeq;
		readSeqLength = tmpReadSeq.length();
	}*/



};
#endif