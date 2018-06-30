// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SHORTANCHORSJALIGNMENT_INFO_H
#define SHORTANCHORSJALIGNMENT_INFO_H

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
#include <sstream>

using namespace std;



class ShortAnchorSJalignment_Info
{
private:
	int shortAnchorLengthMAX;
	int shortAnchorLengthMIN;
	int SJoffset;
public:
	ShortAnchorSJalignment_Info()
	{
		shortAnchorLengthMAX = 10;
		shortAnchorLengthMIN = 8;
		SJoffset = 0;
	}

	void correctShortAnchorSJalignment(	
		const string& inputSAMpath, Index_Info* indexInfo, 
		AlignmentToJunc_supportNum_Info* align2juncInfo,
		const string& output_final_sam_file,
		const string& output_corrected_sam_file_original,
		const string& output_corrected_sam_file)
	{
		ofstream final_sam_ofs(output_final_sam_file.c_str());
		ofstream corrected_sam_original_ofs(output_corrected_sam_file_original.c_str());
		ofstream corrected_sam_ofs(output_corrected_sam_file);
		ifstream sam_ifs(inputSAMpath.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			if(sam_ifs.eof())
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			int oriChrNameInt;
			bool shortAnchorSJ_head_bool, shortAnchorSJ_tail_bool;
			//vector<Jump_Code> oriJumpCode_head, oriJumpCode_tail;
			pair<int,int> shortAnchorSJposPair_head, shortAnchorSJposPair_tail;
			int headSeqLen, tailSeqLen;
			string tmpReadSeq;
			string tmpAlignmentNameStr, tmpAlignmentChrNameStr, 
				tmpAlignmentChrMapPosStr, tmpAlignmentCigarStringStr;
			string  tmpFlagStr, tmpMapQstr, tmpNextChrNameStr, tmpNextChrPosStr,
				tmpTemplateLenStr, tmpReadQualSeq, tmpAlignmentTagStr;
			this->getShortAnchorSJFromSAM_fullInfo(
				tmpAlignStr, indexInfo, oriChrNameInt,
				shortAnchorSJ_head_bool, shortAnchorSJ_tail_bool,
				//oriJumpCode_head, oriJumpCode_tail,
				shortAnchorSJposPair_head, shortAnchorSJposPair_tail,
				headSeqLen, tailSeqLen, tmpReadSeq,
				tmpAlignmentNameStr, tmpAlignmentChrNameStr,
				tmpAlignmentChrMapPosStr, tmpAlignmentCigarStringStr,
				tmpFlagStr, tmpMapQstr, tmpNextChrNameStr, tmpNextChrPosStr,
				tmpTemplateLenStr, tmpReadQualSeq, tmpAlignmentTagStr);
			if(shortAnchorSJ_head_bool)
			{		

			}
			if(shortAnchorSJ_tail_bool)
			{
				cout << endl << "******************************************" << endl;
				cout << "tmpAlignmentNameStr: " << tmpAlignmentNameStr << endl;
				cout << "tmpAlignmentChrNameStr: " << tmpAlignmentChrNameStr << endl;
				cout << "oriChrNameInt: " << oriChrNameInt << endl;
				cout << "tmpAlignmentChrMapPosStr: " << tmpAlignmentChrMapPosStr << endl;
				cout << "tmpAlignmentCigarStringStr: " << tmpAlignmentCigarStringStr << endl;
				cout << "shortAnchorSJ_head_bool: " << shortAnchorSJ_head_bool << endl;
				cout << "shortAnchorSJ_tail_bool: " << shortAnchorSJ_tail_bool << endl;
				cout << "headSeqLen: " << headSeqLen << endl;
				cout << "tailSeqLen: " << tailSeqLen << endl;
				this->correctShortAnchorTail(
					indexInfo, oriChrNameInt,
					shortAnchorSJposPair_tail, 
					tailSeqLen, tmpReadSeq,
					final_sam_ofs,
					corrected_sam_original_ofs,
					corrected_sam_ofs,
					align2juncInfo,
					tmpAlignmentNameStr, tmpAlignmentChrNameStr,
					tmpAlignmentChrMapPosStr, tmpAlignmentCigarStringStr,
					tmpFlagStr, tmpMapQstr, tmpNextChrNameStr, tmpNextChrPosStr,
					tmpTemplateLenStr, tmpReadQualSeq, tmpAlignmentTagStr);
			}
		}
		sam_ifs.close();
		final_sam_ofs.close();
		corrected_sam_original_ofs.close();
		corrected_sam_ofs.close();
	}

	void extractShortAnchorSJ_seqExtension(	
		const string& inputSAMpath, Index_Info* indexInfo, 
		AlignmentToJunc_supportNum_Info* align2juncInfo,
		const string& output_correctShortAnchorSJ_sam_file,
		const string& output_incorrectShortAnchorSJ_sam_file,
		const string& output_correctShortAnchorSJ_seqExtension_file,
		const string& output_incorrectShortAnchorSJ_seqExtension_file)
	{
		ofstream correctShortAnchorSJ_sam_ofs(
			output_correctShortAnchorSJ_sam_file.c_str());
		ofstream incorrectShortAnchorSJ_sam_ofs(
			output_incorrectShortAnchorSJ_sam_file.c_str());
		ofstream correctShortAnchorSJ_seqExtension_ofs(
			output_correctShortAnchorSJ_seqExtension_file.c_str());
		ofstream incorrectShortAnchorSJ_seqExtension_ofs(
			output_incorrectShortAnchorSJ_seqExtension_file.c_str());

		ifstream sam_ifs(inputSAMpath.c_str());
		while(1)
		{
			if(sam_ifs.eof())
				break;
			string tmpAlignStr;
			getline(sam_ifs, tmpAlignStr);
			if(sam_ifs.eof())
				break;
			if(tmpAlignStr.at(0) == '@')
				continue;
			//cout << "start to input SAM string: " << tmpAlignStr << endl;
			int oriChrNameInt;
			bool shortAnchorSJ_head_bool, shortAnchorSJ_tail_bool;
			//vector<Jump_Code> oriJumpCode_head, oriJumpCode_tail;
			pair<int,int> shortAnchorSJposPair_head, shortAnchorSJposPair_tail;
			int headSeqLen, tailSeqLen;
			string tmpReadSeq;
			string tmpAlignmentNameStr, tmpAlignmentChrNameStr, 
				tmpAlignmentChrMapPosStr, tmpAlignmentCigarStringStr;
			this->getShortAnchorSJFromSAM(
				tmpAlignStr, indexInfo, oriChrNameInt,
				shortAnchorSJ_head_bool, shortAnchorSJ_tail_bool,
				//oriJumpCode_head, oriJumpCode_tail,
				shortAnchorSJposPair_head, shortAnchorSJposPair_tail,
				headSeqLen, tailSeqLen, tmpReadSeq,
				tmpAlignmentNameStr, tmpAlignmentChrNameStr,
				tmpAlignmentChrMapPosStr, tmpAlignmentCigarStringStr);
			// cout << endl << "******************************************" << endl;
			// cout << "tmpAlignmentNameStr: " << tmpAlignmentNameStr << endl;
			// cout << "tmpAlignmentChrNameStr: " << tmpAlignmentChrNameStr << endl;
			// cout << "oriChrNameInt: " << oriChrNameInt << endl;
			// cout << "tmpAlignmentChrMapPosStr: " << tmpAlignmentChrMapPosStr << endl;
			// cout << "tmpAlignmentCigarStringStr: " << tmpAlignmentCigarStringStr << endl;
			// cout << "shortAnchorSJ_head_bool: " << shortAnchorSJ_head_bool << endl;
			// cout << "shortAnchorSJ_tail_bool: " << shortAnchorSJ_tail_bool << endl;
			// cout << "headSeqLen: " << headSeqLen << endl;
			// cout << "tailSeqLen: " << tailSeqLen << endl;
			if(shortAnchorSJ_head_bool)
			{	
				/*this->extension_output_head(
					indexInfo, oriChrNameInt,
					oriJumpCode_head, 
					shortAnchorSJposPair_head, 
					headSeqLen, tmpReadSeq,
					correctShortAnchorSJ_sam_ofs,
					incorrectShortAnchorSJ_sam_ofs,
					correctShortAnchorSJ_seqExtension_ofs,
					incorrectShortAnchorSJ_seqExtension_ofs);*/				
			}
			if(shortAnchorSJ_tail_bool)
			{
			cout << endl << "******************************************" << endl;
			cout << "tmpAlignmentNameStr: " << tmpAlignmentNameStr << endl;
			cout << "tmpAlignmentChrNameStr: " << tmpAlignmentChrNameStr << endl;
			cout << "oriChrNameInt: " << oriChrNameInt << endl;
			cout << "tmpAlignmentChrMapPosStr: " << tmpAlignmentChrMapPosStr << endl;
			cout << "tmpAlignmentCigarStringStr: " << tmpAlignmentCigarStringStr << endl;
			cout << "shortAnchorSJ_head_bool: " << shortAnchorSJ_head_bool << endl;
			cout << "shortAnchorSJ_tail_bool: " << shortAnchorSJ_tail_bool << endl;
			cout << "headSeqLen: " << headSeqLen << endl;
			cout << "tailSeqLen: " << tailSeqLen << endl;
				this->extension_output_tail(
					indexInfo, oriChrNameInt,
					//oriJumpCode_tail, 
					shortAnchorSJposPair_tail, 
					tailSeqLen, tmpReadSeq,
					correctShortAnchorSJ_sam_ofs,
					incorrectShortAnchorSJ_sam_ofs,
					correctShortAnchorSJ_seqExtension_ofs,
					incorrectShortAnchorSJ_seqExtension_ofs,
					align2juncInfo,
					tmpAlignmentNameStr, tmpAlignmentChrNameStr,
					tmpAlignmentChrMapPosStr, tmpAlignmentCigarStringStr);
			}
		}
		sam_ifs.close();
		correctShortAnchorSJ_sam_ofs.close();
		incorrectShortAnchorSJ_sam_ofs.close();
		correctShortAnchorSJ_seqExtension_ofs.close();
		incorrectShortAnchorSJ_seqExtension_ofs.close();
	}
	/*
	void extension_output_head(
		Index_Info* indexInfo, 
		int& oriChrNameInt
		vector<Jump_Code>& oriJumpCode_head,
		pair<int,int>& shortAnchorSJposPair_head, 
		int& headSeqLen, const string& tmpReadSeq,
		correctShortAnchorSJ_sam_ofs,
		incorrectShortAnchorSJ_sam_ofs,
		correctShortAnchorSJ_seqExtension_ofs,
		incorrectShortAnchorSJ_seqExtension_ofs,
		AlignmentToJunc_supportNum_Info* align2juncInfo)
	{
		int tmpSJdonerEndPos = shortAnchorSJposPair_head.first;
		int tmpSJacceptorStartPos = shortAnchorSJposPair_head.second;
		bool correctSJ_bool 
			= align2juncInfo->foundInAlignInferJunctionHash(
				oriChrNameInt, 
				tmpSJdonerEndPos, tmpSJacceptorStartPos, 
				SJoffset);
		vector<Jump_Code>& headSeqExtensionJumpCodeVec;
		FixSingleAnchor_NWDP_Info* fixSingleAnchorNwdpInfo = new FixSingleAnchor_NWDP_Info();


	}*/

	string jumpCodeVec2jumpCodeSeq(vector<Jump_Code>& tmpJumpCodeVec)
	{
		string tmpStr;
		for(int tmp = 0; tmp < tmpJumpCodeVec.size(); tmp++)
		{
			int tmpJumpCodeLen = tmpJumpCodeVec[tmp].len;
			string tmpJumpCodeType = tmpJumpCodeVec[tmp].type;
			tmpStr += int_to_str(tmpJumpCodeLen);
			tmpStr += tmpJumpCodeType;
		}
		return tmpStr;
	}

	bool toTryExtensionOrNot(int oriChrNameInt, int tmpSJdonerEndPos,
		int tmpSJacceptorStartPos, AlignmentToJunc_supportNum_Info* align2juncInfo)
	{
		int SJoffset = 5;
		int foundSJ_supNum 
			= align2juncInfo->foundInAlignmentToJuncInfo_withOffset_returnSupNum(
				oriChrNameInt, tmpSJdonerEndPos, tmpSJacceptorStartPos, SJoffset);
		cout << "foundSJ_supNum : " << foundSJ_supNum  << endl;
		if(foundSJ_supNum < 5)
			return true;
		else
			return false;
	}

	bool correctShortAnchorTail(
		Index_Info* indexInfo, 
		int oriChrNameInt,
		pair<int,int>& shortAnchorSJposPair_tail, 
		int tailSeqLen, string& tmpReadSeq,
		ofstream& final_sam_ofs,
		ofstream& corrected_sam_original_ofs,
		ofstream& corrected_sam_ofs,
		AlignmentToJunc_supportNum_Info* align2juncInfo,
		string& tmpAlignmentNameStr,
		string& tmpAlignmentChrNameStr,
		string& tmpAlignmentChrMapPosStr,
		string& tmpAlignmentCigarStringStr,
		string& tmpFlagStr, string& tmpMapQstr, 
		string& tmpNextChrNameStr, string& tmpNextChrPosStr,
		string& tmpTemplateLenStr, string& tmpReadQualSeq, 
		string& tmpAlignmentTagStr)// false -- no need to correct
	{
		int tmpSJdonerEndPos = shortAnchorSJposPair_tail.first;
		int tmpSJacceptorStartPos = shortAnchorSJposPair_tail.second;
		cout << "tmpSJdonerEndPos: " << tmpSJdonerEndPos << endl;
		cout << "tmpAcceptorStartPos: " << tmpSJacceptorStartPos << endl;
		bool toTryExtension_bool = this->toTryExtensionOrNot(oriChrNameInt, tmpSJdonerEndPos, tmpSJacceptorStartPos, align2juncInfo);
		if(toTryExtension_bool)
		{	
			vector<Jump_Code> tailSeqExtensionJumpCodeVec;
			FixSingleAnchor_NWDP_Info* fixSingleAnchorNwdpInfo 
				= new FixSingleAnchor_NWDP_Info();
			int tmpReadSeqLen = tmpReadSeq.length();
			cout << "tmpReadSeqLen: " << tmpReadSeqLen << endl;
			string tmpReadTailSeq = tmpReadSeq.substr(tmpReadSeqLen-tailSeqLen);
			cout << "tmpReadTailSeq: " << tmpReadTailSeq << endl;
			cout << "oriChrNameInt: " << oriChrNameInt << endl;
			string tmpChrSubSeq = indexInfo->returnChromStrSubstr(
				oriChrNameInt, tmpSJdonerEndPos+1, tailSeqLen+3);
			cout << "tmpChrSubSeq: " << tmpChrSubSeq << endl;
			fixSingleAnchorNwdpInfo->doNWDP(//_withMismatchJumpCode(
				tmpReadTailSeq, tmpChrSubSeq);
			fixSingleAnchorNwdpInfo->copyJumpCodeVec2TargetVec(tailSeqExtensionJumpCodeVec);
			cout << "end of fixSingleAnchorNwdpInfo " << endl;
			bool correctTailOrNot_bool = fixSingleAnchorNwdpInfo->fixedOrNot(tmpReadTailSeq.length());
			if(correctTailOrNot_bool)
			{
				corrected_sam_original_ofs << tmpAlignmentNameStr << "\t"
					<< tmpFlagStr << "\t" 
					<< tmpAlignmentChrNameStr << "\t"
					<< tmpAlignmentChrMapPosStr << "\t"
					<< tmpMapQstr << "\t"
					<< tmpAlignmentCigarStringStr << "\t"
					<< tmpNextChrNameStr << "\t"
					<< tmpNextChrPosStr << "\t"
					<< tmpTemplateLenStr << "\t" 
					<< tmpReadSeq << "\t" 
					<< tmpReadQualSeq << "\t" 
					<< tmpAlignmentTagStr
					<< endl;
				string tmpNewCigarStr = this->newCigarString_tail(
					tmpAlignmentCigarStringStr, tailSeqExtensionJumpCodeVec, 
					tmpReadSeqLen, tmpReadTailSeq.length());
				corrected_sam_ofs << tmpAlignmentNameStr << "\t"
					<< tmpFlagStr << "\t" 
					<< tmpAlignmentChrNameStr << "\t"
					<< tmpAlignmentChrMapPosStr << "\t"
					<< tmpMapQstr << "\t"
					<< tmpNewCigarStr << "\t"
					<< tmpNextChrNameStr << "\t"
					<< tmpNextChrPosStr << "\t"
					<< tmpTemplateLenStr << "\t" 
					<< tmpReadSeq << "\t" 
					<< tmpReadQualSeq << "\t" 
					<< tmpAlignmentTagStr
					<< endl;
			}
			else
			{
				final_sam_ofs << tmpAlignmentNameStr << "\t"
					<< tmpFlagStr << "\t" 
					<< tmpAlignmentChrNameStr << "\t"
					<< tmpAlignmentChrMapPosStr << "\t"
					<< tmpMapQstr << "\t"
					<< tmpAlignmentCigarStringStr << "\t"
					<< tmpNextChrNameStr << "\t"
					<< tmpNextChrPosStr << "\t"
					<< tmpTemplateLenStr << "\t" 
					<< tmpReadSeq << "\t" 
					<< tmpReadQualSeq << "\t" 
					<< tmpAlignmentTagStr
					<< endl;				
			}
			delete fixSingleAnchorNwdpInfo;
		}
		else
		{
			final_sam_ofs << tmpAlignmentNameStr << "\t"
				<< tmpFlagStr << "\t" 
				<< tmpAlignmentChrNameStr << "\t"
				<< tmpAlignmentChrMapPosStr << "\t"
				<< tmpMapQstr << "\t"
				<< tmpAlignmentCigarStringStr << "\t"
				<< tmpNextChrNameStr << "\t"
				<< tmpNextChrPosStr << "\t"
				<< tmpTemplateLenStr << "\t" 
				<< tmpReadSeq << "\t" 
				<< tmpReadQualSeq << "\t" 
				<< tmpAlignmentTagStr
				<< endl;

			return false;
		}
	}


	void extension_output_tail(
		Index_Info* indexInfo, 
		int oriChrNameInt,
		//vector<Jump_Code>& oriJumpCode_tail,
		pair<int,int>& shortAnchorSJposPair_tail, 
		int tailSeqLen, string& tmpReadSeq,
		ofstream& correctShortAnchorSJ_sam_ofs,
		ofstream& incorrectShortAnchorSJ_sam_ofs,
		ofstream& correctShortAnchorSJ_seqExtension_ofs,
		ofstream& incorrectShortAnchorSJ_seqExtension_ofs,
		AlignmentToJunc_supportNum_Info* align2juncInfo,
		string& tmpAlignmentNameStr,
		string& tmpAlignmentChrNameStr,
		string& tmpAlignmentChrMapPosStr,
		string& tmpAlignmentCigarStringStr
		)
	{
		int tmpSJdonerEndPos = shortAnchorSJposPair_tail.first;
		int tmpSJacceptorStartPos = shortAnchorSJposPair_tail.second;
		cout << "tmpSJdonerEndPos: " << tmpSJdonerEndPos << endl;
		cout << "tmpAcceptorStartPos: " << tmpSJacceptorStartPos << endl;
		bool correctSJ_bool 
			= align2juncInfo->foundInAlignmentToJuncInfo_withOffset(
				oriChrNameInt, tmpSJdonerEndPos, tmpSJacceptorStartPos, SJoffset);
		cout << "correctSJ_bool: " << correctSJ_bool << endl;
		vector<Jump_Code> tailSeqExtensionJumpCodeVec;

		FixSingleAnchor_NWDP_Info* fixSingleAnchorNwdpInfo 
			= new FixSingleAnchor_NWDP_Info();
		int tmpReadSeqLen = tmpReadSeq.length();
		cout << "tmpReadSeqLen: " << tmpReadSeqLen << endl;
		string tmpReadTailSeq = tmpReadSeq.substr(tmpReadSeqLen-tailSeqLen);
		cout << "tmpReadTailSeq: " << tmpReadTailSeq << endl;
		cout << "oriChrNameInt: " << oriChrNameInt << endl;
		string tmpChrSubSeq = indexInfo->returnChromStrSubstr(
			oriChrNameInt, tmpSJdonerEndPos+1, tailSeqLen+3);
		cout << "tmpChrSubSeq: " << tmpChrSubSeq << endl;
		fixSingleAnchorNwdpInfo->doNWDP_withMismatchJumpCode(
			tmpReadTailSeq, tmpChrSubSeq);
		fixSingleAnchorNwdpInfo->copyJumpCodeVec2TargetVec(tailSeqExtensionJumpCodeVec);
		delete fixSingleAnchorNwdpInfo;
		cout << "end of fixSingleAnchorNwdpInfo " << endl;
		if(correctSJ_bool)
		{
			correctShortAnchorSJ_sam_ofs << tmpAlignmentNameStr << "\t"
				<< tmpAlignmentChrNameStr << "\t"
				<< tmpAlignmentChrMapPosStr << "\t"
				<< tmpAlignmentCigarStringStr << "\t" 
				<< tmpReadSeq << "\t" << endl;
			correctShortAnchorSJ_seqExtension_ofs
				<< this->jumpCodeVec2jumpCodeSeq(tailSeqExtensionJumpCodeVec) << "\t"
				<< "tailSeq: " << tmpReadTailSeq << "\tLen: " << tmpReadTailSeq.length() << "\t" 
				<< tmpAlignmentNameStr << "\t"
				<< tmpAlignmentChrNameStr << "\t"
				<< tmpAlignmentChrMapPosStr << "\t"
				<< tmpAlignmentCigarStringStr << endl;
		}
		else
		{
			incorrectShortAnchorSJ_sam_ofs << tmpAlignmentNameStr << "\t"
				<< tmpAlignmentChrNameStr << "\t"
				<< tmpAlignmentChrMapPosStr << "\t"
				<< tmpAlignmentCigarStringStr << "\t" 
				<< tmpReadSeq << "\t" << endl;
			incorrectShortAnchorSJ_seqExtension_ofs
				<< this->jumpCodeVec2jumpCodeSeq(tailSeqExtensionJumpCodeVec) << "\t"
				<< "tailSeq: " << tmpReadTailSeq << "\tLen: " << tmpReadTailSeq.length()
				<< "\t" 
				<< tmpAlignmentNameStr << "\t"
				<< tmpAlignmentChrNameStr << "\t"
				<< tmpAlignmentChrMapPosStr << "\t"
				<< tmpAlignmentCigarStringStr << endl;
		}
	}

	void getShortAnchorSJFromSAM(
		const string& samStr,
		Index_Info* indexInfo, int& oriChrNameInt,
		bool& shortAnchorSJ_head_bool, bool& shortAnchorSJ_tail_bool, 
		//vector<Jump_Code>& oriJumpCode_head, vector<Jump_Code>& oriJumpCode_tail,
		pair<int,int>& shortAnchorSJposPair_head, pair<int,int>& shortAnchorSJposPair_tail,
		int& headSeqLen, int& tailSeqLen,
		string& tmpReadSeq,
		string& tmpAlignmentNameStr, string& tmpAlignmentChrNameStr,
		string& tmpAlignmentChrMapPosStr, string& tmpAlignmentCigarStringStr)
	{
		string tmpFlagStr, tmpMapQstr, tmpNextChrNameStr, tmpNextChrPosStr, 
			tmpTemplateLenStr, tmpReadQualSeq, tmpAlignmentTagStr;
		this->getShortAnchorSJFromSAM_fullInfo(samStr, indexInfo, oriChrNameInt,
			shortAnchorSJ_head_bool, shortAnchorSJ_tail_bool, 
			shortAnchorSJposPair_head, 	shortAnchorSJposPair_tail,
			headSeqLen,	tailSeqLen,
			tmpReadSeq,
			tmpAlignmentNameStr, tmpAlignmentChrNameStr,
			tmpAlignmentChrMapPosStr, tmpAlignmentCigarStringStr,
			tmpFlagStr, tmpMapQstr, tmpNextChrNameStr, tmpNextChrPosStr,
			tmpTemplateLenStr, tmpReadQualSeq, tmpAlignmentTagStr);
	}

	void getShortAnchorSJFromSAM_fullInfo(
		const string& samStr,
		Index_Info* indexInfo, int& oriChrNameInt,
		bool& shortAnchorSJ_head_bool, bool& shortAnchorSJ_tail_bool, 
		//vector<Jump_Code>& oriJumpCode_head, vector<Jump_Code>& oriJumpCode_tail,
		pair<int,int>& shortAnchorSJposPair_head, pair<int,int>& shortAnchorSJposPair_tail,
		int& headSeqLen, int& tailSeqLen,
		string& tmpReadSeq,
		string& tmpAlignmentNameStr, string& tmpAlignmentChrNameStr,
		string& tmpAlignmentChrMapPosStr, string& tmpAlignmentCigarStringStr, 
		string& tmpFlagStr, string& tmpMapQstr, string& tmpNextChrNameStr, string& tmpNextChrPosStr,
		string& tmpTemplateLenStr, string& tmpReadQualSeq, string& tmpAlignmentTagStr)
	{
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		tmpAlignmentNameStr = samFieldVec[0];
		tmpFlagStr = samFieldVec[1];
		tmpAlignmentChrNameStr = samFieldVec[2];
		oriChrNameInt = indexInfo->convertStringToInt(tmpAlignmentChrNameStr);
		tmpAlignmentChrMapPosStr = samFieldVec[3];
		int mapChrPos = atoi(tmpAlignmentChrMapPosStr.c_str());
		tmpMapQstr = samFieldVec[4];
		tmpAlignmentCigarStringStr = samFieldVec[5];
		tmpNextChrNameStr = samFieldVec[6];
		tmpNextChrPosStr = samFieldVec[7];
		tmpTemplateLenStr = samFieldVec[8];
		tmpReadSeq = samFieldVec[9];
		int tmpReadSeqLen = tmpReadSeq.length();
		tmpReadQualSeq = samFieldVec[10];
		tmpAlignmentTagStr = samStr.substr(startLoc);
		int tmpAlignmentTagStrLength = tmpAlignmentTagStr.length();
		tmpAlignmentTagStr = tmpAlignmentTagStr.substr(0, tmpAlignmentTagStrLength-1);

		vector<Jump_Code> cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(tmpAlignmentCigarStringStr, cigarStringJumpCodeVec);	

		vector< pair<int,int> > tmpSJposPairVec;
		vector< int > tmpSJindexVec_cigarStringJumpCodeVec;
		this->generateSJposVecFromJumpCodeVec(mapChrPos, cigarStringJumpCodeVec, 
			tmpSJposPairVec, tmpSJindexVec_cigarStringJumpCodeVec);

		if(tmpSJindexVec_cigarStringJumpCodeVec.size() > 0)
		{
			int firstSJindexInJumpCodeVec = tmpSJindexVec_cigarStringJumpCodeVec[0];
			int lastSJindexInJumpCodeVec = tmpSJindexVec_cigarStringJumpCodeVec[tmpSJindexVec_cigarStringJumpCodeVec.size()-1];
			int baseLocInFrontOfFirstSJ = this->getEndLocInReadOfSpecificJumpCode(
				cigarStringJumpCodeVec, firstSJindexInJumpCodeVec-1);
			int baseLocJustBehindLastSJ = this->getEndLocInReadOfSpecificJumpCode(
				cigarStringJumpCodeVec, lastSJindexInJumpCodeVec)+1;

			int tmpSeqLenBeforeFirstSJ = baseLocInFrontOfFirstSJ;
			int tmpSeqLenBehindLastSJ = tmpReadSeqLen - baseLocJustBehindLastSJ + 1;
			
			// check head
			if((tmpSeqLenBeforeFirstSJ <= shortAnchorLengthMAX)
				&&(tmpSeqLenBeforeFirstSJ >= shortAnchorLengthMIN))
			{
				shortAnchorSJ_head_bool = true;
				shortAnchorSJposPair_head.first = tmpSJposPairVec[0].first;
				shortAnchorSJposPair_head.second = tmpSJposPairVec[0].second;
				headSeqLen = tmpSeqLenBeforeFirstSJ;
			}
			else
			{
				shortAnchorSJ_head_bool = false;
			}
		
			// check tail
			if((tmpSeqLenBehindLastSJ <= shortAnchorLengthMAX)
				&&(tmpSeqLenBehindLastSJ >= shortAnchorLengthMIN))
			{
				shortAnchorSJ_tail_bool = true;
				shortAnchorSJposPair_tail.first = tmpSJposPairVec[tmpSJindexVec_cigarStringJumpCodeVec.size()-1].first;
				shortAnchorSJposPair_tail.second = tmpSJposPairVec[tmpSJindexVec_cigarStringJumpCodeVec.size()-1].second;			
				tailSeqLen = tmpSeqLenBehindLastSJ;
			}
			else
			{
				shortAnchorSJ_tail_bool = false;
			}
		}
		else
		{
			shortAnchorSJ_head_bool = false;
			shortAnchorSJ_tail_bool = false;
			return;
		}
	}

	void generateSJposVecFromJumpCodeVec(int mapChrPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
		vector< pair<int,int> >& SJposPairVec, vector<int>& SJindexVec_cigarStringJumpCodeVec)
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int lastJumpCodeIndex = tmp-1;
				int currentJumpCodeIndex = tmp;
				int tmpDonerEndPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, lastJumpCodeIndex);
				int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, currentJumpCodeIndex) + 1;
				SJposPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
				SJindexVec_cigarStringJumpCodeVec.push_back(tmp);
			}
		}		
	}	

	int getEndLocInReadOfSpecificJumpCode(vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		int tmpEndLocInRead = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "D")
			{
				//tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				//tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return tmpEndLocInRead;		
	}

	int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		int tmpEndPos = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
	}	

	void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
		
		string candidateJumpCodeType = "SMNID";
		while(1)
		{
			jumpCodeEndPosInCigarStr = 
				jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
			if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
				{break;}
			else
			{
				tmpJumpCodeLength = 
					atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	string newCigarString_tail(string& tmpAlignmentCigarStringStr, vector<Jump_Code>& tailSeqExtensionJumpCodeVec, int readLength, int tailLength)
	{
		vector<Jump_Code> finalJumpCodeVec;

		vector<Jump_Code> tmpJumpCodeVec;
		this->cigarString2jumpCodeVec(tmpAlignmentCigarStringStr, tmpJumpCodeVec);
		int otherPartLength = readLength - tailLength;

		int tmpLength = 0;
		for(int tmpIndex = 0; tmpIndex < tmpJumpCodeVec.size(); tmpIndex++)
		{
			string tmpJumpCodeType = tmpJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = tmpJumpCodeVec[tmpIndex].len;
			if((tmpJumpCodeType == "M")||(tmpJumpCodeType == "I")||(tmpJumpCodeType == "S"))			
			{
				tmpLength += tmpJumpCodeLength;
				if(tmpLength < otherPartLength)
				{
					finalJumpCodeVec.push_back(tmpJumpCodeVec[tmpIndex]);
				}
				else if(tmpLength == otherPartLength)
				{
					finalJumpCodeVec.push_back(tmpJumpCodeVec[tmpIndex]);
					break;
				}				
				else if(tmpLength > otherPartLength)
				{
					int newTmpLength = tmpLength - (tmpLength - otherPartLength);
					Jump_Code newTmpJumpCode(newTmpLength, tmpJumpCodeType);
					finalJumpCodeVec.push_back(newTmpJumpCode);
				}
				else
				{

				}
			}
		}

		int currentFinalJumpCodeVecSize = finalJumpCodeVec.size();
		string tmpLastFinalJumpCodeType = finalJumpCodeVec[currentFinalJumpCodeVecSize-1].type;
		string firstExtensionJumpCodeType = tailSeqExtensionJumpCodeVec[0].type;
		if(tmpLastFinalJumpCodeType == firstExtensionJumpCodeType)
		{
			finalJumpCodeVec[currentFinalJumpCodeVecSize-1].len += tailSeqExtensionJumpCodeVec[0].len;
		}
		else
		{
			finalJumpCodeVec.push_back(tailSeqExtensionJumpCodeVec[0]);
		}

		for(int tmp = 1; tmp < tailSeqExtensionJumpCodeVec.size(); tmp++)
		{
			finalJumpCodeVec.push_back(tailSeqExtensionJumpCodeVec[tmp]);
		}

		string newCigarString;// = this->jumpCodeVec2jumpCodeSeq
		for(int tmp = 0; tmp < finalJumpCodeVec.size(); tmp++)
		{
			newCigarString += finalJumpCodeVec[tmp].toString();
		}
		return newCigarString;
	}
};

#endif
