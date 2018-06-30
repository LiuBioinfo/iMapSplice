// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SNPINTRANSCRIPT_INFO_H
#define SNPINTRANSCRIPT_INFO_H

#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "SNPhash_info.h"

using namespace std;

class SNPinTranscript_Info
{
private:
	vector< vector<int> > SNPlocInTranscriptVecVec; 
	vector< vector<string> > SNPalterBaseInTranscriptVecVec;
public:
	SNPinTranscript_Info()
	{}

	string jumpCodeVec2Str(vector<Jump_Code>& jumpCodeVec)
	{
		string tmpStr;
		for(int tmp = 0; tmp < jumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeStr = (jumpCodeVec[tmp]).toString();
			tmpStr += tmpJumpCodeStr;
		}
		return tmpStr;
	}

	void addSNP2transcript(Transcript_Set* transcriptInfo, 
		SNPhash_Info& snpHashInfo, Index_Info* indexInfo)
	{
		//cout << "start to addSNP2transcript " << endl;
		int transcriptNum = transcriptInfo->returnTranscriptNum();
		//cout << "transcriptNum: " << transcriptNum << endl;
		for(int tmp = 0; tmp < transcriptNum; tmp++)
		{
			//cout << "tmp transcript: " << tmp << endl;
			vector< vector<int> > tmpCandiSNPposVecVecInTranscript;
			vector< vector<string> > tmpCandiSNPbaseVecVecInTranscript;
			string tmpTranscriptChromName = transcriptInfo->returnTranscriptChromName(tmp);
			int tmpTranscriptChromNameInt = indexInfo->convertStringToInt(tmpTranscriptChromName);
			vector<int> tmpTranscriptExonStartPosVec;
			vector<int> tmpTranscriptExonEndPosVec;
			transcriptInfo->copyExonPos2anotherVec(tmp,
				tmpTranscriptExonStartPosVec, tmpTranscriptExonEndPosVec);
			int tmpExonNum = transcriptInfo->returnTranscriptExonNum(tmp);
			//cout << "tmpExonNum: " << tmpExonNum << endl;
			for(int tmpExon = 0; tmpExon < tmpExonNum; tmpExon++)
			{	
				// generate SNP pos and base in each exon
				vector<int> candiSNPposVecInTmpExon;
				vector<string> candiSNPbaseVecInTmpExon;
				int tmpExonStartPos = tmpTranscriptExonStartPosVec[tmpExon];
				int tmpExonEndPos = tmpTranscriptExonEndPosVec[tmpExon];
				snpHashInfo.returnSNPvecWithinRegion(tmpTranscriptChromNameInt, 
					tmpExonStartPos, tmpExonEndPos, 
					candiSNPposVecInTmpExon, candiSNPbaseVecInTmpExon);
				tmpCandiSNPposVecVecInTranscript.push_back(candiSNPposVecInTmpExon);
				tmpCandiSNPbaseVecVecInTranscript.push_back(candiSNPbaseVecInTmpExon);
			}

			vector<int> tmpSNPposVec2add;
			vector<string> tmpSNPbaseVec2add;
			int tmpAccumulatedTranscriptSeqLen = 0;
			for(int tmpVecIndex = 0; tmpVecIndex < tmpCandiSNPposVecVecInTranscript.size(); tmpVecIndex ++)
			{
				// for each exon
				int tmpExonStartPos = tmpTranscriptExonStartPosVec[tmpVecIndex];
				int tmpExonEndPos = tmpTranscriptExonEndPosVec[tmpVecIndex];
				for(int tmpVecIndex_2 = 0; tmpVecIndex_2 < tmpCandiSNPposVecVecInTranscript[tmpVecIndex].size(); tmpVecIndex_2 ++)
				{
					int tmpSNPpos = (tmpCandiSNPposVecVecInTranscript[tmpVecIndex])[tmpVecIndex_2]; // genomic pos for SNP
					int tmpDistance2exonStart = tmpSNPpos - tmpExonStartPos + 1;
					int tmpSNPlocInTranscript = tmpAccumulatedTranscriptSeqLen + tmpDistance2exonStart; // transcript loc for SNP
					string tmpSNPbase = (tmpCandiSNPbaseVecVecInTranscript[tmpVecIndex])[tmpVecIndex_2]; // SNP base
					tmpSNPposVec2add.push_back(tmpSNPlocInTranscript);
					tmpSNPbaseVec2add.push_back(tmpSNPbase);
				}
				int tmpExonLength = tmpExonEndPos - tmpExonStartPos + 1;
				tmpAccumulatedTranscriptSeqLen += tmpExonLength;
			}
			SNPlocInTranscriptVecVec.push_back(tmpSNPposVec2add);
			SNPalterBaseInTranscriptVecVec.push_back(tmpSNPbaseVec2add);
		}
		//cout << "end of addSNP2transcript " << endl;
	}


	void outputSimulatedTransReadSeqAroundSNP(
		string& outputSimulatedReadFilePath, 
		Transcript_Set* transcriptInfo, Index_Info* indexInfo,
		int simulatedReadSeqLength, bool outputFastaReadWithOriOrSNPbaseBool)
	{
		//cout << "start to do outputSimulatedTransReadSeqAroundSNP ..." << endl;
		ofstream read_ofs(outputSimulatedReadFilePath.c_str());
		int transcriptNum = transcriptInfo->returnTranscriptNum();
		//cout << "transcriptNum: " << transcriptNum << endl;
		for(int tmp = 0; tmp < transcriptNum; tmp++)
		{
			//cout << "tmpIndex in transcriptome: " << tmp << endl;
			int tmpSNPnumInTranscript = SNPlocInTranscriptVecVec[tmp].size();
			//cout << "tmpSNPnumInTranscript: " << tmpSNPnumInTranscript << endl;
			if(tmpSNPnumInTranscript == 0)
				continue;
			string tmpTranscriptName = transcriptInfo->returnTranscriptID(tmp);
			string tmpTranscriptSeq = transcriptInfo->returnTranscriptSeq(tmp, indexInfo);
			string tmpUpdatedTranscriptSeq = tmpTranscriptSeq;
			//cout << "tmpTranscriptSeq: " << endl << tmpTranscriptSeq << endl;
			// update transcript seqeunces
			if(!outputFastaReadWithOriOrSNPbaseBool)
			{	
				for(int tmpSNPindex = 0; tmpSNPindex < tmpSNPnumInTranscript; tmpSNPindex ++)
				{
					int tmpSNPloc = (SNPlocInTranscriptVecVec[tmp])[tmpSNPindex];
					string tmpSNPbase = (SNPalterBaseInTranscriptVecVec[tmp])[tmpSNPindex];
					//cout << "tmpSNPbase: " << tmpSNPbase << endl;
					tmpUpdatedTranscriptSeq.replace(tmpSNPloc - 1, 1, tmpSNPbase);
				}
			}
			//cout << "tmpUpdatedTranscriptSeq: " << endl << tmpUpdatedTranscriptSeq << endl;

			string tmpSimulatedRead_aroundSNP_chrNameStr = transcriptInfo->returnTranscriptChromName(tmp);
			//cout << "tmpSimulatedRead_aroundSNP_chrNameStr: " << tmpSimulatedRead_aroundSNP_chrNameStr << endl;
			int tmpTranscript_seq_length = tmpUpdatedTranscriptSeq.length();
			//cout << "tmpTranscript_seq_length: " << tmpTranscript_seq_length << endl;
			// output simulated transReadSeq around SNPs
			for(int tmpSNPindex = 0; tmpSNPindex < tmpSNPnumInTranscript; tmpSNPindex ++)
			{
				//cout << "tmpSNPindex: " << tmpSNPindex << endl;
				int halfLength = (simulatedReadSeqLength - 1) / 2;

				int tmpSNPlocInTranscript = (SNPlocInTranscriptVecVec[tmp])[tmpSNPindex];
				//cout << "tmpSNPlocInTranscript: " << tmpSNPlocInTranscript << endl;
				int tmpSNPgenomicPos = transcriptInfo->returnGenomicPosWithLocInTranscript(tmp, tmpSNPlocInTranscript);
				//cout << "tmpSNPgenomicPos: " << tmpSNPgenomicPos << endl;
				int tmpSimulatedRead_seq_startLocInTranscript_initial = tmpSNPlocInTranscript - halfLength;
				int tmpSimulatedRead_seq_endLocInTranscript_initial = tmpSNPlocInTranscript + halfLength;
				//cout << "tmpSimulatedRead_seq_startLocInTranscript_initial: " << tmpSimulatedRead_seq_startLocInTranscript_initial << endl;
				//cout << "tmpSimulatedRead_seq_endLocInTranscript_initial: " << tmpSimulatedRead_seq_endLocInTranscript_initial << endl;
				int tmpSimulatedRead_seq_startLocInTranscript_final;

				int tmpSimulatedRead_aroundSNP_startPos;
				string tmpSimulatedRead_aroundSNP_cigarString;
				int tmpSimulatedRead_aroundSNP_SNPlocInReadSeq;

				if(simulatedReadSeqLength > tmpTranscript_seq_length)
					break;
				if((tmpSimulatedRead_seq_startLocInTranscript_initial >= 1)
					&&(tmpSimulatedRead_seq_endLocInTranscript_initial <= tmpTranscript_seq_length))
				{
					tmpSimulatedRead_seq_startLocInTranscript_final = tmpSimulatedRead_seq_startLocInTranscript_initial;
					tmpSimulatedRead_aroundSNP_startPos = transcriptInfo->returnGenomicPosWithLocInTranscript(tmp, tmpSimulatedRead_seq_startLocInTranscript_final);
					vector<Jump_Code> tmpJumpCodeVec;
					transcriptInfo->generateGenomicJumpCodeVecWithStartAndEndLocInTranscript(tmp, tmpSimulatedRead_seq_startLocInTranscript_final, 
						tmpSimulatedRead_seq_startLocInTranscript_final + 2 * halfLength, tmpJumpCodeVec);
					tmpSimulatedRead_aroundSNP_cigarString = this->jumpCodeVec2Str(tmpJumpCodeVec);
					tmpSimulatedRead_aroundSNP_SNPlocInReadSeq = halfLength + 1;
				}
				else if(tmpSimulatedRead_seq_startLocInTranscript_initial < 1)
				{
					tmpSimulatedRead_seq_startLocInTranscript_final = 1;
					tmpSimulatedRead_aroundSNP_startPos = transcriptInfo->returnGenomicPosWithLocInTranscript(tmp, 1);
					vector<Jump_Code> tmpJumpCodeVec;
					transcriptInfo->generateGenomicJumpCodeVecWithStartAndEndLocInTranscript(tmp, 1, simulatedReadSeqLength, tmpJumpCodeVec);
					tmpSimulatedRead_aroundSNP_cigarString = this->jumpCodeVec2Str(tmpJumpCodeVec);
					tmpSimulatedRead_aroundSNP_SNPlocInReadSeq = tmpSNPlocInTranscript;
				}
				else if(tmpSimulatedRead_seq_endLocInTranscript_initial > tmpTranscript_seq_length)
				{
					tmpSimulatedRead_seq_startLocInTranscript_final = tmpTranscript_seq_length - simulatedReadSeqLength + 1;
					tmpSimulatedRead_aroundSNP_startPos = transcriptInfo->returnGenomicPosWithLocInTranscript(tmp, tmpSimulatedRead_seq_startLocInTranscript_final);
					vector<Jump_Code> tmpJumpCodeVec;
					transcriptInfo->generateGenomicJumpCodeVecWithStartAndEndLocInTranscript(tmp, tmpSimulatedRead_seq_startLocInTranscript_final, 
						tmpTranscript_seq_length, tmpJumpCodeVec);					
					tmpSimulatedRead_aroundSNP_cigarString = this->jumpCodeVec2Str(tmpJumpCodeVec);
					tmpSimulatedRead_aroundSNP_SNPlocInReadSeq = simulatedReadSeqLength - (tmpTranscript_seq_length - tmpSNPlocInTranscript);
				}
				else
				{}
				string tmpReadLengthStr = int_to_str(simulatedReadSeqLength);
				string tmpSimulatedRead_aroundSNP_startPosStr = int_to_str(tmpSimulatedRead_aroundSNP_startPos);
				string tmpSimulatedRead_aroundSNP_SNPlocInReadSeqStr = int_to_str(tmpSimulatedRead_aroundSNP_SNPlocInReadSeq);
				string tmpSimulatedRead_aroundSNP_genomicPosStr = int_to_str(tmpSNPgenomicPos);
				string tmpSimulatedRead_aroundSNP_SNPlocInTranscriptStr = int_to_str(tmpSNPlocInTranscript);
				// 1.SAM_chrName 2.SAM_startPos 3.SAM_cigarString 4.SNPlocInRead 5.readLength 6. SNPposInChr 7.transcriptName 8.SNPlocInTranscript
				string tmpSimulatedRead_name = ">" + tmpSimulatedRead_aroundSNP_chrNameStr + ":"
					+ tmpSimulatedRead_aroundSNP_startPosStr + ":" + tmpSimulatedRead_aroundSNP_cigarString + ":"
					+ tmpSimulatedRead_aroundSNP_SNPlocInReadSeqStr + ":" + tmpReadLengthStr + ":" 
					+ tmpSimulatedRead_aroundSNP_genomicPosStr + ":" + tmpTranscriptName + ":"
					+ tmpSimulatedRead_aroundSNP_SNPlocInTranscriptStr;

				string tmpSimulatedRead_seq = tmpUpdatedTranscriptSeq.substr(
					tmpSimulatedRead_seq_startLocInTranscript_final - 1, simulatedReadSeqLength);

				read_ofs << tmpSimulatedRead_name << endl 
					<< tmpSimulatedRead_seq << endl; 
			}
		}
		read_ofs.close();
	}
};

#endif