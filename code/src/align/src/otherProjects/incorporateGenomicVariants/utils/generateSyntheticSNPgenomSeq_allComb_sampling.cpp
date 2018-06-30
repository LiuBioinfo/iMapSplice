// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "../general/SNPhash_info.h"

using namespace std;

//typedef set<int> 

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputIndexFolderPath" << endl;
		cout << "#2 inputAnnotationPath" << endl;
		cout << "#3 inputFormattedSNPfile" << endl;
		cout << "#4 outputFilePrefix" << endl;
		cout << "#5 syntheticSeqLength" << endl;
		cout << "#6 samplingNumMax" << endl;
		exit(1);
	}
	int totalSNPnum = 0;
	int aggregativeSNPnum = 0;

	string inputFormattedSNPfile = argv[3];
	string syntheticSeqLengthStr = argv[5];
	string samplingNumMaxStr = argv[6];
	int syntheticSeqLength = atoi(syntheticSeqLengthStr.c_str());
	int halfLength = (syntheticSeqLength - 1) / 2;
	int samplingNumMax = atoi(samplingNumMaxStr.c_str());
	cout << "syntheticSeqLength: " << syntheticSeqLength << endl;
	cout << "halfLength: " << halfLength << endl;  
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_prefix = outputFolderStr;
	string outputFile_SNPlist = output_prefix + "/SNPinAnn.txt";
	string outputFile_syntheticSeqFa = output_prefix + "/SNPinAnn.fa";
	string outputFile_log = output_prefix + "/log";
	ofstream log_ofs(outputFile_log.c_str());

	log_ofs << "Command: \n" << argv[0] << endl << argv[1] << endl << argv[2] 
		<< endl << argv[3] << endl << argv[4] << endl << argv[5] << endl;

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	////////////////// update chr seq with SNPs << endl;
	cout << "start to update chr seq with SNPs" << endl;
	indexInfo->insertSNP2chromStr(inputFormattedSNPfile, log_ofs);
	////////////////// transcript info loading .... //////////
	string transcript_file_path = argv[2];
	cout << "start to load transcriptInfo " << endl;
	Transcript_Set* transcriptInfo = new Transcript_Set();
	ifstream transcript_file_ifs(transcript_file_path.c_str());
	string transcript_type_GAF = "GAF";
	transcriptInfo->extractTranscript(transcript_file_ifs, indexInfo, transcript_type_GAF);

	///////////////////////  load SNPs //////////////////////////////////
	SNPhash_Info tmpSNPhashInfo;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	tmpSNPhashInfo.generateSNPhash_formattedSNPfile(inputFormattedSNPfile, indexInfo);

	/////////////////////// start to extract SNPs in transcriptome ///////////////////////////
	vector< set<int> > SNPposSetVec;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpSNPposSet;
		SNPposSetVec.push_back(tmpSNPposSet);
 	}

 	cout << "start to extract SNP from gene ann" << endl;
	int transcriptNum = transcriptInfo->returnTranscriptNum();
	cout << "transcriptNum: " << transcriptNum << endl;
	for(int tmp = 0; tmp < transcriptNum; tmp++)
	{	
		string tmpTranscriptChromName = transcriptInfo->returnTranscriptChromName(tmp);
		int tmpTranscriptChromNameInt = indexInfo->convertStringToInt(tmpTranscriptChromName);
		vector<int> tmpTranscriptExonStartPosVec;
		vector<int> tmpTranscriptExonEndPosVec;
		transcriptInfo->copyExonPos2anotherVec(tmp,	tmpTranscriptExonStartPosVec, tmpTranscriptExonEndPosVec);
		int tmpExonNum = transcriptInfo->returnTranscriptExonNum(tmp);
		for(int tmpExon = 0; tmpExon < tmpExonNum; tmpExon++)
		{
			// generate SNP pos and base in each exon
			vector<int> candiSNPposVecInTmpExon;
			vector<string> candiSNPbaseVecInTmpExon;
			int tmpExonStartPos = tmpTranscriptExonStartPosVec[tmpExon];
			int tmpExonEndPos = tmpTranscriptExonEndPosVec[tmpExon];
			tmpSNPhashInfo.returnSNPvecWithinRegion(tmpTranscriptChromNameInt, tmpExonStartPos, 
				tmpExonEndPos, candiSNPposVecInTmpExon, candiSNPbaseVecInTmpExon);
			for(int tmpSNP = 0; tmpSNP < candiSNPposVecInTmpExon.size(); tmpSNP++)
			{
				int tmpSNPpos = candiSNPposVecInTmpExon[tmpSNP];
				SNPposSetVec[tmpTranscriptChromNameInt].insert(tmpSNPpos);
			}
		}	
	}

	cout << "start to output SNP and synthetic seq from SNP in gene ann" << endl;
	ofstream SNPinAnn_list_ofs(outputFile_SNPlist.c_str());
	ofstream SNPinAnn_fa_ofs(outputFile_syntheticSeqFa.c_str());
	for(int tmpChr = 0; tmpChr < SNPposSetVec.size(); tmpChr++)
	{
		int tmpChrLength = indexInfo->returnChromLength(tmpChr);
		int tmpStartPosInChr = 1;
		int tmpEndPosInChr = tmpChrLength;
		if(syntheticSeqLength > tmpChrLength)
		{
			cout << "tmpChrLength: " << tmpChrLength << endl;
			cout << "syntheticSeqLength: " << syntheticSeqLength << endl;
			continue;
		}
		for(set<int>::iterator setIntIter = SNPposSetVec[tmpChr].begin();
			setIntIter != SNPposSetVec[tmpChr].end(); setIntIter ++)
		{
			totalSNPnum ++;
			int tmpSNPpos = (*setIntIter);
			string tmpSNPbase;
			bool searchAndReturnSNPinHash_bool = tmpSNPhashInfo.searchAndReturnSNPbase(tmpChr, tmpSNPpos, tmpSNPbase);
			if(!searchAndReturnSNPinHash_bool)
			{
				cout << "SNP inconsitent, exiting ......" << endl;
				exit(1);
			}
			else
			{
				string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
				string tmpRefBase = indexInfo->returnChromStrSubstr(tmpChr, tmpSNPpos, 1);
				SNPinAnn_list_ofs << tmpChrName << "\t" << tmpSNPpos << "\t" << tmpRefBase << "\t" << tmpSNPbase << endl;
				int tmpSyntheticSeq_startPosInChr_ini = tmpSNPpos - halfLength;
				int tmpSyntheticSeq_endPosInChr_ini = tmpSNPpos + halfLength;				
				//int tmpSyntheticSeq_endPosInChr_final;
				int tmpSyntheticSeq_startPosInChr_final;
				int tmpSNPlocInSyntheticSeq_final;
				if((tmpSyntheticSeq_startPosInChr_ini >= 1)&&(tmpSyntheticSeq_endPosInChr_ini <= tmpChrLength))
				{
					tmpSyntheticSeq_startPosInChr_final = tmpSyntheticSeq_startPosInChr_ini;
					tmpSNPlocInSyntheticSeq_final = halfLength + 1;
				}
				else if(tmpSyntheticSeq_startPosInChr_ini < 1)
				{
					tmpSyntheticSeq_startPosInChr_final = 1;
					tmpSNPlocInSyntheticSeq_final = tmpSNPpos;
				}
				else if(tmpSyntheticSeq_endPosInChr_ini > tmpChrLength)
				{
					tmpSyntheticSeq_startPosInChr_final = tmpChrLength - 2 * halfLength;
					tmpSNPlocInSyntheticSeq_final = tmpChrLength - halfLength;
				}
				else
				{}
				int tmpSyntheticSeq_endPosInChr_final = tmpSyntheticSeq_startPosInChr_final + syntheticSeqLength - 1;
				vector<int> tmpCandiSNPposVecInTmpSyntheticSeq_raw;
				vector<string> tmpCandiSNPbaseVecInTmpSyntheticSeq_raw;
				tmpSNPhashInfo.returnSNPvecWithinRegion(tmpChr, 
					tmpSyntheticSeq_startPosInChr_final, tmpSyntheticSeq_endPosInChr_final, 
					tmpCandiSNPposVecInTmpSyntheticSeq_raw, tmpCandiSNPbaseVecInTmpSyntheticSeq_raw);
				vector<int> tmpCandiSNPposVecInTmpSyntheticSeq_other;
				vector<string> tmpCandiSNPbaseVecInTmpSyntheticSeq_other;				
				for(int tmp = 0; tmp < tmpCandiSNPposVecInTmpSyntheticSeq_raw.size(); tmp++)
				{
					if(tmpCandiSNPposVecInTmpSyntheticSeq_raw[tmp] != tmpSNPpos)
					{
						tmpCandiSNPposVecInTmpSyntheticSeq_other.push_back(tmpCandiSNPposVecInTmpSyntheticSeq_raw[tmp]);
						tmpCandiSNPbaseVecInTmpSyntheticSeq_other.push_back(tmpCandiSNPbaseVecInTmpSyntheticSeq_raw[tmp]);
					}
				}
				int tmpSyntheticSeq_SNPnum_other = tmpCandiSNPposVecInTmpSyntheticSeq_other.size();
				if(tmpSyntheticSeq_SNPnum_other > 0)
				{
					vector< vector<bool> > tmpSNPexpressOrNotVecVec;
					int tmp_totalCasesNum = pow(2, tmpSyntheticSeq_SNPnum_other);
					vector< int > tmpAndBaseVec;
					for(int tmp = 0; tmp < tmpSyntheticSeq_SNPnum_other; tmp++)
					{
						int tmpAndBase = pow(2, tmp);
						tmpAndBaseVec.push_back(tmpAndBase);
					}
					if(tmp_totalCasesNum > samplingNumMax)
					{
						tmp_totalCasesNum = samplingNumMax;
						aggregativeSNPnum ++;
					}
					for(int tmpCase = 0; tmpCase < tmp_totalCasesNum; tmpCase++)
					{
						vector<bool> tmpSNPexpressOrNotVec;
						for(int tmpAndBaseIndex = 0; tmpAndBaseIndex < tmpAndBaseVec.size(); tmpAndBaseIndex ++)
						{
							bool tmpSNP_expressOrNot_bool = (tmpCase & tmpAndBaseVec[tmpAndBaseIndex]);
							tmpSNPexpressOrNotVec.push_back(tmpSNP_expressOrNot_bool);
						}
						tmpSNPexpressOrNotVecVec.push_back(tmpSNPexpressOrNotVec);
					}
					string tmpSyntheticGenomeSeq_raw = indexInfo->returnChromStrSubstr(tmpChr, 
						tmpSyntheticSeq_startPosInChr_final, syntheticSeqLength);
					for(int tmpCase = 0; tmpCase < tmp_totalCasesNum; tmpCase++)
					{
						SNPinAnn_fa_ofs << ">" << tmpChrName << ":" << tmpSyntheticSeq_startPosInChr_final << ":" << syntheticSeqLength
							<< "M:" << tmpSNPlocInSyntheticSeq_final << ":" << syntheticSeqLength << ":" << tmpCase + 1 << endl; 
						string tmpSyntheticGenomeSeq_tmpComb = tmpSyntheticGenomeSeq_raw;
						for(int tmpSNPindex = 0; tmpSNPindex < tmpSyntheticSeq_SNPnum_other; tmpSNPindex ++)
						{
							if((tmpSNPexpressOrNotVecVec[tmpCase])[tmpSNPindex])
								tmpSyntheticGenomeSeq_tmpComb.replace(
									tmpCandiSNPposVecInTmpSyntheticSeq_other[tmpSNPindex] - tmpSyntheticSeq_startPosInChr_final, 
									1, tmpCandiSNPbaseVecInTmpSyntheticSeq_other[tmpSNPindex]);
						}
						SNPinAnn_fa_ofs << tmpSyntheticGenomeSeq_tmpComb << endl;
					}
				}
				else
				{
					SNPinAnn_fa_ofs << ">" << tmpChrName << ":" << tmpSyntheticSeq_startPosInChr_final << ":" << syntheticSeqLength 
						<< "M:" << tmpSNPlocInSyntheticSeq_final << ":" << syntheticSeqLength << ":1" << endl;
					SNPinAnn_fa_ofs << indexInfo->returnChromStrSubstr(tmpChr, tmpSyntheticSeq_startPosInChr_final, syntheticSeqLength) << endl;		
				}
			}			
		}
	}
	SNPinAnn_list_ofs.close();
	SNPinAnn_fa_ofs.close();

	log_ofs << "Total SNP #: " << totalSNPnum << endl;
	log_ofs << "Aggregative SNP #: " << aggregativeSNPnum << endl;

	cout << "All jobs done!" << endl;
	log_ofs << "All jobs done!" << endl;
	log_ofs.close();
	transcriptInfo->memoryFree();
	delete transcriptInfo;
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}