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
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

#include "stats_info.h"
#include "constantDefinitions.h"
//#include "general/option_info.h"
#include "general/read_block_test.h"
#include "general/bwtmap_info.h"
#include "general/DoubleAnchorScore.h"
#include "general/sbndm.h"
#include "general/otherFunc.h"
#include "general/index_info.h"
#include "general/enhanced_suffix_array_info.h"
#include "general/annotation_info.h"
#include "phase1/repeatRegion.h"
#include "general/segmentMapping.h"
#include "general/splice_info.h"
#include "general/fixGapRelationParameters.h"
#include "general/read_info.h"
#include "general/seg_info.h"
#include "general/fixDoubleAnchorMatch_info.h"
#include "general/fixDoubleAnchorInsertion_info.h"
#include "general/fixDoubleAnchorDeletion_info.h"
#include "general/fixDoubleAnchorSplice_complicate_info.h"
#include "general/fixDoubleAnchorSplice_info.h"
#include "general/path_info.h"
#include "general/gap_info.h"
#include "general/align_info.h"
#include "general/peAlign_info.h"

#include "phase2/spliceJunctionHash_info.h"
#include "phase2/unmapEnd_info.h"
#include "phase2/unfixedHead.h"
#include "phase2/unfixedTail.h"
#include "phase2/incompleteLongHead.h"
#include "phase2/incompleteLongTail.h"
#include "phase2/sam2junc.h"
#include "fixHeadTail.h"
#include "phase2/fixOneEndUnmapped.h"
#include "fixPhase1.h"
#include "general/readSeqPreProcessing.h"
#include "general/headerSection_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc != 7)
	{
		cout << "Exetuable GlobalIndex LocalIndex inputAlignmentInfo_unpaired outputFolder thread_num fasta_or_fastq" << endl;
		exit(1);
	}

	string globalIndexStr = argv[1];
	string localIndexStr = argv[2];
	string inputAlignInfoStr_unpaired = argv[3];
	string outputFolderStr = argv[4];
	string threadNumStr = argv[5];
	string fasta_or_fastq_str = argv[6];

	int normalRecordNum_fixOneEndUnmapped = 100000;
	bool load2ndLevelIndexBool = true;
	bool Do_extendHeadTail_fixOneEndUnmapped = true;
	bool annotation_provided_bool = false;
	bool Do_annotation_only_bool = false;

	int threads_num = atoi(threadNumStr.c_str());
	bool fasta_or_fastq_bool;
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta")||(fasta_or_fastq_str == "FASTA"))
	{
		fasta_or_fastq_bool = true;
	}
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq")||(fasta_or_fastq_str == "FASTQ"))
	{
		fasta_or_fastq_bool = false;
	}
	else
	{
		cout << "input file format invalid" << endl;
		exit(1);
	}

	Annotation_Info* annotationInfo = new Annotation_Info();

	string mkdirOutputCommand_phase2 = "mkdir -p " + outputFolderStr;
	system(mkdirOutputCommand_phase2.c_str());

	string outputDirStr = outputFolderStr + "/";
	string log_ofs_str = outputDirStr + "/process.log";
	ofstream log_ofs(log_ofs_str.c_str());
	string runtime_ofs_str = outputDirStr + "/runtime.log";
	ofstream runtime_log_ofs(runtime_ofs_str.c_str());
	string settings_ofs_str = outputDirStr + "/settings.log";
	ofstream settings_log_ofs(settings_ofs_str.c_str());

	string indexStr = globalIndexStr + "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, settings_log_ofs);

	settings_log_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	settings_log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	settings_log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	log_ofs << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	time_t nowtime;
	struct tm *local;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load 2nd level index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index starts ......" << endl << endl; 

	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;

	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;

	string secondLevelIndexStr = localIndexStr + "/";

	if(load2ndLevelIndexBool)
	{
		log_ofs << "start to load second-level index ..." << endl;
		
		int secondLevelIndexNO = 0;
		for(int tmpChrNO = 0; tmpChrNO < indexInfo->returnChromNum(); tmpChrNO ++)
		{
			for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->returnSecondLevelIndexPartsNum(tmpChrNO)); tmpSecondLevelIndexNO ++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
				string tmpFileNumStr = tmpFileNumChar;
				
				string inputIndexFileStr = secondLevelIndexStr + "/" + 
					//indexInfo->chrNameStr[tmpChrNO] 
					indexInfo->returnChrNameStr(tmpChrNO)
					+ "/" 
					//+ indexInfo->chrNameStr[tmpChrNO] + "_part."
					+ tmpFileNumStr + "/";//"." + "test3_";


				string secondLevelIndexFileChromStr = inputIndexFileStr + "chrom"; 
				ifstream secondLevelChrom_file_ifs(secondLevelIndexFileChromStr.c_str(), ios::binary);
				string secondLevelIndexFileSaStr = inputIndexFileStr + "SA";
				ifstream secondLevelSA_file_ifs(secondLevelIndexFileSaStr.c_str(), ios::binary);

					string secondLevelIndexFileLcpCompressStr = inputIndexFileStr + "_lcpCompress";	
					ifstream secondLevelLcpCompress_file_ifs(secondLevelIndexFileLcpCompressStr.c_str(), ios::binary);	
					string secondLevelIndexFileChildTabStr = inputIndexFileStr + "childTab";	
					ifstream secondLevelChildTab_file_ifs(secondLevelIndexFileChildTabStr.c_str(), ios::binary);
					string secondLevelIndexFileDetChildStr = inputIndexFileStr + "detChild";	
					ifstream secondLevelDetChild_file_ifs(secondLevelIndexFileDetChildStr.c_str(), ios::binary);					

				int sizeOfIndex = indexInfo->returnSecondLevelIndexNormalSize() + 1;
				char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					tmpSecondLevelChrom[tmpMallocSpace] = '0';
				}
				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
				{
					//(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
				}

				bool No_ATGC_Bool = true;
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					char ch = tmpSecondLevelChrom[tmpMallocSpace];
					if((ch == 'A')||(ch == 'T')||(ch == 'G')||(ch == 'C'))
					{
						No_ATGC_Bool = false;
						break;
					}
				}				
				if(No_ATGC_Bool)
				{
					//(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
				}	

				secondLevelChrom.push_back(tmpSecondLevelChrom);
				
				unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));
				secondLevelSa.push_back(tmpSecondLevelSa);

				BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress.push_back(tmpSecondLevelLcpCompress);
					
				unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab.push_back(tmpSecondLevelChildTab);

				BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild.push_back(tmpSecondLevelDetChild);

				secondLevelChrom_file_ifs.close();
				secondLevelSA_file_ifs.close();

				secondLevelLcpCompress_file_ifs.close();
				secondLevelChildTab_file_ifs.close();
				secondLevelDetChild_file_ifs.close();							

				secondLevelIndexNO ++;
			}
			log_ofs << "finish loading 2nd-level index of " << indexInfo->returnChrNameStr(tmpChrNO) << endl; 
		}
		log_ofs << "finish loading ALL 2nd-level index !" << endl;
		log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
	}


	string OutputSamFile_oneEndMapped = outputDirStr + "oneEndUnmapped.pairedComplete.sam";
	ofstream OutputSamFile_oneEndMapped_ofs(OutputSamFile_oneEndMapped.c_str());

	string OutputAlignInfoFile_oneEndMapped_incomplete = outputDirStr + "oneEndUnmapped.pairedInomplete.alignInfo";
	ofstream OutputAlignInfoFile_oneEndMapped_incomplete_ofs(OutputAlignInfoFile_oneEndMapped_incomplete.c_str());

	string OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore = outputDirStr + "oneEndUnmapped.bothEndsUnmapped_lowScore.sam";
	ofstream OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());

	string OutputSamFile_oneEndMapped_unpair = outputDirStr + "oneEndUnmapped.unpaired.sam";
	ofstream OutputSamFile_oneEndMapped_unpair_ofs(OutputSamFile_oneEndMapped_unpair.c_str());

	string OutputSamFile_oneEndMapped_alignInfo = outputDirStr + "oneEndUnmapped.pairedComplete.sam_alignInfo";
	ofstream OutputSamFile_oneEndMapped_alignInfo_ofs(OutputSamFile_oneEndMapped_alignInfo.c_str());

	int tmpRecordNum_oneEndUnmapped = 0;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	log_ofs << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;
	runtime_log_ofs << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;
	cout << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;
	string oneEndMappedFileStr = inputAlignInfoStr_unpaired;
	ifstream inputRecord_ifs(oneEndMappedFileStr.c_str());
	string line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11;
		
	int normalRecordNum = normalRecordNum_fixOneEndUnmapped; //1000000;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;// = normalRecordNum;	


		for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
		{
			vector<string> line1StrVec(normalRecordNum);
			vector<string> line2StrVec(normalRecordNum);
			vector<string> line3StrVec(normalRecordNum);
			vector<string> line4StrVec(normalRecordNum);
			vector<string> line5StrVec(normalRecordNum);
			vector<string> line6StrVec(normalRecordNum);
			vector<string> line7StrVec(normalRecordNum);
			vector<string> line8StrVec(normalRecordNum);
			vector<string> line9StrVec(normalRecordNum);
			vector<string> line10StrVec(normalRecordNum);

			vector<string> peAlignInfoVec_fixUnpair(normalRecordNum);
			vector<string> peAlignSamVec_fixUnpair(normalRecordNum);
			vector<string> peAlignSamVec_unpair_fixUnpair(normalRecordNum);
			vector<string> peAlignInfoVec_pair_complete_fixUnpair(normalRecordNum);
			vector<string> peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair(normalRecordNum);

			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << endl << "[" << asctime(local) 
				<< "start to input oneEndUnmapped records, turn: " << tmpTurn+1 << endl;
			cout << endl << endl << "[" << asctime(local) 
				<< "start to input oneEndUnmapped records, turn: " << tmpTurn+1 << endl;

			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(inputRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}				
				getline(inputRecord_ifs, line11);
				if(inputRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}					
				getline(inputRecord_ifs, line1);
				getline(inputRecord_ifs, line2);
				getline(inputRecord_ifs, line3);
				getline(inputRecord_ifs, line4);
				getline(inputRecord_ifs, line5);
				getline(inputRecord_ifs, line6);
				getline(inputRecord_ifs, line7);
				getline(inputRecord_ifs, line8);
				getline(inputRecord_ifs, line9);
				getline(inputRecord_ifs, line10);
				//getline(inputRecord_ifs, line11);

				line1StrVec[recordNumTmp] = line1;
				line2StrVec[recordNumTmp] = line2;
				line3StrVec[recordNumTmp] = line3;
				line4StrVec[recordNumTmp] = line4;
				line5StrVec[recordNumTmp] = line5;
				line6StrVec[recordNumTmp] = line6;
				line7StrVec[recordNumTmp] = line7;
				line8StrVec[recordNumTmp] = line8;
				line9StrVec[recordNumTmp] = line9;
				line10StrVec[recordNumTmp] = line10;

			}	
		
			runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local) 
				<< "finish reading record, turn: " << tmpTurn+1 << endl;
			runtime_log_ofs << endl << "[" << asctime(local) 
				<< "start to fix oneEndUnmapped, turn: " << tmpTurn+1 << endl;	
			cout << endl << "[" << asctime(local) 
				<< "finish reading record, turn: " << tmpTurn+1 << endl;
			cout << endl << "[" << asctime(local) 
				<< "start to fix oneEndUnmapped, turn: " << tmpTurn+1 << endl;	
			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				//tmpRecordNum_oneEndUnmapped ++;
				int threadNO = omp_get_thread_num();

				PE_Read_Info peReadInfo;// = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_getline(
					line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], line3StrVec[tmpOpenMP],
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP], line6StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo, indexInfo, fasta_or_fastq_bool);		

				//cout << "start to output readInfo: " << endl;
				//peReadInfo.printPEreadInfo();
				//cout << "finish output readInfo: " << endl;

				FixOneEndUnmappedInfo* fixOneEndUnmappedInfo = new FixOneEndUnmappedInfo();
				fixOneEndUnmappedInfo->fixOneEndUnmapped(peReadInfo, peAlignInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo);

				peAlignInfo->alignmentFilter_fixOneEndUnmapped_SJpenalty(peReadInfo.returnReadSeqLength_1(),
					peReadInfo.returnReadSeqLength_2());

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				bool allUnpairedAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();

				bool unique_bool = peAlignInfo->checkUniqueOrMulti();

				//statsInfo->increNum_fixUnpaired(threadNO, pairExistsBool, allAlignmentCompleteBool, 
				//	allUnpairedAlignmentCompleteBool, unique_bool);

				//string tmpPeAlignSamStr, tmpPeAlignInfoStr, tmpPeAlignSamStr_unpair, tmpPeAlignSamStr_bothEndsUnmapped_lowScore;

				peAlignSamVec_fixUnpair[tmpOpenMP] = "";
				peAlignInfoVec_fixUnpair[tmpOpenMP] = "";
				peAlignSamVec_unpair_fixUnpair[tmpOpenMP] = "";
				peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmpOpenMP] = "";


				if(pairExistsBool && allAlignmentCompleteBool) // some pair exists, all completed, print out paired SAM info
				{
					bool completeAlign_lowScore_bool = 
						peAlignInfo->alignPairScoreTooLow_bool(ALIGNMENT_SCORE_MIN_OUTPUT);
					if(completeAlign_lowScore_bool)
					{
						//statsInfo->increLowScoreComplete_fixUnpair(threadNO, unique_bool);
						peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmpOpenMP] = peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool);					
					}
					else
					{	
						peAlignSamVec_fixUnpair[tmpOpenMP] = peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo,
							fasta_or_fastq_bool);
					}
				}
				else if(pairExistsBool && (!allAlignmentCompleteBool)) // pair exists, incomplete
				{
					peAlignInfoVec_fixUnpair[tmpOpenMP] = peAlignInfo->getTmpAlignInfoForFinalPair(
						 		peReadInfo.returnReadName_1(), peReadInfo.returnReadName_2(), 
								peReadInfo.returnReadSeq_1(), peReadInfo.returnReadSeq_2(),
								peReadInfo.returnReadQual_1(), peReadInfo.returnReadQual_2(),
								fasta_or_fastq_bool);
				}
				else
				{
					peAlignSamVec_unpair_fixUnpair[tmpOpenMP] = peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
						peReadInfo, fasta_or_fastq_bool);				
				}

				// peAlignSamVec_fixUnpair[tmpOpenMP] = tmpPeAlignSamStr;
				// peAlignInfoVec_fixUnpair[tmpOpenMP] = tmpPeAlignInfoStr;
				// peAlignSamVec_unpair_fixUnpair[tmpOpenMP] = tmpPeAlignSamStr_unpair;
				// peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmpOpenMP] = tmpPeAlignSamStr_bothEndsUnmapped_lowScore;
				
				peAlignInfo->memoryFree();
				delete peAlignInfo;

			}

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish fixing oneEndUnmapped, turn: " << tmpTurn+1 << endl;// << endl;
			runtime_log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
			cout << endl << "[" << asctime(local)
				<< "finish fixing oneEndUnmapped, turn: " << tmpTurn+1 << endl;// << endl;
			cout << "start to output ... turn: " << tmpTurn+1 << endl;
						
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				if(peAlignSamVec_fixUnpair[tmp] != "")
				{
					OutputSamFile_oneEndMapped_ofs << peAlignSamVec_fixUnpair[tmp] << endl;
					//PairedReadNum ++;
					string().swap(peAlignSamVec_fixUnpair[tmp]);
				}
				if(peAlignInfoVec_fixUnpair[tmp] != "")
				{
					OutputAlignInfoFile_oneEndMapped_incomplete_ofs << peAlignInfoVec_fixUnpair[tmp] << endl;// << endl;
					string().swap(peAlignInfoVec_fixUnpair[tmp]);
				}
				if(peAlignSamVec_unpair_fixUnpair[tmp] != "")
				{
					OutputSamFile_oneEndMapped_unpair_ofs << peAlignSamVec_unpair_fixUnpair[tmp] << endl;
					//UnpairedReadNum ++;
					string().swap(peAlignSamVec_unpair_fixUnpair[tmp]);
				}
				if(peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] != "")
				{
					OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs << peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] << endl;
					string().swap(peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp]);
				}
			}		
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish output, turn: " << tmpTurn+1 << endl << endl;// << endl;
			cout << endl << "[" << asctime(local)
				<< "finish output, turn: " << tmpTurn+1 << endl << endl;// << endl;
		}		
		inputRecord_ifs.close();

	delete annotationInfo;

	return 0;
}