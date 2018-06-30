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

#include "../../general/extractUnmapAlignment2ReadFile.h"
#include "../../phase1/arrayQueue.h"
#include "../../stats_info.h"
#include "../../constantDefinitions.h"
#include "../../general/option_info.h"
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/otherFunc.h"
#include "../../general/index_info.h"
#include "../../general/enhanced_suffix_array_info.h"
#include "../../general/annotation_info.h"
#include "../../phase1/repeatRegion.h"
#include "../../general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "../../general/splice_info.h"
#include "../../general/fixGapRelationParameters.h"
#include "../../general/read_info.h"
#include "../../general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "../../general/fixDoubleAnchorNWDP_info.h"
#include "../../general/fixDoubleAnchorMatch_info.h"
#include "../../general/fixDoubleAnchorInsertion_info.h"
#include "../../general/fixDoubleAnchorDeletion_info.h"
#include "../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../general/fixDoubleAnchorSplice_info.h"
#include "../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../general/path_info.h"
#include "../../general/gap_info.h"
#include "../../general/align_info.h"
#include "../../general/peAlign_info.h"
#include "../../general/groupSeg_info.h"
#include "../../general/alignInferJunctionHash_info_vec.h"
#include "../../phase2/spliceJunctionHash_info.h"
#include "../../phase2/unmapEnd_info.h"
#include "../../phase2/unfixedHead.h"
#include "../../phase2/unfixedTail.h"
#include "../../phase2/incompleteLongHead.h"
#include "../../phase2/incompleteLongTail.h"
#include "../../phase2/sam2junc.h"
#include "../../fixHeadTail.h"
#include "../../phase2/fixOneEndUnmapped.h"
#include "../../fixPhase1.h"
#include "../../general/readSeqPreProcessing.h"
#include "../../general/headerSection_info.h"
#include "../../general/otherFunc2.h"
#include "general/segMap2aluFamily_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputSingleAluIndexPath inputRead_1 inputRead_2 outputFolder threads_num" << endl;
		exit(1);
	}
	string threads_num_str = argv[5];
	int threads_num = atoi(threads_num_str.c_str());

    string outputDirStr = argv[4];
   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());	
   	string statsStr = outputDirStr + "/stats.txt";
   	ofstream stats_ofs(statsStr.c_str());
   	string progressLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(progressLogStr.c_str());

	string indexStr = argv[1];
	string preIndexArrayPreStr = indexStr;
    preIndexArrayPreStr.append("/");
    indexStr.append("/");

	string preIndexMapLengthArrayStr = preIndexArrayPreStr; preIndexMapLengthArrayStr.append("_MapLength"); ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; preIndexIntervalStartArrayStr.append("_IntervalStart"); ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; preIndexIntervalEndArrayStr.append("_IntervalEnd"); ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
 	log_ofs << "finish loading preIndex ..." << endl;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	log_ofs << "index: " << indexStr << endl;
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	log_ofs << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;
	string SA_file = indexStr; SA_file.append("_SA"); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); 
	string childTab_file = indexStr; childTab_file.append("_childTab"); 
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); 	
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 	
	ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);

	cout << endl << "... start to load enhanced Suffix Array ......" << endl << endl; 
	log_ofs << endl << "... start to load enhanced Suffix Array ......" << endl << endl;
	log_ofs << "start to load SA" << endl;
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load lcpCompress" << endl;
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "start to load childTab " << endl;
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load detChild" << endl;
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "All index files loaded" << endl;
	cout << "... all index loaded ......" << endl << endl; 
	log_ofs << "... all index loaded ......" << endl << endl;

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;
	int normalRecordNum = 400000;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;
	int readPairNum = 0;
	bool fasta_or_fastq_bool = false;
	bool SE_or_PE_bool = false;
	int readTotalNum = 0;

	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readQualSeq1Vec(normalRecordNum);
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readQualSeq2Vec(normalRecordNum);

	vector<int> readSegmentMapScoreVec(normalRecordNum);
	vector<string> readSegMapInfoStrVec(normalRecordNum);
	cout << "... segment mapping process starts ......" << endl << endl; 
	log_ofs << "... segment mapping process starts ......" << endl << endl; 

	string inputPeReadPath_1 = argv[2];
	string inputPeReadPath_2 = argv[3];
	string outputSegmentMapScore = outputDirStr + "/score.txt";
	string outputSegLen = outputDirStr + "/segLen.txt";
	ifstream inputRead_ifs(inputPeReadPath_1.c_str());
	ifstream inputRead_PE_ifs(inputPeReadPath_2.c_str());
	ofstream score_ofs(outputSegmentMapScore.c_str());
	ofstream segLen_ofs(outputSegLen.c_str());
	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();
	for(tmpTurn = 0; ; tmpTurn++)
	{
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		log_ofs << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		realRecordNum = normalRecordNum;
		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		getline(inputRead_ifs, line1); // readName_1
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		getline(inputRead_ifs, line2); // readSeq_1
    		getline(inputRead_ifs, line3);
    		getline(inputRead_ifs, line4);
	    	getline(inputRead_PE_ifs, line1_PE); // readName_2
	    	getline(inputRead_PE_ifs, line2_PE); // readSeq_2
			getline(inputRead_PE_ifs, line3_PE);
	  		getline(inputRead_PE_ifs, line4_PE);	    	
    		readName1Vec[recordNumTmp] = line1.substr(1);
    		line2_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2);		
    		readSeq1Vec[recordNumTmp] = line2_afterProcess;    		
    		readQualSeq1Vec[recordNumTmp] = line4;    				
	    	readName2Vec[recordNumTmp] = line1_PE.substr(1);
	    	line2_PE_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2_PE);
	    	readSeq2Vec[recordNumTmp] = line2_PE_afterProcess;
	   		readQualSeq2Vec[recordNumTmp] = line4_PE;
		}
		readTotalNum += realRecordNum;
		log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		log_ofs << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		log_ofs << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;

		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();
			PE_Read_Info readInfo;
			//cout << endl << "read_name_1: " << readName1Vec[tmpOpenMP] << endl;
			//cout << "read_name_2: " << readName2Vec[tmpOpenMP] << endl;			
			//cout << "start to initiateReadInfo ..." << endl;
			readInfo.initiateReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP], readSeq1Vec[tmpOpenMP], 
				readSeq2Vec[tmpOpenMP], readQualSeq1Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], false, false);
			//cout << "start to do segment Map " << endl;  
			SegMap2aluFamily_Info_PeReadBothDir tmpSegMapInfo;
			tmpSegMapInfo.segmentMap(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray, readInfo);
			//cout << "start to generate_segLenVec_segMapScore" << endl;
			tmpSegMapInfo.generate_segLenVec_segMapScore(readInfo);
			//cout << "start to generateMaxSegMapScore" << endl;
			tmpSegMapInfo.generateMaxSegMapScore();
			int tmpSegMapScoreMax = tmpSegMapInfo.returnMaxSegMapScore();
			//cout << "tmpSegMapScoreMax: " << tmpSegMapScoreMax << endl;
			string tmpSegMapInfoStr = tmpSegMapInfo.returnSegLenVecStr();
			//cout << "tmpSegMapInfoStr: " << tmpSegMapInfoStr << endl;
			readSegmentMapScoreVec[tmpOpenMP] = tmpSegMapScoreMax;
			readSegMapInfoStrVec[tmpOpenMP] = tmpSegMapInfoStr;
		}
		
		log_ofs << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{	
			score_ofs << readName1Vec[tmp] << "\t" << readSegmentMapScoreVec[tmp] << endl;
			segLen_ofs << readName1Vec[tmp] << endl << readSegMapInfoStrVec[tmp] << endl;
		}
		log_ofs << "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
	}

	delete readPreProcessInfo;
	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	score_ofs.close();
	segLen_ofs.close();
	log_ofs.close();
	stats_ofs.close();
	return 0;
}