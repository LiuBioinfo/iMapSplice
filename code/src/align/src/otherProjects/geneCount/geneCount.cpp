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
#include "../../general/gene_count_vec.h"
#include "../../general/transcript2geneMap_info.h"
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
#include "../../general/alignInferJunctionHash_info.h"
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

using namespace std;

int main(int argc, char**argv)
{
    bool checkQualSeqForShortAnchorSeqToTargetMap = false;
    cout << "Attention! checkQualSeqForShortAnchorSeqToTargetMap true or not: " << checkQualSeqForShortAnchorSeqToTargetMap << endl;
    bool checkQualSeqForReadSegSeq = false;
	cout << "Attention! checkQualSeqForReadSegSeq true or not: " << checkQualSeqForReadSegSeq << endl;
    /////////////////   get option from command line ////////////////////
    Option_Info* optionInfo = new Option_Info();
    optionInfo->getOpt_long(argc, argv);
    string outputDirStr = optionInfo->outputFolder_path; //argv[3];
   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());	
   	string settingsLogStr = outputDirStr + "/settings.log";
   	ofstream settings_log_ofs(settingsLogStr.c_str());
   	string progressLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(progressLogStr.c_str());
   	string runtimeLogStr = outputDirStr + "/runtime.log";
   	ofstream runtime_log_ofs(runtimeLogStr.c_str());
   	string statsStr = outputDirStr + "/stats.txt";
   	ofstream stats_ofs(statsStr.c_str());
   	optionInfo->outputOptStr(settings_log_ofs);
   	bool SE_or_PE_bool = optionInfo->SE_or_PE_bool;
	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   switches of seperate processes    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
   	bool realGenome_provided_bool = optionInfo->realGenome_provided_bool;
   	string realGenome_file_path = optionInfo->realGenome_file_path;
   	cout << "transcirpt Mode: ref genome: " << realGenome_file_path << endl;
   	bool annotation_provided_bool = false;
   	bool Do_annotation_only_bool = annotation_provided_bool;
	bool Do_cirRNA = true;
	Do_cirRNA = false;
	bool Do_extendHeadTail_phase1 = true;
	Do_extendHeadTail_phase1 = false;
	int normalRecordNum_1stMapping = 500000;
	int readTotalNum = 0;

	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl; 
	//////////////////////////////////////////////        LOAD INDEX         ////////////////////////////////////////////////////////////////
    string InputReadFile = optionInfo->read_file_path_1;//read sample, exacted from fastq file every time
    string InputReadFile_PE = optionInfo->read_file_path_2;// another end read for pair-end reads
	int threads_num = optionInfo->threads_num;//atoi(threadsNumStr.c_str());
	omp_set_num_threads(threads_num);
	bool InputAsFastq = (!(optionInfo->fasta_or_fastq_bool));
	bool fasta_or_fastq_bool = optionInfo->fasta_or_fastq_bool;
	if((checkQualSeqForShortAnchorSeqToTargetMap && fasta_or_fastq_bool)||(checkQualSeqForReadSegSeq && fasta_or_fastq_bool))
	{
		cout << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		log_ofs << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		cout << "if checkQualSeqForReadSegSeq, fasta_or_fastq_bool must be true" << endl;
		log_ofs << "if checkQualSeqForReadSegSeq, fasta_or_fastq_bool must be true" << endl;
		exit(1);
	}
	/////////////////////////////////////          LOAD INDEX         ////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
    cout << "start to load preIndex ..." << endl;
    log_ofs << "start to load preIndex ..." << endl;
    string preIndexArrayPreStr;
    string indexStr;// = "/data/homes/lxauky/adSA_table/mm9_table/testAll2_table/testAll2Index";
    string chromDirStr;
    string secondLevelIndexStr;
    indexStr = optionInfo->global_index_file_path_prefix; //argv[6];
    preIndexArrayPreStr = indexStr;
    chromDirStr = optionInfo->chromsome_file_path_prefix; //argv[8];
    secondLevelIndexStr = optionInfo->local_index_file_path_prefix; //argv[7];
    preIndexArrayPreStr.append("/");
    indexStr.append("/");
    chromDirStr.append("/");
    secondLevelIndexStr.append("/");
	string preIndexMapLengthArrayStr = preIndexArrayPreStr; preIndexMapLengthArrayStr.append("_MapLength"); ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; preIndexIntervalStartArrayStr.append("_IntervalStart"); ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; preIndexIntervalEndArrayStr.append("_IntervalEnd"); ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
 	cout << "finish loading preIndex ..." << endl;
 	log_ofs << "finish loading preIndex ..." << endl;	
	string SA_file = indexStr; SA_file.append("_SA"); ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	string childTab_file = indexStr; childTab_file.append("_childTab"); ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);	
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, settings_log_ofs);
	settings_log_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	cout << "start to load whole genome" << endl;
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	settings_log_ofs << "chromSize = " <<indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;
	/////////////////////////////////////   start to load real genome index  /////////////////////////////////////	
	string genome_chrom_bit_file = realGenome_file_path; genome_chrom_bit_file += "/_chrom"; 
	ifstream genome_chrom_bit_file_ifs(genome_chrom_bit_file.c_str());
	string genome_parameter_file = realGenome_file_path; genome_parameter_file += "/_parameter"; 
	ifstream genome_parameter_file_ifs(genome_parameter_file.c_str());

	Index_Info*  genomeIndexInfo = new Index_Info(genome_parameter_file_ifs, log_ofs);
	log_ofs << "real genome index: " << realGenome_file_path << endl;
	char *genome_chrom; genome_chrom = (char*)malloc((genomeIndexInfo->returnIndexSize()) * sizeof(char));
	genome_chrom_bit_file_ifs.read((char*)genome_chrom, (genomeIndexInfo->returnIndexSize()) * sizeof(char));
	genomeIndexInfo->readGenome(genome_chrom);
	log_ofs << "real genome size = " << genomeIndexInfo->returnChromStringLength() << endl;
	genomeIndexInfo->initiate();
	genomeIndexInfo->initiateChrNameIndexArray(1000);	
	log_ofs << "finish loading real genome" << endl;

	////////////////// transcript info loading .... //////////
	log_ofs << "start to load transcript info" << endl;
	cout << "start to load transcriptInfo " << endl;
	cout << "transcriptInfo: " << optionInfo->transcript_file_path << endl;
	Transcript_Set* transcriptInfo = new Transcript_Set();
	string transcript_file_path = optionInfo->transcript_file_path;
	ifstream transcript_file_ifs(transcript_file_path.c_str());
	string transcript_type = optionInfo->transcript_type;
	transcriptInfo->extractTranscript(transcript_file_ifs, genomeIndexInfo, transcript_type);

	///////////////// gene info loading ...... ///////////////
	log_ofs << "start to load gene info and generate transcriptome2geneMapInfo" << endl;
	cout << "start to load gene info and generate transcriptome2geneMapInfo" << endl;
	string transcript2geneMapFile = optionInfo->local_index_file_path_prefix;
	Transcript2geneMap_Info* transcript2geneMapInfo = new Transcript2geneMap_Info();
	transcript2geneMapInfo->initiate_transcript2geneMap(transcript2geneMapFile);
	int geneNum = transcript2geneMapInfo->returnGeneNum();
	////////////////// initiating geneCount_merged info ..... ///////////
	log_ofs << "start to initiating geneCountInfo_merged ..... " << endl;
	cout << "start to initiating transcriptCountInfo_merged ..... " << endl;
	Gene_Count* geneCountInfo_merged = new Gene_Count();
	geneCountInfo_merged->initiate_withGeneNum(geneNum);
	////////////////// initiating geneCount_vec info .... ///////////////
	log_ofs << "start to initiating geneCountInfo vec  ..... " << endl;
	cout << "start to initiating geneCountInfo vec ..... " << endl;
	Gene_Count_Vec* geneCountInfo_vec = new Gene_Count_Vec();
	geneCountInfo_vec->initiate(threads_num, geneNum);

	/////////////////////////////////////   start to load annotation  /////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load annotation file ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load annotation file ......" << endl << endl; 	
	string annotation_file_path = optionInfo->annotation_file_path;
	ifstream annotation_ifs(annotation_file_path.c_str());
	Annotation_Info* annotationInfo = new Annotation_Info();
	if(annotation_provided_bool)
		annotationInfo->initiateAndReadAnnotationFile(indexInfo, annotation_ifs);
	/////////////////////////////////////   finish loading annotation  /////////////////////////////////////
	cout << "start to load index files" << endl;
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "All index files loaded" << endl;
	cout << "All index files loaded" << endl;
	//////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl;
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	Stats_Info* statsInfo = new Stats_Info();
	statsInfo->initiate_stats_info_PE(threads_num);
	string geneCount_file_str = outputDirStr + "/geneCount.txt";
	string geneCount_otherStats_file_str = outputDirStr + "/geneCount_otherStats.txt";
	ifstream inputRead_ifs(InputReadFile.c_str());
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;
	int normalRecordNum = normalRecordNum_1stMapping; //1000000;//1500000;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;
	int readPairNum = 0;

	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readQualSeq1Vec(normalRecordNum);
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readQualSeq2Vec(normalRecordNum);
	vector< RepeatRegion_Info* > repeatRegionInfoVec;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
		repeatRegionInfoVec.push_back(repeatRegionInfo);
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... gene count process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... gene count process starts ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... gene count process starts ......" << endl << endl; 

	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();
	for(tmpTurn = 0; ; tmpTurn++) // used to control turns to process
	{
		if(EndOfRecord)
			break;		
		int recordNum = normalRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
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
    		readName1Vec[recordNumTmp] = line1.substr(1);
    		getline(inputRead_ifs, line2); // readSeq_1
    		line2_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2);		
    		readSeq1Vec[recordNumTmp] = line2_afterProcess;    		
    		if(InputAsFastq)
    		{
    			getline(inputRead_ifs, line3);
    			getline(inputRead_ifs, line4);
    			readQualSeq1Vec[recordNumTmp] = line4;
    		}
    		getline(inputRead_PE_ifs, line1_PE); // readName_2
    		readName2Vec[recordNumTmp] = line1_PE.substr(1);
    		getline(inputRead_PE_ifs, line2_PE);
    		line2_PE_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2_PE);
    		readSeq2Vec[recordNumTmp] = line2_PE_afterProcess;
    		if(InputAsFastq)
    		{
    			getline(inputRead_PE_ifs, line3_PE);
    			getline(inputRead_PE_ifs, line4_PE);
    			readQualSeq2Vec[recordNumTmp] = line4_PE;
    		}
		}
		readTotalNum += realRecordNum;
		runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		runtime_log_ofs << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to fix reads, turn: " << tmpTurn + 1 << endl;

		omp_set_num_threads(threads_num);
		#pragma omp parallel for
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();
			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP], readSeq1Vec[tmpOpenMP], 
				readSeq2Vec[tmpOpenMP], readQualSeq1Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], 
				fasta_or_fastq_bool, SE_or_PE_bool);
			//cout << endl << "read_name: " << readName1Vec[tmpOpenMP] << endl;
    		FixPhase1Info fixPhase1Info;
			fixPhase1Info.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, SE_or_PE_bool);

			// cout << "\nsegInfo_Nor1: " << endl;
			// fixPhase1Info.coutSegInfo_Nor1(indexInfo);
			// cout << "segInfo_Rcm1: " << endl;
			// fixPhase1Info.coutSegInfo_Rcm1(indexInfo);
			// cout << "segInfo_Nor2: " << endl;
			// fixPhase1Info.coutSegInfo_Nor2(indexInfo);
			// cout << "segInfo_Rcm2: " << endl;
			// fixPhase1Info.coutSegInfo_Rcm2(indexInfo);
			//fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo, annotation_provided_bool, 
			//	Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, SE_or_PE_bool);
			// fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
			// 	annotation_provided_bool, Do_annotation_only_bool, annotationInfo, SE_or_PE_bool);			
			
			#ifdef MAP_INFO
			cout << "output map info: " << endl;
			fixPhase1Info.coutDebugInfo(readInfo, indexInfo, SE_or_PE_bool);
			cout << "finish output map info ..." << endl;
			#endif
			
			//PE_Read_Alignment_Info peAlignInfo_transcript;
			//peAlignInfo_transcript.initiatePeAlignInfo(fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
			//	fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, indexInfo, SE_or_PE_bool);
			
			#ifdef MAP_INFO
			cout << "peAlignInfo: \n" << peAlignInfo.returnPeAlignInfoStr() << endl << endl;
			#endif
			//peAlignInfo_transcript.removeDuplicateMismatch(false);	
			//peAlignInfo_transcript.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();

			//cout << endl << "tmpReadName: " << readName1Vec[tmpOpenMP] << endl;
			//cout << "tmpAlignInfo: " << endl <<	peAlignInfo_transcript.getTmpAlignInfoForFinalPair(
			//	readInfo.returnReadName_1(), readInfo.returnReadName_2(), readInfo.returnReadSeq_1(), readInfo.returnReadSeq_2(),
			//	readInfo.returnReadQual_1(), readInfo.returnReadQual_2(), false, 49) << endl;;
			//cout << "start to do doGeneCount_segMapOnTranscript ....." << endl;
			fixPhase1Info.doGeneCount_segMapOnTranscript(geneCountInfo_vec, threadNO, 
				transcript2geneMapInfo, readInfo.returnReadLength_end1() + readInfo.returnReadLength_end2(), indexInfo);
			//peAlignInfo_transcript.doGeneCount(geneCountInfo_vec, threadNO, transcriptInfo, transcript2geneMapInfo);
			fixPhase1Info.memoryFree();
			//peAlignInfo_transcript.memoryFree();
		} // read file end
		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "finish fixing reads, turn: " << tmpTurn + 1 << endl;
	}
	cout << "finish gene count on total reads" << endl;
	cout << "readTotalNum: " << readTotalNum << endl;
	settings_log_ofs << "readTotalNum: " << readTotalNum << endl;

	geneCountInfo_vec->merge2oneGeneCount(geneCountInfo_merged);
	//geneCountInfo_merged->output(geneSetInfo, geneCount_ofs, geneCount_otherStats_ofs);
	transcript2geneMapInfo->outputGeneCount(geneCountInfo_merged, geneCount_file_str, geneCount_otherStats_file_str);

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	free(preIndexMapLengthArray); free(preIndexIntervalStartArray); free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress); free(childTab); free(verifyChild); free(chrom);
	annotation_ifs.close();
	transcript_file_ifs.close();
	delete geneCountInfo_merged;
	geneCountInfo_vec->memoryFree();
	delete geneCountInfo_vec;
	delete transcript2geneMapInfo;
	transcriptInfo->memoryFree();
	delete transcriptInfo;
	delete annotationInfo;	
	delete indexInfo;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  
    return 0;
} //end main
