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

#ifdef CAL_TIME
clock_t read_file_begin, read_file_end, read_file_end2, read_file_end3, 
		align_begin, align_end, align_cost = 0,
    	overall_begin, overall_end, overall_cost = 0,
    	input_begin, input_end, input_cost = 0,
    	output_begin, output_end, output_cost = 0,
    	segMap_begin, segMap_end, segMap_cost = 0,
    	getPath_begin, getPath_end, getPath_cost = 0,
    	fixGap_begin, fixGap_end, fixGap_cost = 0,
    	
    	fixGap_1_begin, fixGap_1_end, fixGap_1_cost = 0,
    	fixGap_2_begin, fixGap_2_end, fixGap_2_cost = 0,

    	getChrLocation_begin, getChrLocation_end, getChrLocation_cost = 0,

    	getInterval_begin, getInterval_end, getInterval_cost = 0,
    	getFirstInterval_begin, getFirstInterval_end, getFirstInterval_cost = 0,   
    	getLcp_begin_1, getLcp_end_1, getLcp_cost_1 = 0,
    	getLcp_begin_2, getLcp_end_2, getLcp_cost_2 = 0,
       	getLcp_begin_3, getLcp_end_3, getLcp_cost_3 = 0,
       	 	
       	score_string_begin, score_string_end, score_string_cost = 0,
       	checkTwoStringMatch_begin, checkTwoStringMatch_end, checkTwoStringMatch_cost = 0,

    	searchPrefix_begin, searchPrefix_end, searchPrefix_cost = 0,
    	searchInSA_begin, searchInSA_end, searchInSA_cost = 0,
    	insert2SegInfo_begin, insert2SegInfo_end, insert2SegInfo_cost = 0,

    	readFromFile_begin, readFromFile_end, readFromFile_cost = 0,
    	readPreProcess_begin, readPreProcess_end, readPreProcess_cost = 0,

    	getPEalignInfo_begin, getPEalignInfo_end, getPEalignInfo_cost = 0,
    	selectBestAlign_begin, selectBestAlign_end, selectBestAlign_cost = 0,
    	getSamFormat_begin, getSamFormat_end, getSamFormat_cost = 0,
		insertSamFormat_begin, insertSamFormat_end, insertSamFormat_cost = 0,
		getReadInfo_begin, getReadInfo_end, getReadInfo_cost = 0,
		freeMem_begin, freeMem_end, freeMem_cost = 0;
#endif 


//#include "switch.h"
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
//#include "../../general/alignmentToJunc.h"

//#define PreIndexSize 268435456

using namespace std;  


unsigned int PairedReadNum = 0, BothUnmappedReadNum = 0, BothUnmappedReadNum_mappedToRepeatRegion = 0, UnpairedReadNum = 0;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int main(int argc, char**argv)
{
    bool checkQualSeqForShortAnchorSeqToTargetMap = false;
    cout << "Attention! checkQualSeqForShortAnchorSeqToTargetMap true or not: " 
    	<< checkQualSeqForShortAnchorSeqToTargetMap << endl;
    bool checkQualSeqForReadSegSeq = false;
	cout << "Attention! checkQualSeqForReadSegSeq true or not: " 
		<< checkQualSeqForReadSegSeq << endl;

	#ifdef CAL_TIME
    overall_begin = clock();
    #endif
    /////////////////   get option from command line ////////////////////
    Option_Info* optionInfo = new Option_Info();
    optionInfo->getOpt_long(argc, argv);

    //exit(1);
	//////////////////////////////////////////////////
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

	bool Do_Phase1_Only = true;//optionInfo->Do_phase1_only_bool;
	//Do_Phase1_Only = true;

	bool outputAlignInfoAndSamForAllPairedAlignmentBool = false;
	//outputAlignInfoAndSamForAllPairedAlignmentBool = true;

	bool removeAllIntermediateFilesBool = false;
	//removeAllIntermediateFilesBool = true;

	bool DoSam2JuncBool = false;
	DoSam2JuncBool = true;

	bool load2ndLevelIndexBool = false;
	load2ndLevelIndexBool = true;

	bool load2ndLevelIndexBool_compressedSize = false;
	load2ndLevelIndexBool_compressedSize = true;

	bool DoRemappingOnUnmapEndReadsBool = false;
	DoRemappingOnUnmapEndReadsBool = true;

	bool DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	DoRemappingOnUnfixedHeadTailAlignmentBool = true;

	bool outputDirectlyBool_Phase1Only = true;
	outputDirectlyBool_Phase1Only = false;

	bool Do_cirRNA = true;
	Do_cirRNA = false;

	//bool Do_extendHeadTail = true;
	//Do_extendHeadTail = false;

	bool Do_extendHeadTail_phase1 = true;
	//Do_extendHeadTail_phase1 = false;

	bool Do_extendHeadTail_fixOneEndUnmapped = true;
	//Do_extendHeadTail_fixOneEndUnmapped = false;

	bool Do_extendHeadTail_fixHeadTail = true;
	//Do_extendHeadTail_fixHeadTail = false;	

	
 
	bool Do_fixHeadTail_remapping = true;
	//Do_fixHeadTail_remapping = false;

	bool Do_fixHeadTail_greedyMapping = true;
	//Do_fixHeadTail_greedyMapping = false;

	bool Do_fixHeadTail_remappingAndTargetMapping = true;
	//Do_fixHeadTail_remappingAndTargetMapping = false;

	bool Do_fixHeadTail_remappingAgain = true;
	Do_fixHeadTail_remappingAgain = false;


	if(Do_Phase1_Only)
	{
		settings_log_ofs << "Do_Phase1 only!" << endl;
		DoSam2JuncBool = false;
		load2ndLevelIndexBool = false;
		load2ndLevelIndexBool_compressedSize = false;
		DoRemappingOnUnmapEndReadsBool = false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	}	
	else
	{
		settings_log_ofs << "Do_Phase1_Phase2! " << endl;
		DoSam2JuncBool = true;//false;
		load2ndLevelIndexBool = true;//false;
		load2ndLevelIndexBool_compressedSize = true;//false;
		DoRemappingOnUnmapEndReadsBool = true;//false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = true;//false;
	}

	int normalRecordNum_1stMapping = 10000;//2000000;
	int normalRecordNum_fixOneEndUnmapped = 2000000;
	int normalRecordNum_fixHeadTail = 2000000;

	int readTotalNum = 0;

	optionInfo->outputSwitchInfo(Do_Phase1_Only, outputAlignInfoAndSamForAllPairedAlignmentBool,
		removeAllIntermediateFilesBool, Do_cirRNA, outputDirectlyBool_Phase1Only, 
		normalRecordNum_1stMapping, normalRecordNum_fixOneEndUnmapped,
		normalRecordNum_fixHeadTail, Do_extendHeadTail_phase1, 
		Do_extendHeadTail_fixOneEndUnmapped, Do_extendHeadTail_fixHeadTail, 
		Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
		Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, settings_log_ofs);
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl; 

	//////////////////////////////////////////////        LOAD INDEX         ////////////////////////////////////////////////////////////////

    string InputReadFile = optionInfo->read_file_path_1;//read sample, exacted from fastq file every time
    string InputReadFile_PE = optionInfo->read_file_path_2;// another end read for pair-end reads

	//string threadsNumStr = argv[4];
	int threads_num = optionInfo->threads_num;//atoi(threadsNumStr.c_str());

	#ifdef CAL_TIME	
	threads_num = 1;
	#endif

	#ifdef MAP_INFO
	threads_num = 1;
	#endif

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
	/////////////////////////////////////          LOAD INDEX         ////////////////////////////////////////////////////////////////
	#ifdef CAL_TIME
	read_file_begin = clock();
	#endif
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 

    //cout << "start to load preIndex ..." << endl;
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
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	settings_log_ofs << "chromSize = " <<indexInfo->returnChromStringLength() << endl;
	log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	log_ofs << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;

	/////////////////////////////////////   start to load real genome index  /////////////////////////////////////
	
	string genome_chrom_bit_file = realGenome_file_path; genome_chrom_bit_file += "/_chrom"; 
	ifstream genome_chrom_bit_file_ifs(genome_chrom_bit_file.c_str());
	string genome_parameter_file = realGenome_file_path; genome_parameter_file += "/_parameter"; 
	ifstream genome_parameter_file_ifs(genome_parameter_file.c_str());

	Index_Info*  genomeIndexInfo = new Index_Info(genome_parameter_file_ifs, log_ofs);
	log_ofs << "real genome index: " << realGenome_file_path << endl;
	log_ofs << "start to load whole real genome " << endl;
	char *genome_chrom; genome_chrom = (char*)malloc((genomeIndexInfo->returnIndexSize()) * sizeof(char));
	genome_chrom_bit_file_ifs.read((char*)genome_chrom, (genomeIndexInfo->returnIndexSize()) * sizeof(char));
	genomeIndexInfo->readGenome(genome_chrom);
	log_ofs << "real genome size = " << genomeIndexInfo->returnChromStringLength() << endl;
	log_ofs << "start to initiate and load every chromosome for ref genome " << endl;
	genomeIndexInfo->initiate();
	log_ofs << "start to initiate real genome chrNameIndexArray" << endl;
	genomeIndexInfo->initiateChrNameIndexArray(1000);	
	log_ofs << "finish loading real genome" << endl;

	////////////////// transcript info loading .... /////
	log_ofs << "start to load transcript info" << endl;
	cout << "start to load transcriptInfo " << endl;
	settings_log_ofs << "transcriptInfo: " << optionInfo->transcript_file_path << endl;
	cout << "transcriptInfo: " << optionInfo->transcript_file_path << endl;
	Transcript_Set* transcriptInfo = new Transcript_Set();
	string transcript_file_path = optionInfo->transcript_file_path;
	ifstream transcript_file_ifs(transcript_file_path.c_str());
	string transcript_type = optionInfo->transcript_type;
	transcriptInfo->extractTranscript(transcript_file_ifs, genomeIndexInfo, transcript_type);
	/////////////////////////////////////   start to load annotation  /////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load annotation file ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load annotation file ......" << endl << endl; 	
	string annotation_file_path = optionInfo->annotation_file_path;
	ifstream annotation_ifs(annotation_file_path.c_str());
	Annotation_Info* annotationInfo = new Annotation_Info();
	if(annotation_provided_bool)
	{
		//annotation_ifs.open(annotation_file_path);
		annotationInfo->initiateAndReadAnnotationFile(indexInfo, annotation_ifs);
	}
	/////////////////////////////////////   finish loading annotation  /////////////////////////////////////
	log_ofs << "start to load SA" << endl;
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load lcpCompress" << endl;
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));	
	log_ofs << "start to load childTab " << endl;
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load detChild" << endl;
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "All index files loaded" << endl;
	runtime_log_ofs << "All index files loaded" << endl;
	#ifdef CAL_TIME
	read_file_end = clock();
	double read_file_time = (double)(read_file_end - read_file_begin)/CLOCKS_PER_SEC;
	log_ofs << "read_file cpu time = " << read_file_time << endl;
	#endif
	//////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... whole genome index loaded ......" << endl << endl;

	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	
	HeaderSection_Info* headerSectionInfo = new HeaderSection_Info(genomeIndexInfo);	

	Stats_Info* statsInfo = new Stats_Info();
	statsInfo->initiate_stats_info_PE(threads_num);

	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
	/* align main*/
   	string tmpHeadSectionInfo = outputDirStr + "/headSectionInfo";
   	ofstream tmpHeadSectionInfo_ofs(tmpHeadSectionInfo.c_str());    	
   	tmpHeadSectionInfo_ofs << headerSectionInfo->returnHeaderSectionInfoStr() << endl;
	tmpHeadSectionInfo_ofs.close();

	string mkdirOutputCommand_phase1 = "mkdir -p " + outputDirStr + "/phase1_output";
	system(mkdirOutputCommand_phase1.c_str());

	string mkdirOutputCommand_repeatRegionFile = mkdirOutputCommand_phase1 + "/repeat_region";
   	system(mkdirOutputCommand_repeatRegionFile.c_str());
	string repeatRegionFile = outputDirStr + "/phase1_output/repeat_region/repeatRegion";
	ofstream repeatRegionFile_ofs(repeatRegionFile.c_str());

	string mkdirOutputCommand_tmpAlignCompleteRead = mkdirOutputCommand_phase1 + "/completePair";
	system(mkdirOutputCommand_tmpAlignCompleteRead.c_str());
	string tmpAlignCompleteRead = outputDirStr + "/phase1_output/completePair/completePair.sam";
	ofstream tmpAlignCompleteRead_ofs(tmpAlignCompleteRead.c_str());
	string tmpAlignCompleteRead_alignInfo = outputDirStr + "/phase1_output/completePair/completePair.sam_alignInfo";
	ofstream tmpAlignCompleteRead_alignInfo_ofs(tmpAlignCompleteRead_alignInfo.c_str());

	#ifdef CHECK_MULTI
	string tmpAlignCompleteRead_multi = outputDirStr + "/phase1_output/completePair/completePair.sam.multi";
	ofstream tmpAlignCompleteRead_multi_ofs(tmpAlignCompleteRead_multi.c_str());
	string tmpAlignCompleteRead_unique = outputDirStr + "/phase1_output/completePair/completePair.sam.unique";
	ofstream tmpAlignCompleteRead_unique_ofs(tmpAlignCompleteRead_unique.c_str());
	#endif

	string mkdirOutputCommand_tmpAlignOneEndUnmapped = mkdirOutputCommand_phase1 + "/oneEndUnmapped";
	system(mkdirOutputCommand_tmpAlignOneEndUnmapped.c_str());
	string tmpAlignOneEndUnmapped = outputDirStr + "/phase1_output/oneEndUnmapped/oneEndUnmapped";
	if(Do_Phase1_Only)
		tmpAlignOneEndUnmapped += ".sam";	
	else
		tmpAlignOneEndUnmapped += ".alignInfo";
	ofstream tmpAlignOneEndUnmapped_ofs(tmpAlignOneEndUnmapped.c_str());


	//string mkdirOutputCommand_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = mkdirOutputCommand_phase1 + "/bothEndsUnmapped_mappedToRepeatRegion";
	//system(mkdirOutputCommand_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
	string tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile 
		= outputDirStr + "/phase1_output/repeat_region/bothEndsUnmapped_mappedToRepeatRegion.sam"; 
	ofstream tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs(
		tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());

	string mkdirOutputCommand_tmpAlignBothEndsUnmapped = mkdirOutputCommand_phase1 + "/bothEndsUnmapped";
	system(mkdirOutputCommand_tmpAlignBothEndsUnmapped.c_str());
	string tmpAlignBothEndsUnmapped = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped.sam";
	ofstream tmpAlignBothEndsUnmapped_ofs(tmpAlignBothEndsUnmapped.c_str());
	string tmpAlignBothEndsUnmapped_lowScore = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped_lowScore.sam";
	ofstream tmpAlignBothEndsUnmapped_lowScore_ofs(tmpAlignBothEndsUnmapped_lowScore.c_str());

	string mkdirOutputCommand_tmpAlignIncompletePair = mkdirOutputCommand_phase1 + "/incomplete";
	system(mkdirOutputCommand_tmpAlignIncompletePair.c_str());
	string tmpAlignIncompletePair = outputDirStr + "/phase1_output/incomplete/incomplete.alignInfo"; 
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());
	string tmpAlignIncompletePair_SAM = outputDirStr + "/phase1_output/incomplete/incompletePair.sam"; 
	ofstream tmpAlignIncompletePair_SAM_ofs(tmpAlignIncompletePair_SAM.c_str());	

	string tmpIntermediateJunctionFile = outputDirStr + "/phase2_output/inter.junc";
	/*
		vector< ofstream* > tmpAlignCompleteRead_ofs_vec; 
		vector< ofstream* > tmpAlignOneEndUnmapped_ofs_vec;
		vector< ofstream* > tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs_vec;
		vector< ofstream* > tmpAlignBothEndsUnmapped_ofs_vec;
		vector< ofstream* > tmpAlignIncompletePair_ofs_vec;
		vector< ofstream* > tmpAlignIncompletePair_SAM_ofs_vec;
		vector< ofstream* > tmpAlignCompleteRead_alignInfo_ofs_vec;
	
	#ifdef CHECK_MULTI
	vector< ofstream* > tmpAlignCompleteRead_ofs_vec_multi; 
	vector< ofstream* > tmpAlignCompleteRead_ofs_vec_unique; 
	#endif

	if(outputDirectlyBool_Phase1Only)
	{  
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			string tmp_tmpAlignCompleteRead = tmpAlignCompleteRead + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignCompleteRead_ofs = new ofstream(tmp_tmpAlignCompleteRead.c_str());
			tmpAlignCompleteRead_ofs_vec.push_back(tmp_tmpAlignCompleteRead_ofs);

			string tmp_tmpAlignOneEndUnmapped = tmpAlignOneEndUnmapped + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignOneEndUnmapped_ofs = new ofstream(tmp_tmpAlignOneEndUnmapped.c_str());
			tmpAlignOneEndUnmapped_ofs_vec.push_back(tmp_tmpAlignOneEndUnmapped_ofs);

			string tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs = new ofstream(tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
			tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs_vec.push_back(tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs);

			string tmp_tmpAlignBothEndsUnmapped = tmpAlignBothEndsUnmapped + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignBothEndsUnmapped_ofs = new ofstream(tmp_tmpAlignBothEndsUnmapped.c_str());
			tmpAlignBothEndsUnmapped_ofs_vec.push_back(tmp_tmpAlignBothEndsUnmapped_ofs);

			string tmp_tmpAlignIncompletePair = tmpAlignIncompletePair + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignIncompletePair_ofs = new ofstream(tmp_tmpAlignIncompletePair.c_str());
			tmpAlignIncompletePair_ofs_vec.push_back(tmp_tmpAlignIncompletePair_ofs);

			string tmp_tmpAlignIncompletePair_SAM = tmpAlignIncompletePair_SAM + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignIncompletePair_SAM_ofs = new ofstream(tmp_tmpAlignIncompletePair_SAM.c_str());
			tmpAlignIncompletePair_SAM_ofs_vec.push_back(tmp_tmpAlignIncompletePair_SAM_ofs);

			string tmp_tmpAlignCompleteRead_alignInfo = tmpAlignCompleteRead_alignInfo + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignCompleteRead_alignInfo_ofs = new ofstream(tmp_tmpAlignCompleteRead_alignInfo.c_str());
			tmpAlignCompleteRead_alignInfo_ofs_vec.push_back(tmp_tmpAlignCompleteRead_alignInfo_ofs);

			#ifdef CHECK_MULTI
			string tmp_tmpAlignCompleteRead_multi = tmpAlignCompleteRead_multi + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignCompleteRead_multi_ofs = new ofstream(tmp_tmpAlignCompleteRead_multi.c_str());
			tmpAlignCompleteRead_ofs_vec_multi.push_back(tmp_tmpAlignCompleteRead_multi_ofs);
			
			string tmp_tmpAlignCompleteRead_unique = tmpAlignCompleteRead_unique + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignCompleteRead_unique_ofs = new ofstream(tmp_tmpAlignCompleteRead_unique.c_str());
			tmpAlignCompleteRead_ofs_vec_unique.push_back(tmp_tmpAlignCompleteRead_unique_ofs);
			#endif 
		}
	}*/

	//*(tmpAlignCompleteRead_ofs_vec[1]) << "testing" << endl;

	ifstream inputRead_ifs(InputReadFile.c_str());
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;

	int normalRecordNum = normalRecordNum_1stMapping; //1000000;//1500000;

	#ifdef DEBUG_INFO
	normalRecordNum = 1;
	#endif

	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;

	int readPairNum = 0;

	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readQualSeq1Vec(normalRecordNum);
	//vector<string> readSeq1Vec_RC(normalRecordNum);
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readQualSeq2Vec(normalRecordNum);

	vector<string> PeAlignSamStrVec_complete(normalRecordNum);

	vector<string> PeAlignInfoStrVec_inCompletePair(normalRecordNum);

	vector<string> PeAlignInfoStrVec_oneEndUnmapped(normalRecordNum);

	vector<string> PeAlignSamStrVec_bothEndsUnmapped(normalRecordNum);

	vector<string> PeAlignSamStrVec_bothEndsUnmapped_lowScore(normalRecordNum);	

	vector<string> PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion(normalRecordNum);

	vector<string> PeAlignSamStrVec_inCompletePair(normalRecordNum);

	vector<string> PeAlignInfoStrVec_completePaired(normalRecordNum);
	vector< RepeatRegion_Info* > repeatRegionInfoVec;

	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
		repeatRegionInfoVec.push_back(repeatRegionInfo);
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
	int readLengthMax_tmp = 0;

	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();

	int perfectMatch_pair = 0;

	for(tmpTurn = 0; 
		//tmpTurn <= 300     //used to control # of rounds to process
		; tmpTurn++)
	{
		#ifdef CAL_TIME		
		input_begin = clock();
		#endif
		if(EndOfRecord)
			break;		
		int recordNum = normalRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
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
    		//readSeq1Vec[recordNumTmp] = line2;

    		int readLength_1 = line2.length();

			if(readLength_1 > readLengthMax_tmp)
				readLengthMax_tmp = readLength_1;
			
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
    		//readSeq2Vec[recordNumTmp] = line2_PE;

    		int readLength_2 = line2_PE.length();

    		if(readLength_2 > readLengthMax_tmp)
    			readLengthMax_tmp = readLength_2;
    		
    		if(InputAsFastq)
    		{
    			getline(inputRead_PE_ifs, line3_PE);
    			getline(inputRead_PE_ifs, line4_PE);
    			readQualSeq2Vec[recordNumTmp] = line4_PE;
    		}
		}
		readTotalNum += realRecordNum;
		#ifdef CAL_TIME	
		input_end = clock();
		input_cost = input_cost + input_end - input_begin;
		#endif
		runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;
		#ifdef CAL_TIME	
		align_begin = clock();
		#endif

		omp_set_num_threads(threads_num);

		#pragma omp parallel for
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			/*readPairNum ++;
			if(readPairNum == 831)  
			{

			}
			else
			{
				continue;
			}*/
		
			int threadNO = omp_get_thread_num();

			#ifdef CAL_TIME
			getReadInfo_begin = clock();
			#endif

			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP],
				readSeq1Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP],
				readQualSeq1Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], fasta_or_fastq_bool, SE_or_PE_bool);
			//cout << "read_name: " << readName1Vec[tmpOpenMP] << endl;
			#ifdef CAL_TIME 
			getReadInfo_end = clock();
			getReadInfo_cost = getReadInfo_cost + getReadInfo_end - getReadInfo_begin;

    		segMap_begin = clock();
    		#endif
    		
			//FixPhase1Info* fixPhase1Info = new FixPhase1Info();
    		FixPhase1Info fixPhase1Info;

			fixPhase1Info.fixPhase1_segInfo(
				sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, SE_or_PE_bool);

			#ifdef CAL_TIME
			segMap_end = clock();
			segMap_cost = segMap_cost + segMap_end - segMap_begin;
			#endif

			#ifdef CAL_TIME
			getPath_begin = clock();
			#endif

			fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, 
				MAX_SPLICE_DISTANCE_PHASE1, SE_or_PE_bool);

			#ifdef CAL_TIME
			getPath_end = clock();
			getPath_cost = getPath_cost + getPath_end - getPath_begin;

			fixGap_begin = clock();
			#endif

			fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, SE_or_PE_bool);
			
			#ifdef MAP_INFO
			cout << "output map info: " << endl;
			fixPhase1Info.coutDebugInfo(readInfo, indexInfo, SE_or_PE_bool);
			cout << "finish output map info ..." << endl;
			#endif

			#ifdef CAL_TIME
			fixGap_end = clock();
			fixGap_cost = fixGap_cost + fixGap_end - fixGap_begin;
			
			getPEalignInfo_begin = clock();
			#endif

			PE_Read_Alignment_Info peAlignInfo_transcript;
			peAlignInfo_transcript.initiatePeAlignInfo(
				fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
				fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, indexInfo, SE_or_PE_bool);
			#ifdef MAP_INFO
			cout << "transcript_peAlignInfo: \n" << peAlignInfo_transcript.returnPeAlignInfoStr() << endl << endl;
			cout << "finish outputing transcript_peAlignInfo" << endl;
			#endif			
			PE_Read_Alignment_Info peAlignInfo;
			peAlignInfo.convertTranscriptAlignInfo2GenomeAlignInfo(peAlignInfo_transcript,
				transcriptInfo, genomeIndexInfo);
			#ifdef MAP_INFO
			cout << "peAlignInfo: \n" << peAlignInfo.returnPeAlignInfoStr() << endl << endl;
			#endif
			#ifdef CAL_TIME
			getPEalignInfo_end = clock();
			getPEalignInfo_cost = getPEalignInfo_cost + getPEalignInfo_end - getPEalignInfo_begin;	
			#endif
			#ifdef CAL_TIME
			selectBestAlign_begin = clock();
			#endif

			peAlignInfo.removeDuplicateMismatch(false);
			if(Do_Phase1_Only)
			{			
				peAlignInfo.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
			}
			else
			{
				peAlignInfo.alignmentFilter_fixPhase1_SJpenalty(readInfo.returnReadSeqLength_1(),
					readInfo.returnReadSeqLength_2());
			}
			if(!outputDirectlyBool_Phase1Only)
			{	
				PeAlignSamStrVec_complete[tmpOpenMP] = "";			
				PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = "";
				PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = "";
				PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = "";
				PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmpOpenMP] = "";
				PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmpOpenMP] = "";
				PeAlignSamStrVec_inCompletePair[tmpOpenMP] = "";
			}
			#ifdef CAL_TIME
			selectBestAlign_end = clock();
			selectBestAlign_cost = selectBestAlign_cost + selectBestAlign_end - selectBestAlign_begin;	
			#endif		
			
			#ifdef CAL_TIME
			getSamFormat_begin = clock();
			#endif

			peAlignInfo.output_phase1(
				outputDirectlyBool_Phase1Only,
				Do_Phase1_Only,
				PeAlignSamStrVec_complete,
				PeAlignInfoStrVec_inCompletePair,
				PeAlignInfoStrVec_oneEndUnmapped,
				PeAlignSamStrVec_bothEndsUnmapped,
				PeAlignSamStrVec_bothEndsUnmapped_lowScore,
				PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
				PeAlignSamStrVec_inCompletePair,
				//PeAlignInfoStrVec_completePaired,
				repeatRegionInfoVec,
				readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo,
				fasta_or_fastq_bool);


			#ifdef CAL_TIME
			getSamFormat_end = clock();
			getSamFormat_cost = getSamFormat_cost + getSamFormat_end - getSamFormat_begin;

			freeMem_begin = clock();
			#endif

			fixPhase1Info.memoryFree();
			//delete fixPhase1Info;
			//delete readInfo; 
			peAlignInfo_transcript.memoryFree();
			peAlignInfo.memoryFree();
			//delete peAlignInfo;

			#ifdef CAL_TIME
			freeMem_end = clock();
			freeMem_cost = freeMem_cost + freeMem_end - freeMem_begin;
			#endif

		} // read file end
		
		#ifdef CAL_TIME
		align_end = clock();
		align_cost = align_cost + align_end - align_begin;
		
		output_begin = clock();
		#endif
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
		
		if(!outputDirectlyBool_Phase1Only)
		{
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{	
				if(PeAlignSamStrVec_complete[tmp] != "")
				{
					tmpAlignCompleteRead_ofs << PeAlignSamStrVec_complete[tmp] << endl;
					PairedReadNum ++;
					//PeAlignSamStrVec_complete[tmp] = "";
				}			

				if(PeAlignInfoStrVec_inCompletePair[tmp] != "")
				{
					tmpAlignIncompletePair_ofs << PeAlignInfoStrVec_inCompletePair[tmp] << endl;
					//PeAlignInfoStrVec_inCompletePair[tmp] = "";
				}
				
				if(PeAlignInfoStrVec_oneEndUnmapped[tmp] != "")
				{
					tmpAlignOneEndUnmapped_ofs << PeAlignInfoStrVec_oneEndUnmapped[tmp] << endl;
					//PeAlignInfoStrVec_oneEndUnmapped[tmp] = "";
				}
				
				if(PeAlignSamStrVec_bothEndsUnmapped[tmp] != "")
				{
					tmpAlignBothEndsUnmapped_ofs << PeAlignSamStrVec_bothEndsUnmapped[tmp] << endl;
					//PeAlignSamStrVec_bothEndsUnmapped[tmp] = "";
					BothUnmappedReadNum ++;
				}
				if(PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] != "")
				{
					tmpAlignBothEndsUnmapped_lowScore_ofs << PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] << endl;
				}
				if(PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] != "")
				{
					tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs << PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] << endl;
					BothUnmappedReadNum_mappedToRepeatRegion ++;
				}

				if(PeAlignSamStrVec_inCompletePair[tmp] != "")
				{
					tmpAlignIncompletePair_SAM_ofs << PeAlignSamStrVec_inCompletePair[tmp] << endl;
				}

			}
		}
		#ifdef CAL_TIME
		output_end = clock();
		output_cost = output_cost + output_end - output_begin;
		#endif
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
	}
	//log_ofs << "perfectMatch_pair #: " << perfectMatch_pair << endl;
	repeatRegionFile_ofs << "Repeat Region Info: size = " << repeatRegionInfoVec.size() << endl;
	for(int tmpThread = 0; tmpThread < threads_num; tmpThread++)
	{
		repeatRegionInfoVec[tmpThread]->outputRepeatRegion(tmpThread+1, indexInfo, sa, 100, repeatRegionFile_ofs);
		//repeatRegionInfoVec[tmpThread]->outputRepeatRegion(esaInfo, tmpThread+1, indexInfo, //sa, 
		//	100, repeatRegionFile_ofs);
	}


	settings_log_ofs << "readTotalNum: " << readTotalNum << endl;

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	tmpAlignCompleteRead_ofs.close();
	//tmpAlignIncompletePair_ofs.close();
	tmpAlignOneEndUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs.close();
	tmpAlignIncompletePair_SAM_ofs.close();
	if(Do_Phase1_Only)
	{
		tmpAlignIncompletePair_ofs.close();
	}

	//fclose(fp_in);

	free(preIndexMapLengthArray); free(preIndexIntervalStartArray); free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);
	free(verifyChild);

	free(chrom);
	//esaInfo->freeMem();
	//delete esaInfo;
	
	#ifdef CAL_TIME
	overall_end = clock();
	#endif
	
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 

	log_ofs << endl << "**********************************" << endl << "**********************************";
	runtime_log_ofs << endl << "**********************************" << endl << "**********************************";


	#ifdef CAL_TIME
	//if(threads_num == 1)
	//{	
		double overall_time = (double)(overall_end - overall_begin)/CLOCKS_PER_SEC;
		double input_time = (double)input_cost/CLOCKS_PER_SEC;
		double align_time = (double)align_cost/CLOCKS_PER_SEC;
		double ouput_time = (double)output_cost/CLOCKS_PER_SEC;

		double getReadInfo_time = (double)getReadInfo_cost/CLOCKS_PER_SEC;
		double segMap_time = (double)segMap_cost/CLOCKS_PER_SEC;
		double getPath_time = (double)getPath_cost/CLOCKS_PER_SEC;
		double fixGap_time = (double)fixGap_cost/CLOCKS_PER_SEC;

		double fixGap_1_time = (double)fixGap_1_cost/CLOCKS_PER_SEC;
		double fixGap_2_time = (double)fixGap_2_cost/CLOCKS_PER_SEC;

		double getChrLocation_time = (double)getChrLocation_cost/CLOCKS_PER_SEC;

		double getPEalignInfo_time = (double)getPEalignInfo_cost/CLOCKS_PER_SEC;
		double selectBestAlign_time = (double)selectBestAlign_cost/CLOCKS_PER_SEC;
		double getSamFormat_time = (double)getSamFormat_cost/CLOCKS_PER_SEC;
		double freeMem_time = (double)freeMem_cost/CLOCKS_PER_SEC;

		double getInterval_time = (double)getInterval_cost/CLOCKS_PER_SEC;
		double getFirstInterval_time = (double)getFirstInterval_cost/CLOCKS_PER_SEC;
		double getLcp_1_time = (double)getLcp_cost_1/CLOCKS_PER_SEC;
		double getLcp_2_time = (double)getLcp_cost_2/CLOCKS_PER_SEC;
		double getLcp_3_time = (double)getLcp_cost_3/CLOCKS_PER_SEC;

		double searchPrefix_time = (double)searchPrefix_cost/CLOCKS_PER_SEC;
		double searchInSA_time = (double)searchInSA_cost/CLOCKS_PER_SEC;		
		double insert2SegInfo_time = (double)insert2SegInfo_cost/CLOCKS_PER_SEC;

		double readFromFile_time = (double)readFromFile_cost/CLOCKS_PER_SEC;
		double readPreProcess_time = (double)readPreProcess_cost/CLOCKS_PER_SEC;

		double score_string_time = (double)score_string_cost/CLOCKS_PER_SEC;
		double checkTwoStringMatch_time = (double)checkTwoStringMatch_cost/CLOCKS_PER_SEC;

		log_ofs << endl << "overall_time = " << overall_time << endl;
		log_ofs << endl << "input_time = " << input_time << endl;
		log_ofs << endl << "align_time = " << align_time << endl << endl;
		
		log_ofs << endl << "getReadInfo_time = " << getReadInfo_time << endl;
		log_ofs << endl << "segMap_time = " << segMap_time << endl;
		log_ofs << endl << "getPath_time = " << getPath_time << endl;
		log_ofs << endl << "fixGap_time = " << fixGap_time << endl;
		
		log_ofs << endl << "fixGap_1_time = " << fixGap_1_time << endl;
		log_ofs << endl << "fixGap_2_time = " << fixGap_2_time << endl;		

		log_ofs << endl << "getChrLocation_time = " << getChrLocation_time << endl;

		log_ofs << endl << "getPEalignInfo_time = " << getPEalignInfo_time << endl;
		log_ofs << endl << "selectBestAlign_time = " << selectBestAlign_time << endl;
		log_ofs << endl << "getSamFormat_time = " << getSamFormat_time << endl;
		log_ofs << endl << "freeMem_time = " << freeMem_time << endl;
		log_ofs << endl << endl << "ouput_time = " << ouput_time << endl;

		log_ofs << endl << "getInterval_time = " << getInterval_time << endl;
		log_ofs << endl << "getFirstInterval_time = " << getFirstInterval_time << endl;

		log_ofs << endl << "getLcp_time_1 = " << getLcp_1_time << endl;
		log_ofs << endl << "getLcp_time_2 = " << getLcp_2_time << endl;
		log_ofs << endl << "getLcp_time_3 = " << getLcp_3_time << endl;

		log_ofs << endl << "searchPrefix_time = " << searchPrefix_time << endl;
		log_ofs << endl << "searchInSA_time = " << searchInSA_time << endl;
		log_ofs << endl << "insert2SegInfo_time = " << insert2SegInfo_time << endl;

		log_ofs << endl << "readFromFile_time = " << readFromFile_time << endl;
		log_ofs << endl << "readPreProcess_time = " << readPreProcess_time << endl;		

		log_ofs << endl << "score_string_time = " << score_string_time << endl;
		log_ofs << endl << "checkTwoStringMatch_time = " << checkTwoStringMatch_time << endl;

	//}
	//else
	//{
	if(threads_num != 1)
		log_ofs << "!!! set thread=1 to get the time-cost information ... !!!" << endl;
		//exit(1);
	//}
	#endif

	log_ofs << endl << "**********************************" << endl << "**********************************" << endl;
	//cout << endl << "totalReadNum = " << read_num << endl;
	/*
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load 2nd level index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index starts ......" << endl << endl; 

	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;

	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;

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
		//loadIndex_end = clock(); 
	}

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ; 		
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ; 	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On one end unmapped Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << endl << "[" << asctime(local) 
		<< "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 	 
	log_ofs << endl << endl << "[" << asctime(local) 
		<< "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 
	runtime_log_ofs << endl << endl << "[" << asctime(local) 
		<< "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 


	string mkdirOutputCommand_phase2 = "mkdir -p " + outputDirStr + "/phase2_output";
	system(mkdirOutputCommand_phase2.c_str());

	string OutputSamFile_oneEndMapped = outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam";
	ofstream OutputSamFile_oneEndMapped_ofs(OutputSamFile_oneEndMapped.c_str());

	string OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore = outputDirStr + "/phase2_output/oneEndUnmapped.bothEndsUnmapped_lowScore.sam";
	ofstream OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());

	string OutputSamFile_oneEndMapped_unpair = outputDirStr + "/phase2_output/oneEndUnmapped.unpaired.sam";
	ofstream OutputSamFile_oneEndMapped_unpair_ofs(OutputSamFile_oneEndMapped_unpair.c_str());

	string OutputSamFile_oneEndMapped_alignInfo = outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam_alignInfo";
	ofstream OutputSamFile_oneEndMapped_alignInfo_ofs(OutputSamFile_oneEndMapped_alignInfo.c_str());	


	int tmpRecordNum_oneEndUnmapped = 0;

	if(DoRemappingOnUnmapEndReadsBool)
	{
		nowtime = time(NULL);
		local = localtime(&nowtime);
		log_ofs << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;
		runtime_log_ofs << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;

		string oneEndMappedFileStr = tmpAlignOneEndUnmapped;

		ifstream inputRecord_ifs(oneEndMappedFileStr.c_str());

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;
		
		int normalRecordNum = normalRecordNum_fixOneEndUnmapped; //1000000;

		bool EndOfRecord = false;

		int tmpTurn = 0;

		int realRecordNum;// = normalRecordNum;

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

		//getline(inputRecord_ifs, line11);

		for(tmpTurn = 0; tmpTurn < TurnNum; tmpTurn++)
		{
			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << endl << "[" << asctime(local) 
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
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local) 
				<< "finish reading record, turn: " << tmpTurn+1 << endl;
			runtime_log_ofs << endl << "[" << asctime(local) 
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
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
					MAX_SPLICE_DISTANCE_PHASE2,
					checkQualSeqForReadSegSeq);

				//peAlignInfo->doNoPairing();
				//peAlignInfo->chooseAllPairedAlignment();
				//peAlignInfo->chooseBestAlignment_selectRandomOneIfMulti();
				//peAlignInfo->chooseBestAlignment_selectAllIfMultiBest();
				//peAlignInfo->chooseBestAlignment_selectRandomOneIfMultiBest_filterOutNoncanonical();
				//peAlignInfo->chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical();
				peAlignInfo->alignmentFilter_fixOneEndUnmapped_SJpenalty(peReadInfo.returnReadSeqLength_1(),
					peReadInfo.returnReadSeqLength_2());

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				bool allUnpairedAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();

				#ifdef NO_FILTERING
				pairExistsBool = false;
				allAlignmentCompleteBool = false;
				allUnpairedAlignmentCompleteBool = false;
				#endif

				bool unique_bool = peAlignInfo->checkUniqueOrMulti();

				statsInfo->increNum_fixUnpaired(threadNO, pairExistsBool, allAlignmentCompleteBool, 
					allUnpairedAlignmentCompleteBool, unique_bool);

				string tmpPeAlignSamStr, tmpPeAlignInfoStr, tmpPeAlignSamStr_unpair, tmpPeAlignSamStr_bothEndsUnmapped_lowScore;

				if(pairExistsBool && allAlignmentCompleteBool) // some pair exists, all completed, print out paired SAM info
				{
					int alignment_score_min_output 
						= peReadInfo.returnAlignmentScoreMinOutput_withComplement(
							ALIGNMENT_SCORE_MIN_OUTPUT_COMPLEMENT);
					bool completeAlign_lowScore_bool
						= peAlignInfo->alignPairScoreTooLow_bool(alignment_score_min_output); 
					if(completeAlign_lowScore_bool)
					{
						statsInfo->increLowScoreComplete_fixUnpair(threadNO, unique_bool);
						tmpPeAlignSamStr = "";
						tmpPeAlignInfoStr = "";
						tmpPeAlignSamStr_unpair = "";
						tmpPeAlignSamStr_bothEndsUnmapped_lowScore = peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool);					
					}
					else
					{	
						tmpPeAlignSamStr = peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo,
							fasta_or_fastq_bool);
						tmpPeAlignInfoStr = "";
						tmpPeAlignSamStr_unpair = "";
						tmpPeAlignSamStr_bothEndsUnmapped_lowScore = "";
					}
				}
				else if(pairExistsBool && (!allAlignmentCompleteBool)) // pair exists, incomplete
				{
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = peAlignInfo->getTmpAlignInfoForFinalPair(
						 		peReadInfo.returnReadName_1(), peReadInfo.returnReadName_2(), 
								peReadInfo.returnReadSeq_1(), peReadInfo.returnReadSeq_2(),
								peReadInfo.returnReadQual_1(), peReadInfo.returnReadQual_2(),
								fasta_or_fastq_bool);
					tmpPeAlignSamStr_unpair = "";
					tmpPeAlignSamStr_bothEndsUnmapped_lowScore = "";
				}
				else
				{
					tmpPeAlignSamStr = "";
					tmpPeAlignInfoStr = "";
					tmpPeAlignSamStr_unpair = peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
						peReadInfo, fasta_or_fastq_bool);
					tmpPeAlignSamStr_bothEndsUnmapped_lowScore = "";					
				}

				peAlignSamVec_fixUnpair[tmpOpenMP] = tmpPeAlignSamStr;
				peAlignInfoVec_fixUnpair[tmpOpenMP] = tmpPeAlignInfoStr;
				peAlignSamVec_unpair_fixUnpair[tmpOpenMP] = tmpPeAlignSamStr_unpair;
				peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmpOpenMP] = tmpPeAlignSamStr_bothEndsUnmapped_lowScore;
				//delete(peReadInfo);
				peAlignInfo->memoryFree();
				delete(peAlignInfo);

			}

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish fixing oneEndUnmapped, turn: " << tmpTurn+1 << endl;// << endl;
			runtime_log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
			
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				if(peAlignSamVec_fixUnpair[tmp] != "")
				{
					OutputSamFile_oneEndMapped_ofs << peAlignSamVec_fixUnpair[tmp] << endl;
					PairedReadNum ++;
				}
				if(peAlignInfoVec_fixUnpair[tmp] != "")
				{
					tmpAlignIncompletePair_ofs << peAlignInfoVec_fixUnpair[tmp] << endl;// << endl;
				}
				if(peAlignSamVec_unpair_fixUnpair[tmp] != "")
				{
					OutputSamFile_oneEndMapped_unpair_ofs << peAlignSamVec_unpair_fixUnpair[tmp] << endl;
					UnpairedReadNum ++;
				}
				if(peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] != "")
				{
					OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs << peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] << endl;
				}
			}	
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish output, turn: " << tmpTurn+1 << endl << endl;// << endl;

		}		
		inputRecord_ifs.close();
	}

	OutputSamFile_oneEndMapped_ofs.close();
	OutputSamFile_oneEndMapped_unpair_ofs.close();
	tmpAlignIncompletePair_ofs.close();
	OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs.close();
	//tmpAlignInfoForDebugFile_oneEndMapped_ofs.close();

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl; 
	
	log_ofs << endl << "**********************************************************************************" << endl << endl;
	runtime_log_ofs << endl << "**********************************************************************************" << endl << endl;

	///////////////////////////////////////  merging incomplete alignment files ///////////////////////////////////////////
	if(outputDirectlyBool_Phase1Only)
	{
		log_ofs << "start to merge incomplete alignment files ..." << endl;
		string cat_cmd_tmpAlignIncompletePair = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " " + tmpAlignIncompletePair + "." + int_to_str(tmp);
		}
		cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " " + tmpAlignIncompletePair;
		cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " > " + tmpAlignIncompletePair + ".all";
		system(cat_cmd_tmpAlignIncompletePair.c_str());
		log_ofs << "finish merging incomplete alignment files ..." << endl;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Sam 2 Junc   ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 

	string juncfile;
	string juncInsFile;
	if(DoSam2JuncBool)
	{
		juncfile = tmpIntermediateJunctionFile;
		juncInsFile = tmpIntermediateJunctionFile + "ins";
		int readLength_max = readLengthMax_tmp; //to Debug
		log_ofs << "... readLengthMax = " << readLength_max << endl;
		int read_width = readLength_max;
		string chrom_dir = chromDirStr;
		int min_intron = 50;
		int max_intron = 300000;
		int min_anchor = 1;
		chrom_dir.append("/");
		double max_rank = 0;

		vector<string> comb_mapreads_files;
		comb_mapreads_files.push_back(tmpAlignCompleteRead);
		comb_mapreads_files.push_back(tmpAlignIncompletePair_SAM);
		comb_mapreads_files.push_back(OutputSamFile_oneEndMapped);

		Covert2JuncComb(juncfile.c_str(), comb_mapreads_files, read_width, chrom_dir, max_rank, min_intron, max_intron, min_anchor);

	}
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On unfixed head/tail Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 

	int junctionNum = 0;
	
	string InputSpliceJunction = tmpIntermediateJunctionFile;
	FILE *fp_spliceJunction = fopen(InputSpliceJunction.c_str(), "r");
	//if(annotation_provided_bool)
	//	InputSpliceJunction = annotation_file_path;

	string OutputSamFile_fixHeadTail_complete_pair = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_pair_ofs(OutputSamFile_fixHeadTail_complete_pair.c_str());
	
	#ifdef CHECK_MULTI	
	string OutputSamFile_fixHeadTail_complete_pair_multi = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam.multi";
	ofstream OutputSamFile_fixHeadTail_complete_pair_mulit_ofs(OutputSamFile_fixHeadTail_complete_pair_multi.c_str());
	string OutputSamFile_fixHeadTail_complete_pair_unique = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam.unique";
	ofstream OutputSamFile_fixHeadTail_complete_pair_unique_ofs(OutputSamFile_fixHeadTail_complete_pair_unique.c_str());
	#endif
	
	string OutputSamFile_fixHeadTail_complete_pair_alignInfo = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam_alignInfo";
	ofstream OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs(OutputSamFile_fixHeadTail_complete_pair_alignInfo.c_str());

	string OutputSamFile_fixHeadTail_incomplete_pair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_ofs(OutputSamFile_fixHeadTail_incomplete_pair.c_str());
	
	#ifdef CHECK_MULTI	
	string OutputSamFile_fixHeadTail_incomplete_pair_multi = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam.multi";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_mulit_ofs(OutputSamFile_fixHeadTail_incomplete_pair_multi.c_str());
	string OutputSamFile_fixHeadTail_incomplete_pair_unique = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam.unique";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_unique_ofs(OutputSamFile_fixHeadTail_incomplete_pair_unique.c_str());
	#endif		
	
	string OutputSamFile_fixHeadTail_incomplete_pair_alignInfo = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam_alignInfo";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs(OutputSamFile_fixHeadTail_incomplete_pair_alignInfo.c_str());	

	string OutputSamFile_fixHeadTail_complete_unpair = outputDirStr + "/phase2_output/fixHeadTail_complete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_unpair_ofs(OutputSamFile_fixHeadTail_complete_unpair.c_str());

	string OutputSamFile_fixHeadTail_incomplete_unpair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_unpair_ofs(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());	

	string OutputSamFile_fixHeadTail_pair_lowScore = outputDirStr + "/phase2_output/fixHeadTail_pair_lowScore.sam";
	ofstream OutputSamFile_fixHeadTail_pair_lowScore_ofs(OutputSamFile_fixHeadTail_pair_lowScore.c_str());	

	SJhash_Info* SJ = new SJhash_Info();
	if(DoRemappingOnUnfixedHeadTailAlignmentBool)
	{
    	log_ofs << "start to build spliceJunction Hash" << endl;

    	bool spliceJunctionHashExists = true;

		string entryString;
		int tabLocation1, tabLocation2, tabLocation3, tabLocation4, tabLocation5;
		char entry[500];
		int chrInt;
		int spliceStartPos;
		int spliceEndPos;
		string chrIntString;
		string spliceStartPosString;
		string spliceEndPosString;

		//SJhash_Info SJ;
		
		/////////////////////////////////////////////////////////////////////////////
		///////////////////// position Hash  ////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		//SJ->initiateSJintHash(indexInfo->chromNum);
		/////////////////////////////////////////////////////////////////////////////
		//////////////////////  string hash /////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());


		fgets(entry, sizeof(entry), fp_spliceJunction);
		while(!feof(fp_spliceJunction))
		{
			fgets(entry, sizeof(entry), fp_spliceJunction);
			if(feof(fp_spliceJunction))
				break;
			junctionNum ++;
			//cout << "entryString: " << entryString << endl;
			entryString = entry;
			tabLocation1 = entryString.find('\t', 0);
			tabLocation2 = entryString.find('\t', tabLocation1+1);
			tabLocation3 = entryString.find('\t', tabLocation2+1);
			chrIntString = entryString.substr(0, tabLocation1);
			spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
			spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
			chrInt = indexInfo->convertStringToInt(chrIntString);
			spliceStartPos = atoi(spliceStartPosString.c_str());
			spliceEndPos = atoi(spliceEndPosString.c_str());

			///////////////////////////////////////////////////////////////////////////
			//////////// position Hash  ///////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////
			//SJ->insert2SJintHash(chrInt, spliceStartPos, spliceEndPos);
			///////////////////////////////////////////////////////////////////////////
			////////////  string hash /////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////			
			SJ->insert2AreaAndStringHash(chrInt, spliceStartPos, spliceEndPos, indexInfo);
		
		}

		fclose(fp_spliceJunction);

		if(junctionNum == 0)
		{
			spliceJunctionHashExists = false;
		}

		log_ofs << "junctionNum = " << junctionNum << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		

		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;

		
		string headTailSoftClippingFile = tmpAlignIncompletePair;// + ".all";

		if(outputDirectlyBool_Phase1Only)
			headTailSoftClippingFile += ".all";

		//FILE *fp_HeadTail = fopen(headTailSoftClippingFile.c_str(), "r");
		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());

		string line1, line2, line3, line4, line5, line6, line7, 
			line8, line9, line10, line11;

		int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;

		bool EndOfRecord = false;

		int tmpTurn = 0;

		int realRecordNum;// = normalRecordNum;

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
		//vector<string> peAlignInfoVec(normalRecordNum);
		vector<string> peAlignSamVec_complete_pair(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_pair(normalRecordNum);

		vector<string> peAlignSamVec_complete_unpair(normalRecordNum);
		vector<string> peAlignSamVec_incomplete_unpair(normalRecordNum);

		vector<string> peAlignSamVec_pair_lowScore(normalRecordNum);	
		//vector<string> peAlignSamVec_incomplete_unpair(normalRecordNum);

		//getline(inputUnfixedHeadTailRecord_ifs, line11);

		for(tmpTurn = 0; ; tmpTurn++)
		{
			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			//cout << "start to read long-head file record" << endl;
			//cout << "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}				
				getline(inputUnfixedHeadTailRecord_ifs, line11);
				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}		
				getline(inputUnfixedHeadTailRecord_ifs, line1);
				getline(inputUnfixedHeadTailRecord_ifs, line2);
				getline(inputUnfixedHeadTailRecord_ifs, line3);
				getline(inputUnfixedHeadTailRecord_ifs, line4);
				getline(inputUnfixedHeadTailRecord_ifs, line5);
				getline(inputUnfixedHeadTailRecord_ifs, line6);
				getline(inputUnfixedHeadTailRecord_ifs, line7);
				getline(inputUnfixedHeadTailRecord_ifs, line8);
				getline(inputUnfixedHeadTailRecord_ifs, line9);
				getline(inputUnfixedHeadTailRecord_ifs, line10);
				//getline(inputUnfixedHeadTailRecord_ifs, line11);

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
		
			//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;

			//cout << "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;

			//cout << "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;

			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				int threadNO = omp_get_thread_num();

				PE_Read_Info peReadInfo;// = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixIncompleteAlignment_getline(
					line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], line3StrVec[tmpOpenMP],
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP], line6StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo, 
					indexInfo, fasta_or_fastq_bool);		

				#ifdef MAP_INFO
				cout << "read_name: " << endl << peReadInfo.returnReadName_1() << endl 
					<< peReadInfo.returnReadName_2() << endl;
				cout << "start fixHeadTail: " << endl;
				cout << "PeAlignInfo:" << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				#endif
				FixHeadTailInfo* fixHeadTailInfo = new FixHeadTailInfo();
				if(Do_fixHeadTail_remapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail);	
				}
				#ifdef MAP_INFO
				cout << "after remapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				#endif
				if(Do_fixHeadTail_greedyMapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo);	
				}
				#ifdef MAP_INFO
				cout << "after greedyMapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				#endif
				if(Do_fixHeadTail_remappingAndTargetMapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail);	
				}
				#ifdef MAP_INFO
				cout << "after remapping And Target Mapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;				
				#endif
				if(Do_fixHeadTail_remappingAgain)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail);	
				}

				fixHeadTailInfo->fixHeadTail_extend2end(peReadInfo, peAlignInfo, indexInfo);
				#ifdef MAP_INFO
				cout << "after extending: " << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;								
				#endif
				// remove duplicate mismatch
				peAlignInfo->removeDuplicateMismatch();
				
				peAlignInfo->chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();							

				string tmpPeAlignSamStr_complete_pair, //tmpPeAlignSamStr_unpair, 
					tmpPeAlignSamStr_incomplete_pair, tmpPeAlignSamStr_complete_unpair, tmpPeAlignSamStr_incomplete_unpair,
					tmpPeAlignSamStr_pair_lowScore;// tmpPeAlignSamStr_incomplete_unpair;

				if(pairExistsBool) // some pair exists, all completed, print out paired SAM info
				{
					bool allFinalPairAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();	
					bool unique_bool = peAlignInfo->checkUniqueOrMulti();
					statsInfo->increPairedNum_fixHeadTail(threadNO, allFinalPairAlignmentCompleteBool, unique_bool);
					
					int alignment_score_min_output 
						= peReadInfo.returnAlignmentScoreMinOutput_withComplement(
							ALIGNMENT_SCORE_MIN_OUTPUT_COMPLEMENT);
					bool align_lowScore_bool 
						= peAlignInfo->alignPairScoreTooLow_bool(alignment_score_min_output); 

					if(align_lowScore_bool)
					{
						statsInfo->increLowScoreComplete_fixHeadTail(threadNO, allFinalPairAlignmentCompleteBool, unique_bool);
						tmpPeAlignSamStr_complete_pair = "";
						tmpPeAlignSamStr_incomplete_pair = "";
						tmpPeAlignSamStr_complete_unpair = "";
						tmpPeAlignSamStr_incomplete_unpair = "";
						tmpPeAlignSamStr_pair_lowScore 
							= peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool);
					}
					else if(allFinalPairAlignmentCompleteBool)
					{
						tmpPeAlignSamStr_complete_pair 
							= peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(
							peReadInfo, fasta_or_fastq_bool);
						tmpPeAlignSamStr_incomplete_pair = "";
						tmpPeAlignSamStr_complete_unpair = "";
						tmpPeAlignSamStr_incomplete_unpair = "";
						tmpPeAlignSamStr_pair_lowScore = "";
					}
					else
					{
						tmpPeAlignSamStr_complete_pair = "";
						tmpPeAlignSamStr_incomplete_pair 
							= peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(
								peReadInfo, fasta_or_fastq_bool);
						tmpPeAlignSamStr_complete_unpair = "";
						tmpPeAlignSamStr_incomplete_unpair = "";	
						tmpPeAlignSamStr_pair_lowScore = "";			
					}
				}
				else //if((!pairExistsBool) && (allAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
				{
					bool allUnpairAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();
					statsInfo->increUnpairedNum_fixHeadTail(threadNO, allUnpairAlignmentCompleteBool);
					if(allUnpairAlignmentCompleteBool)
					{
						tmpPeAlignSamStr_complete_pair = "";
						tmpPeAlignSamStr_incomplete_pair = "";
						tmpPeAlignSamStr_complete_unpair = peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
							peReadInfo, fasta_or_fastq_bool);
						tmpPeAlignSamStr_incomplete_unpair = "";	
						tmpPeAlignSamStr_pair_lowScore = "";
					}
					else
					{
						tmpPeAlignSamStr_complete_pair = "";
						tmpPeAlignSamStr_incomplete_pair = "";						
						tmpPeAlignSamStr_complete_unpair = "";		
						tmpPeAlignSamStr_incomplete_unpair = peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
							peReadInfo, fasta_or_fastq_bool);	
						tmpPeAlignSamStr_pair_lowScore = "";
					}
					//tmpPeAlignSamStr_incomplete_unpair = "";
				}

				peAlignSamVec_complete_pair[tmpOpenMP] = tmpPeAlignSamStr_complete_pair;
				peAlignSamVec_incomplete_pair[tmpOpenMP] = tmpPeAlignSamStr_incomplete_pair;

				peAlignSamVec_complete_unpair[tmpOpenMP] = tmpPeAlignSamStr_complete_unpair;
				peAlignSamVec_incomplete_unpair[tmpOpenMP] = tmpPeAlignSamStr_incomplete_unpair;
				peAlignSamVec_pair_lowScore[tmpOpenMP] = tmpPeAlignSamStr_pair_lowScore;

				//peAlignSamVec_incomplete_unpair[tmpOpenMP] = tmpPeAlignSamStr_incomplete_unpair;

				peAlignInfo->memoryFree();
				//delete(peReadInfo);
				delete peAlignInfo;

			}
			//cout << "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;

			//cout << "start to output ... turn: " << tmpTurn+1 << endl;
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				if(peAlignSamVec_complete_pair[tmp] != "")
				{
					OutputSamFile_fixHeadTail_complete_pair_ofs << peAlignSamVec_complete_pair[tmp] << endl;
					//PairedReadNum ++;
				}

				if(peAlignSamVec_incomplete_pair[tmp] != "")
				{	
					OutputSamFile_fixHeadTail_incomplete_pair_ofs << peAlignSamVec_incomplete_pair[tmp] << endl;
					//PairedReadNum ++;
				}

				if(peAlignSamVec_complete_unpair[tmp] != "")
				{
					OutputSamFile_fixHeadTail_complete_unpair_ofs << peAlignSamVec_complete_unpair[tmp] << endl;
					//UnpairedReadNum ++;
				}
				if(peAlignSamVec_incomplete_unpair[tmp] != "")	
				{
					OutputSamFile_fixHeadTail_incomplete_unpair_ofs << peAlignSamVec_incomplete_unpair[tmp] << endl;
					//UnpairedReadNum ++;
				}
				if(peAlignSamVec_pair_lowScore[tmp] != "")
				{
					OutputSamFile_fixHeadTail_pair_lowScore_ofs << peAlignSamVec_pair_lowScore[tmp] << endl;
				}
			}			
			//cout << "finish output, turn: " << tmpTurn+1 << endl << endl;

		}
		inputUnfixedHeadTailRecord_ifs.close();
	}
	delete(SJ);

	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_pair_lowScore_ofs.close();
	
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);


	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	*/
	statsInfo->getPhase1Stats();
	statsInfo->getFixUnpairedStats();
	statsInfo->getFixHeadTailStats();
	statsInfo->outputAllStats(log_ofs, readTotalNum);
	statsInfo->outputAllStats(stats_ofs, readTotalNum);

	statsInfo->outputFinalStats(stats_ofs, Do_Phase1_Only, readTotalNum);	

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  

	string finalOutputSam = outputDirStr + "/output.sam";
	//if(Do_Phase1_Only)
	//{
		string cat_cmd = "cat "
			+ tmpHeadSectionInfo
			+ " " + tmpAlignCompleteRead  
			//+ " " + tmpAlignOneEndUnmapped
			+ " " + tmpAlignIncompletePair_SAM 
			+ " " + tmpAlignOneEndUnmapped
			+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile
			+ " " + tmpAlignBothEndsUnmapped_lowScore
			+ " " + tmpAlignBothEndsUnmapped
			+ " > " + finalOutputSam;
		system(cat_cmd.c_str()); 
	/*}
	else
	{
		string cat_cmd = "cat " 
			+ tmpHeadSectionInfo
			+ " " + tmpAlignCompleteRead 
			+ " " + OutputSamFile_oneEndMapped 
			+ " " + OutputSamFile_fixHeadTail_complete_pair
			+ " " + OutputSamFile_fixHeadTail_incomplete_pair
			+ " " + OutputSamFile_fixHeadTail_complete_unpair
			+ " " + OutputSamFile_fixHeadTail_incomplete_unpair
			+ " " + OutputSamFile_oneEndMapped_unpair
			+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile
			+ " " + tmpAlignBothEndsUnmapped_lowScore
			+ " " + OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore
			+ " " + OutputSamFile_fixHeadTail_pair_lowScore
			+ " " + tmpAlignBothEndsUnmapped 
			+ " > " + finalOutputSam;
		system(cat_cmd.c_str()); 			
	}*/


	if(removeAllIntermediateFilesBool)
	{
		// if(!Do_Phase1_Only)
		// {	
		// 	remove(InputSpliceJunction.c_str());
		// 	remove(juncInsFile.c_str());
			remove(tmpAlignIncompletePair.c_str());
			// remove(OutputSamFile_fixHeadTail_complete_pair.c_str());
			// remove(OutputSamFile_fixHeadTail_incomplete_pair.c_str());
			// remove(OutputSamFile_fixHeadTail_complete_unpair.c_str());
			// remove(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());		
		
			//remove(OutputSamFile_oneEndMapped.c_str());
			//remove(OutputSamFile_oneEndMapped_unpair.c_str());
		//}
		remove(tmpAlignOneEndUnmapped.c_str());
		remove(tmpAlignIncompletePair_SAM.c_str());
		remove(tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
		remove(tmpAlignBothEndsUnmapped.c_str());
		remove(tmpAlignIncompletePair.c_str());
		remove(tmpAlignCompleteRead_alignInfo.c_str());
		//remove(OutputSamFile_oneEndMapped_alignInfo.c_str());
		remove(tmpAlignCompleteRead.c_str());
		remove(tmpHeadSectionInfo.c_str());

		remove(tmpAlignBothEndsUnmapped_lowScore.c_str());
		//remove(OutputSamFile_fixHeadTail_pair_lowScore.c_str());
		//remove(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());
	}

	annotation_ifs.close();
	delete annotationInfo;	
	delete indexInfo;

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  

    return 0;
} //end main
