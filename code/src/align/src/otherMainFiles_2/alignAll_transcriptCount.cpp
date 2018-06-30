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
#include "general/extractUnmapAlignment2ReadFile.h"
#include "phase1/arrayQueue.h"
#include "stats_info.h"
#include "constantDefinitions.h"
#include "general/option_info.h"
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
//#include "segmentMapping_secondLevel.h"
#include "general/splice_info.h"
#include "general/fixGapRelationParameters.h"
#include "general/read_info.h"
#include "general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "general/transcript_count_vec.h"
#include "general/fixDoubleAnchorNWDP_info.h"
#include "general/fixDoubleAnchorMatch_info.h"
#include "general/fixDoubleAnchorInsertion_info.h"
#include "general/fixDoubleAnchorDeletion_info.h"
#include "general/fixDoubleAnchorSplice_complicate_info.h"
#include "general/fixDoubleAnchorSplice_info.h"
#include "general/fixDoubleAnchorCirRNA_info.h"
#include "general/path_info.h"
#include "general/gap_info.h"
#include "general/align_info.h"
#include "general/peAlign_info.h"
#include "general/groupSeg_info.h"
#include "general/alignInferJunctionHash_info.h"
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
#include "general/otherFunc2.h"
#include "general/alignmentToJunc.h"

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

	int normalRecordNum_1stMapping = 500000;//2000000;
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
	log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	log_ofs << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;
	log_ofs << "finish loading chromosomes" << endl;

	/////////////////////////////////////   start to load real genome index  /////////////////////////////////////
	
	string genome_chrom_bit_file = realGenome_file_path; genome_chrom_bit_file += "/_chrom"; 
	ifstream genome_chrom_bit_file_ifs(genome_chrom_bit_file.c_str());
	string genome_parameter_file = realGenome_file_path; genome_parameter_file += "/_parameter"; 
	ifstream genome_parameter_file_ifs(genome_parameter_file.c_str());

	Index_Info*  genomeIndexInfo = new Index_Info(genome_parameter_file_ifs, log_ofs);
	log_ofs << "real genome index: " << realGenome_file_path << endl;
	cout << "start to load whole real genome " << endl;
	log_ofs << "start to load whole real genome " << endl;
	char *genome_chrom; genome_chrom = (char*)malloc((genomeIndexInfo->returnIndexSize()) * sizeof(char));
	genome_chrom_bit_file_ifs.read((char*)genome_chrom, (genomeIndexInfo->returnIndexSize()) * sizeof(char));
	genomeIndexInfo->readGenome(genome_chrom);
	log_ofs << "real genome size = " << genomeIndexInfo->returnChromStringLength() << endl;
	cout << "start to initiate and load every chromosome for ref genome " << endl;
	log_ofs << "start to initiate and load every chromosome for ref genome " << endl;
	genomeIndexInfo->initiate();
	log_ofs << "start to initiate real genome chrNameIndexArray" << endl;
	genomeIndexInfo->initiateChrNameIndexArray(1000);	
	log_ofs << "finish loading real genome" << endl;

	////////////////// transcript info loading .... //////////
	log_ofs << "start to load transcript info" << endl;
	cout << "start to load transcriptInfo " << endl;
	settings_log_ofs << "transcriptInfo: " << optionInfo->transcript_file_path << endl;
	cout << "transcriptInfo: " << optionInfo->transcript_file_path << endl;
	Transcript_Set* transcriptInfo = new Transcript_Set();
	string transcript_file_path = optionInfo->transcript_file_path;
	ifstream transcript_file_ifs(transcript_file_path.c_str());
	string transcript_type = optionInfo->transcript_type;
	transcriptInfo->extractTranscript(transcript_file_ifs, genomeIndexInfo, transcript_type);

	////////////////// initiating transcriptCount_merged info ..... ///////////
	log_ofs << "start to initiating transcriptCountInfo_merged ..... " << endl;
	cout << "start to initiating transcriptCountInfo_merged ..... " << endl;
	Transcript_Count* transcriptCountInfo_merged = new Transcript_Count();
	transcriptCountInfo_merged->initiate_withTranscriptSetInfo(transcriptInfo);
	////////////////// initiating transcriptCount_vec info .... ///////////////
	log_ofs << "start to initiating transcriptCountInfo vec  ..... " << endl;
	cout << "start to initiating transcriptCountInfo vec ..... " << endl;
	Transcript_Count_Vec* transcriptCountInfo_vec = new Transcript_Count_Vec();
	transcriptCountInfo_vec->initiate(threads_num, transcriptInfo);

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

	string transcriptCount_file_str = outputDirStr + "/transcriptCount.txt";
	ofstream transcriptCount_ofs(transcriptCount_file_str.c_str());
	string transcriptCount_otherStats_file_str = outputDirStr + "/transcriptCount_otherStats.txt";
	ofstream transcriptCount_otherStats_ofs(transcriptCount_otherStats_file_str.c_str());

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
			//#ifdef MAP_INFO
			//cout << "transcript_peAlignInfo: \n" << peAlignInfo_transcript.returnPeAlignInfoStr() << endl << endl;
			//cout << "finish outputing transcript_peAlignInfo" << endl;
			//#endif			
			// PE_Read_Alignment_Info peAlignInfo;
			// peAlignInfo.convertTranscriptAlignInfo2GenomeAlignInfo(peAlignInfo_transcript,
			// 	transcriptInfo, genomeIndexInfo);
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

			//cout << "start to choose the best alignment" << endl;
			peAlignInfo_transcript.removeDuplicateMismatch(false);
			//peAlignInfo.removeDuplicateMismatch(false);
			
			peAlignInfo_transcript.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
			//cout << "end of choosing the best alignment" << endl;

			//peAlignInfo.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();

			// if(!outputDirectlyBool_Phase1Only)
			// {	
			// 	PeAlignSamStrVec_complete[tmpOpenMP] = "";			
			// 	PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = "";
			// 	PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = "";
			// 	PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = "";
			// 	PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmpOpenMP] = "";
			// 	PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmpOpenMP] = "";
			// 	PeAlignSamStrVec_inCompletePair[tmpOpenMP] = "";
			// }
			#ifdef CAL_TIME
			selectBestAlign_end = clock();
			selectBestAlign_cost = selectBestAlign_cost + selectBestAlign_end - selectBestAlign_begin;	
			#endif		
			
			#ifdef CAL_TIME
			getSamFormat_begin = clock();
			#endif
			//cout << "start to do transcriptCount " << endl;

			peAlignInfo_transcript.doTranscriptCount(transcriptCountInfo_vec,
				threadNO, transcriptInfo);

			//cout << "end of doing transcriptCount " << endl;
			// peAlignInfo_transcript.output_phase1(
			// 	outputDirectlyBool_Phase1Only,
			// 	Do_Phase1_Only,
			// 	PeAlignSamStrVec_complete,
			// 	PeAlignInfoStrVec_inCompletePair,
			// 	PeAlignInfoStrVec_oneEndUnmapped,
			// 	PeAlignSamStrVec_bothEndsUnmapped,
			// 	PeAlignSamStrVec_bothEndsUnmapped_lowScore,
			// 	PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
			// 	PeAlignSamStrVec_inCompletePair,
			// 	//PeAlignInfoStrVec_completePaired,
			// 	repeatRegionInfoVec,
			// 	readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo,
			// 	fasta_or_fastq_bool);
			// peAlignInfo.output_phase1(
			// 	outputDirectlyBool_Phase1Only,
			// 	Do_Phase1_Only,
			// 	PeAlignSamStrVec_complete,
			// 	PeAlignInfoStrVec_inCompletePair,
			// 	PeAlignInfoStrVec_oneEndUnmapped,
			// 	PeAlignSamStrVec_bothEndsUnmapped,
			// 	PeAlignSamStrVec_bothEndsUnmapped_lowScore,
			// 	PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
			// 	PeAlignSamStrVec_inCompletePair,
			// 	//PeAlignInfoStrVec_completePaired,
			// 	repeatRegionInfoVec,
			// 	readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo,
			// 	fasta_or_fastq_bool);

			#ifdef CAL_TIME
			getSamFormat_end = clock();
			getSamFormat_cost = getSamFormat_cost + getSamFormat_end - getSamFormat_begin;

			freeMem_begin = clock();
			#endif

			fixPhase1Info.memoryFree();
			//delete fixPhase1Info;
			//delete readInfo; 
			peAlignInfo_transcript.memoryFree();
			//peAlignInfo.memoryFree();
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
		
		// if(!outputDirectlyBool_Phase1Only)
		// {
		// 	for(int tmp = 0; tmp < realRecordNum; tmp++)
		// 	{	
		// 		if(PeAlignSamStrVec_complete[tmp] != "")
		// 		{
		// 			tmpAlignCompleteRead_ofs << PeAlignSamStrVec_complete[tmp] << endl;
		// 			PairedReadNum ++;
		// 			//PeAlignSamStrVec_complete[tmp] = "";
		// 		}			

		// 		if(PeAlignInfoStrVec_inCompletePair[tmp] != "")
		// 		{
		// 			tmpAlignIncompletePair_ofs << PeAlignInfoStrVec_inCompletePair[tmp] << endl;
		// 			//PeAlignInfoStrVec_inCompletePair[tmp] = "";
		// 		}
				
		// 		if(PeAlignInfoStrVec_oneEndUnmapped[tmp] != "")
		// 		{
		// 			tmpAlignOneEndUnmapped_ofs << PeAlignInfoStrVec_oneEndUnmapped[tmp] << endl;
		// 			//PeAlignInfoStrVec_oneEndUnmapped[tmp] = "";
		// 		}
				
		// 		if(PeAlignSamStrVec_bothEndsUnmapped[tmp] != "")
		// 		{
		// 			tmpAlignBothEndsUnmapped_ofs << PeAlignSamStrVec_bothEndsUnmapped[tmp] << endl;
		// 			//PeAlignSamStrVec_bothEndsUnmapped[tmp] = "";
		// 			BothUnmappedReadNum ++;
		// 		}
		// 		if(PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] != "")
		// 		{
		// 			tmpAlignBothEndsUnmapped_lowScore_ofs << PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] << endl;
		// 		}
		// 		if(PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] != "")
		// 		{
		// 			tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs << PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] << endl;
		// 			BothUnmappedReadNum_mappedToRepeatRegion ++;
		// 		}

		// 		if(PeAlignSamStrVec_inCompletePair[tmp] != "")
		// 		{
		// 			tmpAlignIncompletePair_SAM_ofs << PeAlignSamStrVec_inCompletePair[tmp] << endl;
		// 		}

		// 	}
		// }
		#ifdef CAL_TIME
		output_end = clock();
		output_cost = output_cost + output_end - output_begin;
		#endif
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) 
			<< "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
	}
	// repeatRegionFile_ofs << "Repeat Region Info: size = " << repeatRegionInfoVec.size() << endl;
	// for(int tmpThread = 0; tmpThread < threads_num; tmpThread++)
	// {
	// 	repeatRegionInfoVec[tmpThread]->outputRepeatRegion(tmpThread+1, indexInfo, sa, 100, repeatRegionFile_ofs);

	// }

	cout << "readTotalNum: " << readTotalNum << endl;
	settings_log_ofs << "readTotalNum: " << readTotalNum << endl;

	transcriptCountInfo_vec->merge2oneTranscriptCount(transcriptCountInfo_merged);
	transcriptCountInfo_merged->output(
		transcriptInfo, transcriptCount_ofs, transcriptCount_otherStats_ofs);

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

	free(preIndexMapLengthArray); free(preIndexIntervalStartArray); free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);
	free(verifyChild);

	free(chrom);
	
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
	if(threads_num != 1)
		log_ofs << "!!! set thread=1 to get the time-cost information ... !!!" << endl;

	#endif

	log_ofs << endl << "**********************************" << endl << "**********************************" << endl;

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


	if(removeAllIntermediateFilesBool)
	{

		remove(tmpAlignIncompletePair.c_str());
		remove(tmpAlignOneEndUnmapped.c_str());
		remove(tmpAlignIncompletePair_SAM.c_str());
		remove(tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
		remove(tmpAlignBothEndsUnmapped.c_str());
		remove(tmpAlignIncompletePair.c_str());
		remove(tmpAlignCompleteRead_alignInfo.c_str());
		remove(tmpAlignCompleteRead.c_str());
		remove(tmpHeadSectionInfo.c_str());
		remove(tmpAlignBothEndsUnmapped_lowScore.c_str());
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
