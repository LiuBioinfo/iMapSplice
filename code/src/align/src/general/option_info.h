// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef OPTION_INFO_H
#define OPTION_INFO_H

#include <string.h>
#include <string>

#include <getopt.h>

using namespace std;

class Option_Info
{
public:
	string read_file_path_1;
	//bool file_1_fasta_or_fastq_bool;
	string read_file_path_2;
	//bool file_2_fasta_or_fastq_bool;
	string read_file_path_SE;

	int threads_num;
	string global_index_file_path_prefix;
	string local_index_file_path_prefix;
	string chromsome_file_path_prefix;
	bool fasta_or_fastq_bool;
	string outputFolder_path;
	bool Do_phase1_only_bool;

	bool SE_or_PE_bool;

	int mappedLength_perc_min;

	int mappedLength_base_min;

	string annotation_file_path;
	bool annotation_provided_bool;

	string spliceJunctionAlignInferHash_file_path;
	bool spliceJunctionAlignInferHash_provided_bool;

	string realGenome_file_path;
	bool realGenome_provided_bool;

	string transcript_file_path;
	bool transcript_provided_bool;

	string transcript_type;
	bool transcript_type_assigned_bool;
	string optionStr;

	string inputSamFilePath;
	string alreadyMappedAlignmentFile;
	bool extractUnmapAlignment2ReadFile_bool;

	string SNPfilePath;
	bool SNP_provided_bool;

	string SNP_seq_index_path;
	bool SNP_seq_index_provided_bool;
	Option_Info()
	{
		optionStr = "W:1:2:S:G:L:C:T:F:O:M:I:Z:P:Q:R:U:Y:AH";	
		SE_or_PE_bool = false;
		Do_phase1_only_bool = false;
		annotation_provided_bool = false;
		spliceJunctionAlignInferHash_provided_bool = false;
		realGenome_provided_bool = false;
		transcript_provided_bool = false;
		transcript_type_assigned_bool = false;
		extractUnmapAlignment2ReadFile_bool = false;
		mappedLength_perc_min = 0;
		mappedLength_base_min = 0;
		SNP_provided_bool = false;
		string transcript_type = "NULL";
		SNP_seq_index_provided_bool = false;
	}

	string returnInputSamFilePath()
	{
		return inputSamFilePath;
	}

	bool return_SNP_provided_bool()
	{
		return SNP_provided_bool;
	}

	bool return_SNP_seq_index_provided_bool()
	{
		return SNP_seq_index_provided_bool;
	}

	bool return_extractUnmapAlignment2ReadFile_bool()
	{
		return extractUnmapAlignment2ReadFile_bool;
	}

	bool returnSEorPE_bool()
	{
		// fix me
		return SE_or_PE_bool;
	}

	void getOpt_long(int argc, char**argv)
	{
		struct option long_options[] {
			{"SNPfile", 1, NULL, 'P'},
			{"SNP_seq_index_path", 1, NULL, 'Q'},			
			{"read_singleEnd_path", 1, NULL, 'S'},
			{"input_SAM_file_path", 1, NULL, 'W'},
			{"read_end1_path", 1, NULL, '1'},
			{"read_end2_path", 1, NULL, '2'},
			{"global_index_path", 1, NULL, 'G'},
			{"local_index_path", 1, NULL, 'L'},
			{"ref_chr_path", 1, NULL, 'C'},
			{"threads", 1, NULL, 'T'},
			{"read_format", 1, NULL, 'F'},
			{"output_path", 1, NULL, 'O'},
			{"min_mapped_perc", 1, NULL, 'M'},
			{"annotation_path", 1, NULL, 'I'},
			{"spliceJunctionAlignInferHash_path", 1, NULL, 'Z'},
			{"real_genome_path", 1, NULL, 'R'},
			{"transcript_path", 1, NULL, 'U'},
			{"Do_phase1_only", 0, NULL, 'A'},
			{"help", 0, NULL, 'H'},
			{"transcript_type", 1, NULL, 'Y'},
			{NULL, 0, NULL, 0},
		};		
		char ch;
		string threadsStr, formatStr, min_mapped_perc_str, optargStr;

		bool option_S_set = false;
		bool option_1_set = false;
		bool option_2_set = false;
		bool option_O_set = false;
		bool option_T_set = false;
		bool option_F_set = false;
		bool option_G_set = false;
		bool option_L_set = false;
		bool option_C_set = false;

		cout << "command line: " << endl;
		while( (ch = getopt_long(argc, argv, optionStr.c_str(), long_options, NULL)) != -1 )
		{
		    //printf("optind:%d\n",optind);
		    //printf("optarg:%s\n",optarg);
		    //printf("ch:%c\n",ch);
		    switch(ch)
		    {
		    	case 'Y':
		    		transcript_type = optarg;
		    		if((transcript_type == "GAF")||(transcript_type == "BEER"))
		    		{
		    			transcript_type_assigned_bool = true;
		    		}
		    		else
		    		{
		    			cout << "transcript_type invalid! Must be GAF or BEER" << endl;
		    			exit(1);
		    		}
		    		break;
		    	case 'Q':
		    		SNP_seq_index_path = optarg;
		    		SNP_seq_index_provided_bool = true;
		    		break;		    		
		    	case 'P':
		    		SNPfilePath = optarg;
		    		SNP_provided_bool = true;
		    		break;
		    	case 'U':
		    		transcript_file_path = optarg;
		    		transcript_provided_bool = true;
		    		break;
		    	case 'R':
		    		realGenome_file_path = optarg;
		    		realGenome_provided_bool = true;
		    		break;
		    	case 'I':
		    		annotation_file_path = optarg;
		    		annotation_provided_bool = true;
		    		break;
		    	case 'Z':
		    		spliceJunctionAlignInferHash_file_path = optarg;
		    		spliceJunctionAlignInferHash_provided_bool = true;
		    		break;
		    	case 'M':
			        min_mapped_perc_str = optarg;
			        mappedLength_perc_min = atoi(min_mapped_perc_str.c_str());		    		
		    		break;
		    	case 'H':
		    		cout << this->optionInfoHelpStr() << endl;
		    		exit(1);
		    		break;
		      	case 'A':
		      		Do_phase1_only_bool = true;
		      		break;
		      	case 'S':
		      		read_file_path_SE = optarg;
		      		option_S_set = true;
		      		break;
		      	case 'W':
		      		optargStr = optarg;
			        read_file_path_1 = optargStr + ".unmapped.read.1.fq";
			        option_1_set = true;
		        	read_file_path_2 = optargStr + ".unmapped.read.2.fq";
			        option_2_set = true;
			        alreadyMappedAlignmentFile = optargStr + ".mapped.alignment.sam";
			        extractUnmapAlignment2ReadFile_bool = true; 
			        inputSamFilePath = optargStr;  	      	
		      		break;
			    case '1':
			       //printf("option 1:'%s'\n",optarg);
			       read_file_path_1 = optarg;
			       option_1_set = true;
			       break;
			    case '2':
		         	//printf("option 2:'%s'\n",optarg);
		        	read_file_path_2 = optarg;
			        option_2_set = true;
			        break;
			    case 'G':
			        //printf("option G:'%s'\n",optarg);
			        global_index_file_path_prefix = optarg;
			        option_G_set = true;
			        break;
			    case 'L':
			        //printf("option L:'%s'\n",optarg);
			        local_index_file_path_prefix = optarg;
			        option_L_set = true;
			        break;
			    case 'C':
			        //printf("option C:'%s'\n",optarg);
			        chromsome_file_path_prefix = optarg;
			        option_C_set = true;
			        break;
			    case 'T':
			        //printf("option T:'%s'\n",optarg);
			        threadsStr = optarg;
			        threads_num = atoi(threadsStr.c_str());
			        option_T_set = true;
			        break;
			    case 'O':
			        //printf("option O:'%s'\n",optarg);
			        outputFolder_path = optarg;
			        option_O_set = true;
			        break;
			    case 'F':
			        //printf("option F:'%s'\n",optarg);
			        formatStr = optarg;
			        if((formatStr == "fasta")||(formatStr == "Fasta")||(formatStr == "FASTA"))
			        {
			        	fasta_or_fastq_bool = true;
			        }
			        else if((formatStr == "fastq")||(formatStr == "Fastq")||(formatStr == "FASTQ"))
			        {
			        	fasta_or_fastq_bool = false;
			        }
			        else
			        {
			        	cout << "input format error! Please claim '-F fasta' or '-F fastq' !" << endl;
			        	exit(1);
			        }
			        option_F_set = true;
			        break;
			    default:
			        printf("other option:%c\n",ch);
		    }
		    //printf("optopt+%c\n",optopt);
		}

		if(option_1_set && option_2_set && (!option_S_set))
		{	
			fasta_or_fastq_bool = this->checkInputFileFormat(read_file_path_1, read_file_path_2);
			SE_or_PE_bool = false;
		}
		else if((!option_1_set) && (!option_2_set) && option_S_set)
		{
			fasta_or_fastq_bool = this->checkInputFileFormat(read_file_path_SE);
			SE_or_PE_bool = true;
		}
		else
		{
			cout << "invalid read_path set" << endl;
			cout << "option 1 & option 2 or option S should be set" << endl;			
		}

		if(!option_G_set)
		{
			cout << "option G unset" << endl;
			//exit(1);
		}
		if(!option_L_set)
		{
			cout << "option L unset" << endl;
			//exit(1);
		}
		// if(!option_C_set)
		// {
		// 	cout << "option C unset" << endl;
		// 	//exit(1);
		// }
		if(!option_T_set)
		{
			cout << "option T unset" << endl;
			//exit(1);
		}
		if(!option_O_set)
		{
			cout << "option O unset" << endl;
			//exit(1);
		}
		if( (!((option_1_set && option_2_set && (!option_S_set))||((!option_1_set) && (!option_2_set) && option_S_set)))
			||(!option_O_set)||(!option_T_set)||
			//(!option_F_set)||
			(!option_G_set)||(!option_L_set)//||(!option_C_set)
			//(!option_1_set)||(!option_1_set)|| 
			)
		{
			exit(1);
		}

		if(transcript_type_assigned_bool && transcript_provided_bool && realGenome_provided_bool)
		{}
		else if((!transcript_type_assigned_bool)&&(!transcript_provided_bool)&&(!realGenome_provided_bool))
		{}
		else
		{
			cout << "transcriptome prebuilt mode failed! " << endl;
			cout << "transcript_provided_bool: " << transcript_provided_bool << endl;
			cout << "realGenome_provided_bool: " << realGenome_provided_bool << endl;
			cout << "transcript_type_assigned_bool: " << transcript_type_assigned_bool << endl;
			cout << "all the above options should be set 0 or 1 !" << endl;
			exit(1);
		}
	}

	/*
	void getOpt(int argc, char**argv)
	{
		char ch;
		string threadsStr, formatStr, min_mapped_perc_str;

		bool option_1_set = false;
		bool option_2_set = false;
		bool option_O_set = false;
		bool option_T_set = false;
		bool option_F_set = false;
		bool option_G_set = false;
		bool option_L_set = false;
		bool option_C_set = false;

		cout << "command line: " << endl;
		while( (ch = getopt(argc, argv, optionStr.c_str())) != -1 )
		{
		    //printf("optind:%d\n",optind);
		    //printf("optarg:%s\n",optarg);
		    //printf("ch:%c\n",ch);
		    switch(ch)
		    {
		    	case 'U':
		    		transcript_file_path = optarg;
		    		transcript_provided_bool = true;
		    		break;
		    	case 'R':
		    		realGenome_file_path = optarg;
		    		realGenome_provided_bool = true;
		    		break;
		    	case 'I':
		    		annotation_file_path = optarg;
		    		annotation_provided_bool = true;
		    		break;
		    	case 'M':
			        min_mapped_perc_str = optarg;
			        mappedLength_perc_min = atoi(min_mapped_perc_str.c_str());		    		
		    		break;
		    	case 'H':
		    		cout << this->optionInfoHelpStr() << endl;
		    		exit(1);
		    		break;
		      	case 'A':
		      		Do_phase1_only_bool = true;
		      		break;
			    case '1':
			       //printf("option 1:'%s'\n",optarg);
			       read_file_path_1 = optarg;
			       option_1_set = true;
			       break;
			    case '2':
		         	//printf("option 2:'%s'\n",optarg);
		        	read_file_path_2 = optarg;
			        option_2_set = true;
			        break;
			    case 'G':
			        //printf("option G:'%s'\n",optarg);
			        global_index_file_path_prefix = optarg;
			        option_G_set = true;
			        break;
			    case 'L':
			        //printf("option L:'%s'\n",optarg);
			        local_index_file_path_prefix = optarg;
			        option_L_set = true;
			        break;
			    case 'C':
			        //printf("option C:'%s'\n",optarg);
			        chromsome_file_path_prefix = optarg;
			        option_C_set = true;
			        break;
			    case 'T':
			        //printf("option T:'%s'\n",optarg);
			        threadsStr = optarg;
			        threads_num = atoi(threadsStr.c_str());
			        option_T_set = true;
			        break;
			    case 'O':
			        //printf("option O:'%s'\n",optarg);
			        outputFolder_path = optarg;
			        option_O_set = true;
			        break;
			    case 'F':
			        //printf("option F:'%s'\n",optarg);
			        formatStr = optarg;
			        if((formatStr == "fasta")||(formatStr == "Fasta")||(formatStr == "FASTA"))
			        {
			        	fasta_or_fastq_bool = true;
			        }
			        else if((formatStr == "fastq")||(formatStr == "Fastq")||(formatStr == "FASTQ"))
			        {
			        	fasta_or_fastq_bool = false;
			        }
			        else
			        {
			        	cout << "input format error! Please claim '-F fasta' or '-F fastq' !" << endl;
			        	exit(1);
			        }
			        option_F_set = true;
			        break;
			    default:
			        printf("other option:%c\n",ch);
		    }
		    //printf("optopt+%c\n",optopt);
		}

		if(!option_1_set)
		{
			cout << "option 1 unset" << endl;
			//exit(1);
		}
		if(!option_2_set)
		{
			cout << "option 2 unset" << endl;
			//exit(1);
		}
		if(!option_G_set)
		{
			cout << "option G unset" << endl;
			//exit(1);
		}
		if(!option_L_set)
		{
			cout << "option L unset" << endl;
			//exit(1);
		}
		if(!option_C_set)
		{
			cout << "option C unset" << endl;
			//exit(1);
		}
		if(!option_T_set)
		{
			cout << "option T unset" << endl;
			//exit(1);
		}
		if(!option_O_set)
		{
			cout << "option O unset" << endl;
			//exit(1);
		}
		if( (!option_1_set)||(!option_2_set)||(!option_O_set)||(!option_T_set)||
			//(!option_F_set)
			(!option_G_set)||(!option_L_set)||(!option_C_set)
			//(!option_1_set)||(!option_1_set)|| 
			)
		{
			exit(1);
		}

	}*/

	bool checkInputFileFormat(const string& read_file_path_1)
	{
		int read_file_path_1_len = read_file_path_1.length();
		bool file_1_fasta_or_fastq_bool;
		// check input file name;
		if(read_file_path_1_len > 6)
		{
			string tmpSuffix_len3 = read_file_path_1.substr(read_file_path_1_len-3);
			string tmpSuffix_len6 = read_file_path_1.substr(read_file_path_1_len-6);
			if((tmpSuffix_len3 == ".fa")||(tmpSuffix_len3 == ".FA")||(tmpSuffix_len3 == ".Fa")
				||(tmpSuffix_len6 == ".fasta")||(tmpSuffix_len6 == ".FASTA")||(tmpSuffix_len6 == ".Fasta"))
			{
				file_1_fasta_or_fastq_bool = true;
			}
			else if((tmpSuffix_len3 == ".fq")||(tmpSuffix_len3 == ".FQ")||(tmpSuffix_len3 == ".Fq")
				||(tmpSuffix_len6 == ".fastq")||(tmpSuffix_len6 == ".FASTQ")||(tmpSuffix_len6 == ".Fastq"))
			{
				file_1_fasta_or_fastq_bool = false;
			}
			else
			{
				cout << read_file_path_1 << " is not a regular fasta or fastq format file name !\n Please use -H to see detailed information or use -F to specify it." << endl;
				exit(1);
			}
		}
		else if(read_file_path_1_len > 3)
		{
			string tmpSuffix_len3 = read_file_path_1.substr(read_file_path_1_len-3);
			if((tmpSuffix_len3 == ".fa")||(tmpSuffix_len3 == ".FA")||(tmpSuffix_len3 == ".Fa"))
			{
				file_1_fasta_or_fastq_bool = true;
			}
			else if((tmpSuffix_len3 == ".fq")||(tmpSuffix_len3 == ".FQ")||(tmpSuffix_len3 == ".Fq"))
			{
				file_1_fasta_or_fastq_bool = false;
			}
			else
			{
				cout << read_file_path_1 << " is not a regular fasta or fastq format file name !\n Please use -H to see detailed information or use -F to specify it." << endl;
				exit(1);
			}
		}
		else
		{
			cout << read_file_path_1 << " is not a regular fasta or fastq format file name !\n Please use -H to see detailed information or use -F to specify it." << endl;
			exit(1);			
		}
		return file_1_fasta_or_fastq_bool;				
	}

	bool checkInputFileFormat(const string& read_file_path_1, const string& read_file_path_2) // true -- fasta, false -- fastq
	{
		bool file_1_format_bool = this->checkInputFileFormat(read_file_path_1);
		bool file_2_format_bool = this->checkInputFileFormat(read_file_path_2);
		if(file_1_format_bool != file_2_format_bool)
		{
			cout << read_file_path_1 << " and " << read_file_path_2 << " are in different formats, please check it." << endl;
			exit(1);
		}
		else
		{
			return file_1_format_bool;
		}
	}

	string optionInfoHelpStr()
	{
		string optionStr = "Processing Paired-end reads:\ncommand-line main arguments:\n-1 <string> path_to_read_file_end1\n-2 <string> path_to_read_file_end2\n-G <string> path_to_global_index\n-L <string> path_to_local_index\n-C <string> path_to_chromosome_files_folder\n-T <int> threads_num\n-O <string> path_to_output_folder\noptional arguments:\n-A Do_phase1_only\n-I annotation file\n-F <string> input_format_fasta_or_fastq\n-H command line help information\n";
		optionStr = optionStr + "Processing Single-end reads:\ncommand-line main arguments:\n-S <string> path_to_read_file_SE\n-G <string> path_to_global_index\n-L <string> path_to_local_index\n-C <string> path_to_chromosome_files_folder\n-T <int> threads_num\n-O <string> path_to_output_folder\noptional arguments:\n-A Do_phase1_only\n-I annotation file\n-F <string> input_format_fasta_or_fastq\n-H command line help information\n";
		
		return optionStr;
	}

	void outputOptStr(ofstream& output_ofs)
	{
		output_ofs << "read_file_path_1: " << read_file_path_1 << "\nread_file_path_2: " << read_file_path_2
			<< "\nread_file_path_SE: " << read_file_path_SE 
			<< "\nthreads_num: " << threads_num << "\ninput_format_fasta_or_fastq: " << fasta_or_fastq_bool 
			<< "\nglobal_index_file_path_prefix: " << global_index_file_path_prefix << "\nlocal_index_file_path_prefix: "
			<< local_index_file_path_prefix << "\nchromsome_file_path_prefix: " << chromsome_file_path_prefix 
			<< "\noutputFolder_path: " << outputFolder_path << "\nmin_mapped_perc: " << mappedLength_perc_min 
			<< "\n\nannotation_provided_bool: " << annotation_provided_bool << "\nannotation_file_path: " 
			<< annotation_file_path << endl << "\nspliceJunctionAlignInferHash_provided_bool: "
			<< spliceJunctionAlignInferHash_provided_bool << "\nspliceJunctionAlignInferHash_file_path: "
			<< spliceJunctionAlignInferHash_file_path << endl;
	}

	void outputSwitchInfo(bool Do_Phase1_Only, bool outputAlignInfoAndSamForAllPairedAlignmentBool,
		bool removeAllIntermediateFilesBool, bool Do_cirRNA, bool outputDirectlyBool_Phase1Only, 
		int normalRecordNum_1stMapping, int normalRecordNum_fixOneEndUnmapped,
		int normalRecordNum_fixHeadTail, bool Do_extendHeadTail_phase1, 
		bool Do_extendHeadTail_fixOneEndUnmapped, bool Do_extendHeadTail_fixHeadTail, 
		bool Do_fixHeadTail_remapping, bool Do_fixHeadTail_greedyMapping,
		bool Do_fixHeadTail_remappingAndTargetMapping, bool Do_fixHeadTail_remappingAgain,
		ofstream& log_ofs)
	{
		log_ofs << endl << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
		log_ofs << endl << "Do_Phase1_Only: " << Do_Phase1_Only << endl << "outputAlignInfoAndSamForAllPairedAlignmentBool: " 
			<< outputAlignInfoAndSamForAllPairedAlignmentBool << endl << "removeAllIntermediateFilesBool: " 
			<< removeAllIntermediateFilesBool << endl << "Do_cirRNA: " << Do_cirRNA << endl
			<< "outputDirectlyBool_Phase1Only: " << outputDirectlyBool_Phase1Only << endl; ;
		log_ofs << endl << "normalRecordNum_1stMapping: " << normalRecordNum_1stMapping << endl;
		log_ofs << "normalRecordNum_fixOneEndUnmapped: " << normalRecordNum_fixOneEndUnmapped << endl;
		log_ofs << "normalRecordNum_fixHeadTail: " << normalRecordNum_fixHeadTail << endl;

		log_ofs << endl << "Do_extendHeadTail_phase1: " << Do_extendHeadTail_phase1 << endl;
		log_ofs << endl << "Do_extendHeadTail_fixOneEndUnmapped: " << Do_extendHeadTail_fixOneEndUnmapped << endl;
		log_ofs << endl << "Do_extendHeadTail_fixHeadTail: " << Do_extendHeadTail_fixHeadTail << endl; 

		log_ofs << endl << "Do_fixHeadTail_remapping: " << Do_fixHeadTail_remapping << endl;
		log_ofs << endl << "Do_fixHeadTail_greedyMapping: " << Do_fixHeadTail_greedyMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAndTargetMapping: " << Do_fixHeadTail_remappingAndTargetMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAgain: " << Do_fixHeadTail_remappingAgain << endl;

		log_ofs << endl << "minValSegLength: " << minValSegLength << endl;
		log_ofs << endl << "min_anchor_length: " << min_anchor_length << endl;
		log_ofs << endl << "CONFIDENT_SEG_LENGTH_FIX_LONG_END " << CONFIDENT_SEG_LENGTH_FIX_LONG_END << endl;
		log_ofs << endl << "DETECT_NONCANONICAL_SJ: " << DETECT_NONCANONICAL_SJ << endl;
	}

	void outputSwitchInfo(bool Do_Phase1_Only, bool outputAlignInfoAndSamForAllPairedAlignmentBool,
		bool removeAllIntermediateFilesBool, bool Do_cirRNA, bool outputDirectlyBool_Phase1Only, 
		//int normalRecordNum_1stMapping, 
		//int normalRecordNum_fixOneEndUnmapped,
		//int normalRecordNum_fixHeadTail, 
		bool Do_extendHeadTail_phase1, 
		bool Do_extendHeadTail_fixOneEndUnmapped, bool Do_extendHeadTail_fixHeadTail, 
		bool Do_fixHeadTail_remapping, bool Do_fixHeadTail_greedyMapping,
		bool Do_fixHeadTail_remappingAndTargetMapping, bool Do_fixHeadTail_remappingAgain,
		ofstream& log_ofs)
	{
		log_ofs << endl << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
		log_ofs << endl << "Do_Phase1_Only: " << Do_Phase1_Only << endl << "outputAlignInfoAndSamForAllPairedAlignmentBool: " 
			<< outputAlignInfoAndSamForAllPairedAlignmentBool << endl << "removeAllIntermediateFilesBool: " 
			<< removeAllIntermediateFilesBool << endl << "Do_cirRNA: " << Do_cirRNA << endl
			<< "outputDirectlyBool_Phase1Only: " << outputDirectlyBool_Phase1Only << endl; ;
		//log_ofs << endl << "normalRecordNum_1stMapping: " << normalRecordNum_1stMapping << endl;
		//log_ofs << "normalRecordNum_fixOneEndUnmapped: " << normalRecordNum_fixOneEndUnmapped << endl;
		//log_ofs << "normalRecordNum_fixHeadTail: " << normalRecordNum_fixHeadTail << endl;

		log_ofs << endl << "Do_extendHeadTail_phase1: " << Do_extendHeadTail_phase1 << endl;
		log_ofs << endl << "Do_extendHeadTail_fixOneEndUnmapped: " << Do_extendHeadTail_fixOneEndUnmapped << endl;
		log_ofs << endl << "Do_extendHeadTail_fixHeadTail: " << Do_extendHeadTail_fixHeadTail << endl; 

		log_ofs << endl << "Do_fixHeadTail_remapping: " << Do_fixHeadTail_remapping << endl;
		log_ofs << endl << "Do_fixHeadTail_greedyMapping: " << Do_fixHeadTail_greedyMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAndTargetMapping: " << Do_fixHeadTail_remappingAndTargetMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAgain: " << Do_fixHeadTail_remappingAgain << endl;

		log_ofs << endl << "minValSegLength: " << minValSegLength << endl;
		log_ofs << endl << "min_anchor_length: " << min_anchor_length << endl;
		log_ofs << endl << "CONFIDENT_SEG_LENGTH_FIX_LONG_END " << CONFIDENT_SEG_LENGTH_FIX_LONG_END << endl;
		log_ofs << endl << "DETECT_NONCANONICAL_SJ: " << DETECT_NONCANONICAL_SJ << endl;		
	}
};

#endif