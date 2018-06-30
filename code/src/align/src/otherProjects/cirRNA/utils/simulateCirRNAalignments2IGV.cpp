// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// used to simulate fake-alignments to show backSplice junctions in IGV
// generate from backSplcie junctions,
// output SAM format alignments file
// FORMAT: ChrName ChrStartPos ChrEndPos

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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputJuncFile outputFolder" << endl;
		exit(1);
	}
	int simulatedAlignmentReadSeqLength = 50;
	int anchorSeqAroundEachSpliceSite = simulatedAlignmentReadSeqLength / 2;

	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());

	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	//                classify SJs 2 backSplice and normalSplice                        ///
	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	string output_backSpliceJunc = outputFolderStr + "backSplice.junc";
	ofstream backSpliceJunc_ofs(output_backSpliceJunc.c_str());

	string output_backSpliceJunc_canonical = outputFolderStr + "backSplice_canonical.junc";
	ofstream backSpliceJunc_canonical_ofs(output_backSpliceJunc_canonical.c_str());

	string output_backSpliceJunc_semicanonical = outputFolderStr + "backSplice_semicanonical.junc";
	ofstream backSpliceJunc_semicanonical_ofs(output_backSpliceJunc_semicanonical.c_str());

	string output_backSpliceJunc_noncanonical = outputFolderStr + "backSplice_noncanonical.junc";
	ofstream backSpliceJunc_noncanonical_ofs(output_backSpliceJunc_noncanonical.c_str());

	string output_normalSpliceJunc = outputFolderStr + "normalSplice.junc";
	ofstream normalSpliceJunc_ofs(output_normalSpliceJunc.c_str());

	string output_backSpliceJunc_simAlignment = outputFolderStr + "simulatedBackSplice.sam";
	ofstream backSpliceJunc_simSAM_ofs(output_backSpliceJunc_simAlignment.c_str());

	string output_backSpliceJunc_simAlignment_canonical = outputFolderStr + "simulatedBackSplice_canonical.sam";
	ofstream backSpliceJunc_simSAM_canonical_ofs(output_backSpliceJunc_simAlignment_canonical.c_str());

	string output_backSpliceJunc_simAlignment_semicanonical = outputFolderStr + "simulatedBackSplice_semicanonical.sam";
	ofstream backSpliceJunc_simSAM_semicanonical_ofs(output_backSpliceJunc_simAlignment_semicanonical.c_str());

	string output_backSpliceJunc_simAlignment_noncanonical = outputFolderStr + "simulatedBackSplice_noncanonical.sam";
	ofstream backSpliceJunc_simSAM_noncanonical_ofs(output_backSpliceJunc_simAlignment_noncanonical.c_str());	

	string output_normalSpliceJunc_simAlignment = outputFolderStr + "simulatedNormalSplice.sam";
	ofstream normalSpliceJunc_simSAM_ofs(output_normalSpliceJunc_simAlignment.c_str());

	string inputJuncFile = argv[2];
	ifstream junc_ifs(inputJuncFile.c_str());
	while(!(junc_ifs.eof()))
	{	
		string juncStr;
		getline(junc_ifs, juncStr);
		 if(junc_ifs.eof()||(juncStr == ""))
		 	break;

		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));
		string tmpSJchrName = juncFieldVec[0];
		string tmpSJstartPosStr = juncFieldVec[1];
		string tmpSJendPosStr = juncFieldVec[2];
		string tmpSJnameStr = juncFieldVec[3];
		string tmpSJsupportNumStr = juncFieldVec[4];
		int tmpSJchrNameInt = indexInfo->convertStringToInt(tmpSJchrName);
		int tmpSJstartPosInt = atoi(tmpSJstartPosStr.c_str());
		int tmpSJendPosInt = atoi(tmpSJendPosStr.c_str());
		int tmpSJsupportNumInt = atoi(tmpSJsupportNumStr.c_str());
		//cout << "tmpSJstartPosInt: " << tmpSJstartPosInt << endl;
		//cout << "tmpSJendPosInt: " << tmpSJendPosInt << endl;
		if(tmpSJstartPosInt < tmpSJendPosInt)
		{
			normalSpliceJunc_ofs << juncStr << endl;
			int startPos_1stSegment = tmpSJstartPosInt - anchorSeqAroundEachSpliceSite + 1;
			int endPos_1stSegment = tmpSJstartPosInt;
			int startPos_2ndSegment = tmpSJendPosInt;
			int endPos_2ndSegment = tmpSJendPosInt + anchorSeqAroundEachSpliceSite - 1;
			string seq_1stSegment = indexInfo->returnChromStrSubstr(
				tmpSJchrNameInt, startPos_1stSegment, anchorSeqAroundEachSpliceSite);
			string seq_2ndSegment = indexInfo->returnChromStrSubstr(
				tmpSJchrNameInt, startPos_2ndSegment, anchorSeqAroundEachSpliceSite);
			string simulatedAlignmentReadSeq = seq_1stSegment + seq_2ndSegment;
			for(int tmp = 0; tmp < tmpSJsupportNumInt; tmp++)
			{
				string tmpReadSam = "seq." + tmpSJnameStr + "_" + int_to_str(tmp+1)
					+ "\t0\t" + tmpSJchrName + "\t" + int_to_str(startPos_1stSegment) + "\t255\t"
					+ int_to_str(anchorSeqAroundEachSpliceSite) + "M" 
					+ int_to_str(startPos_2ndSegment - endPos_1stSegment - 1) + "N"
					+ int_to_str(anchorSeqAroundEachSpliceSite) + "M" + "\t"
					+ "*\t0\t0\t" + simulatedAlignmentReadSeq + "\t*\tNM:i:0\tIH:i:1\tHI:i:1";
				normalSpliceJunc_simSAM_ofs << tmpReadSam << endl;
			}
		}
		else if(tmpSJstartPosInt > tmpSJendPosInt)
		{
			int tmpBackSpliceSJflankStringCase = indexInfo->returnFlankStringCase(
				tmpSJchrNameInt, tmpSJstartPosInt, tmpSJendPosInt);
			if(tmpBackSpliceSJflankStringCase >= 5)
				backSpliceJunc_canonical_ofs << juncStr << endl;
			else if(tmpBackSpliceSJflankStringCase > 0)
				backSpliceJunc_semicanonical_ofs << juncStr << endl;
			else
			{
				backSpliceJunc_noncanonical_ofs << juncStr << endl;
			}


			backSpliceJunc_ofs << juncStr << endl;
			int startPos_1stSegment = tmpSJendPosInt - 1 - anchorSeqAroundEachSpliceSite + 1;
			int endPos_1stSegment = tmpSJendPosInt - 1;
			int startPos_2ndSegment = tmpSJstartPosInt + 1;
			int endPos_2ndSegment = tmpSJstartPosInt + 1 + anchorSeqAroundEachSpliceSite - 1;
			string seq_1stSegment = indexInfo->returnChromStrSubstr(
				tmpSJchrNameInt, startPos_1stSegment, anchorSeqAroundEachSpliceSite);
			string seq_2ndSegment = indexInfo->returnChromStrSubstr(
				tmpSJchrNameInt, startPos_2ndSegment, anchorSeqAroundEachSpliceSite);
			string simulatedAlignmentReadSeq = seq_1stSegment + seq_2ndSegment;
			for(int tmp = 0; tmp < tmpSJsupportNumInt; tmp++)
			{
				string tmpReadSam = "seq." + tmpSJnameStr + "_" + int_to_str(tmp+1)
					+ "\t0\t" + tmpSJchrName + "\t" + int_to_str(startPos_1stSegment) + "\t255\t"
					+ int_to_str(anchorSeqAroundEachSpliceSite) + "M" 
					+ int_to_str(startPos_2ndSegment - endPos_1stSegment - 1) + "N"
					+ int_to_str(anchorSeqAroundEachSpliceSite) + "M" + "\t"
					+ "*\t0\t0\t" + simulatedAlignmentReadSeq + "\t*\tNM:i:0\tIH:i:1\tHI:i:1";
				backSpliceJunc_simSAM_ofs << tmpReadSam << endl;
				if(tmpBackSpliceSJflankStringCase >= 5)
					backSpliceJunc_simSAM_canonical_ofs << tmpReadSam << endl;
				else if(tmpBackSpliceSJflankStringCase > 0)
					backSpliceJunc_simSAM_semicanonical_ofs << tmpReadSam << endl;
				else
					backSpliceJunc_simSAM_noncanonical_ofs << tmpReadSam << endl;
			}
		}
		else
		{
			cout << "error ! tmpSJstartPosInt == tmpSJendPosInt" << endl;
			exit(1);
		}
	}
	backSpliceJunc_ofs.close();
	normalSpliceJunc_ofs.close();
	backSpliceJunc_simSAM_ofs.close();
	return 0;
}