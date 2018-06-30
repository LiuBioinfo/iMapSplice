// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// used to get softclipped base numbers distribution
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
//#include <omp.h>
#include "../general/read_block_test.h"
#include "../general/bwtmap_info.h"
#include "../general/DoubleAnchorScore.h"
#include "../general/sbndm.h"
#include "../general/splice_info.h"

using namespace std;

bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}

bool primaryOrNot(int tmpFlag)
{
	if(tmpFlag & 0x100)
		return false;
	else
		return true;
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

int main(int argc, char** argv)
{
	if(argc < 6)
	{
		cout << "Executable readLength_max read_num(NOT pair) outputFolder Aligner_1_name SAM_file_1 (Aligner_2_name SAM_file_2 ...)" << endl;
		exit(1);
	}
	string readLengthMaxStr = argv[1];
	int readLengthMax = atoi(readLengthMaxStr.c_str());
	string readNumStr = argv[2];
	unsigned int readNum = atoi(readNumStr.c_str());

	cout << "readLengthMax: " << readLengthMax << endl;
	cout << "readNum: " << readNum << endl;

	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_softClipStats_file = outputFolderStr + "softClipped_stats";
	//string output_mappedBase_totalRead_file = outputFolderStr + "mappedBase_totalRead";
	string output_unmappedBase_totalRead_file = outputFolderStr + "unmappedBase_totalRead";	
	string output_softClippedBase_mappedRead_file = outputFolderStr + "softClippedBase_mappedRead";
	//string output_mappedBase_mappedRead_file = outputFolderStr + "mappedBase_mappedRead";
	string output_softClippedBase_totalRead_file = outputFolderStr + "softClippedBase_totalRead";

	ofstream softClippedStats_ofs(output_softClipStats_file.c_str());
	//ofstream mappedBase_totalRead_ofs(output_mappedBase_totalRead_file.c_str());
	ofstream unmappedBase_totalRead_ofs(output_unmappedBase_totalRead_file.c_str());
	ofstream softClippedBase_mappedRead_ofs(output_softClippedBase_mappedRead_file.c_str());
	//ofstream mappedBase_mappedRead_ofs(output_mappedBase_mappedRead_file.c_str());
	ofstream softClippedBase_totalRead_ofs(output_softClippedBase_totalRead_file.c_str());

	cout << "start to process " << endl;

	vector<string> alignerNameVec;
	vector<string> alignmentFileVec;
	for(int tmp = 4; tmp < argc; tmp += 2)
	{
		string tmpAlignerName = argv[tmp];
		alignerNameVec.push_back(tmpAlignerName);
		string tmpAlignmentFileStr = argv[tmp+1];
		alignmentFileVec.push_back(tmpAlignmentFileStr);
	}

	int alignmentFileVecSize = alignmentFileVec.size();
	unsigned int* clippedBase_num = (unsigned int*)malloc((readLengthMax) * alignmentFileVecSize * sizeof(unsigned int));
	//unsigned int* mappedBase_num = (unsigned int*)malloc((readLengthMax) * alignmentFileVecSize * sizeof(unsigned int));	
	for(int tmp = 0; tmp < readLengthMax * alignmentFileVecSize; tmp++)
	{
		clippedBase_num[tmp] = 0;
		//mappedBase_num[tmp] = 0;
	}
	int mappedRead_num[alignmentFileVecSize];
	int unmappedRead_num[alignmentFileVecSize];
	for(int tmpAligner = 0; tmpAligner < alignmentFileVecSize; tmpAligner++)
	{	
		cout << "tmpAligner: " << tmpAligner << endl;
		//mappedRead_num[tmpAligner] = 0;
		/// count primary alignment number for each mapped length
		string tmpInputSAMfile = alignmentFileVec[tmpAligner];
		ifstream tmpInputSAM_ifs(tmpInputSAMfile.c_str());

		unsigned int primaryAlignmentNum = 0;
		while(!(tmpInputSAM_ifs.eof()))
		{
			string samStr;
			getline(tmpInputSAM_ifs, samStr);
			 if(tmpInputSAM_ifs.eof()||(samStr == ""))
			 	break;
			if(samStr.at(0) == '@')
				continue;
			//totalAlignNum ++;
			//cout << "samStr: " << samStr << endl;
			vector<string> samFieldVec;
			int startLoc = 0;
			for(int tmp = 0; tmp < 6; tmp++)
			{
				int tabLoc = samStr.find("\t", startLoc);
				string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
				samFieldVec.push_back(tmpSamField);
				startLoc = tabLoc + 1;
			}
			samFieldVec.push_back(samStr.substr(startLoc));
			string flagStr = samFieldVec[1];
			int tmpFlag = atoi(flagStr.c_str());
			string mapChrNameStr = samFieldVec[2];
			string mapChrPosStr = samFieldVec[3];
			int mapChrPos = atoi(mapChrPosStr.c_str());
			string cigarString = samFieldVec[5];			
	
			//cout << "tmpFlag" << tmpFlag << endl;
			bool mappedOrNot_bool = mappedOrNot(tmpFlag);
			bool primaryOrNot_bool = primaryOrNot(tmpFlag);
			//cout << "primary: " << primaryOrNot_bool << endl;
			//cout << "mapped: " << mappedOrNot_bool << endl;
			if(!mappedOrNot_bool)
			{
				//unmappedAlignmentNum ++;
			}
			else
			{
				if(primaryOrNot_bool)
				{
					primaryAlignmentNum ++;
					int headSoftClipLength = 0;
					int tailSoftClipLength = 0;
					vector<Jump_Code> tmpJumpCodeVec;
					//cout << "cigarString " << cigarString << endl;
					cigarString2jumpCodeVec(cigarString, tmpJumpCodeVec);
					int tmpJumpCodeVecSize = tmpJumpCodeVec.size();
					//cout << "tmpJumpCodeVecSize: " << tmpJumpCodeVecSize << endl;
					string firstJumpCodeType = tmpJumpCodeVec[0].type;
					int firstJumpCodeLength = tmpJumpCodeVec[0].len;
					string lastJumpCodeType = tmpJumpCodeVec[tmpJumpCodeVecSize-1].type;
					int lastJumpCodeLength = tmpJumpCodeVec[tmpJumpCodeVecSize-1].len;
					//cout << "firstJumpCodeLength: " << firstJumpCodeLength << endl;
					//cout << "lastJumpCodeLength: " << lastJumpCodeLength << endl;
					if(firstJumpCodeType == "S")
					{
						headSoftClipLength = firstJumpCodeLength;
					}
					if(lastJumpCodeType == "S")
					{
						tailSoftClipLength = lastJumpCodeLength;
					}
					//cout << "headSoftClipLength: " << headSoftClipLength << endl;
					//cout << "tailSoftClipLength: " << tailSoftClipLength << endl;
					// for(int tmp = 1; tmp <= headSoftClipLength; tmp++)
					// {
					// 	int tmpBasePosition = tmp;
					// 	int tmpBaseIndex = tmp - 1;
					// 	clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] = clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] + 1;
					// }
					// //cout << "s" << endl;
					// for(int tmp = 1; tmp <= tailSoftClipLength; tmp++)
					// {
					// 	int tmpBasePosition = readLengthMax - tmp + 1;
					// 	int tmpBaseIndex = tmpBasePosition - 1;
					// 	clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] = clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] + 1;
					// }
					// //cout << "s2" << endl;
					// int startMappedBasePosition = headSoftClipLength + 1;
					// int endMappedBaseposition = readLengthMax - tailSoftClipLength;
					// for(int tmpBasePosition = startMappedBasePosition; tmpBasePosition <= endMappedBaseposition; tmpBasePosition ++)
					// {
					// 	int tmpBaseIndex = tmpBasePosition - 1;
					// 	mappedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] = mappedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] + 1;
					// }
					int totalSoftClipLength_eachRead = headSoftClipLength + tailSoftClipLength;
					for(int tmp = 1; tmp <= totalSoftClipLength_eachRead; tmp++)
					{
						int tmpBasePosition = tmp;
						int tmpBaseIndex = tmp - 1;
						clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] = clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] + 1;
					}
				}
				else
				{
					//secondaryAlignNum ++;
				}
			}
		}
		mappedRead_num[tmpAligner] = primaryAlignmentNum;
		unmappedRead_num[tmpAligner] = readNum - mappedRead_num[tmpAligner];
		softClippedStats_ofs << "aligner: " << alignerNameVec[tmpAligner] << endl;
		double mapped_perc = ((double)mappedRead_num[tmpAligner]/(double)readNum)*100;
		double unmapped_perc = 100 - mapped_perc;
		softClippedStats_ofs << "mapped: " << mappedRead_num[tmpAligner] << " --- " << mapped_perc << "%" << endl;
		softClippedStats_ofs << "unmapped: " << readNum - mappedRead_num[tmpAligner] << " --- " << unmapped_perc << "%" << endl;
	}

	for(int tmpBasePosition = 1; tmpBasePosition <= readLengthMax; tmpBasePosition++)
	{
		int tmpBaseIndex = tmpBasePosition - 1;
		// mappedBase_totalRead_ofs << tmpBasePosition;
		// unmappedBase_totalRead_ofs << tmpBasePosition;
		// softClippedBase_mappedRead_ofs << tmpBasePosition;
		// mappedBase_mappedRead_ofs << tmpBasePosition;
		// softClippedBase_totalRead_ofs << tmpBasePosition;
		for(int tmpAligner = 0; tmpAligner < alignmentFileVecSize; tmpAligner++)
		{	
			double softClip_totalRead_perc = ((double)clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] / (double)readNum) * 100;
			//double mappedBase_totalRead_perc = ((double)mappedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] / (double)readNum) * 100;
			double softClip_mappedRead_perc = ((double)clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] / (double)mappedRead_num[tmpAligner]) * 100;
			//mappedBase_totalRead_ofs << tmpBasePosition << "\t" << mappedBase_totalRead_perc << "\t" << alignerNameVec[tmpAligner] << endl;
			double unmap_totalRead_perc = ((double)(clippedBase_num[tmpBaseIndex + tmpAligner * readLengthMax] + unmappedRead_num[tmpAligner]) / (double)readNum) * 100;
			
			softClippedBase_mappedRead_ofs << tmpBasePosition << "\t" << softClip_mappedRead_perc << "\t" << alignerNameVec[tmpAligner] << endl;
			//mappedBase_mappedRead_ofs << tmpBasePosition << "\t" << (100 - softClip_mappedRead_perc) << "\t" << alignerNameVec[tmpAligner] << endl;
			softClippedBase_totalRead_ofs << tmpBasePosition << "\t" << softClip_totalRead_perc << "\t" << alignerNameVec[tmpAligner] << endl;
			unmappedBase_totalRead_ofs << tmpBasePosition << "\t" << unmap_totalRead_perc << "\t" << alignerNameVec[tmpAligner] << endl;
		}
		// mappedBase_totalRead_ofs << endl;
		// unmappedBase_totalRead_ofs << endl;
		// softClippedBase_mappedRead_ofs << endl;
		// mappedBase_mappedRead_ofs << endl;
		// softClippedBase_totalRead_ofs << endl;
	}

	// for(int tmpAligner = 0; tmpAligner < alignmentFileVecSize; tmpAligner ++)
	// {
	// 	unsigned unmappedBaseNumTotal = 0;
	// 	for(int tmpBasePosition = 1; tmpBasePosition <= readLengthMax; tmpBasePosition ++)
	// 	{
	// 		int tmpBaseIndex = tmpBasePosition - 1;
	// 		unmappedBaseNumTotal = unmappedBaseNumTotal + (readNum - mappedBase_num[tmpBaseIndex + tmpAligner*readLengthMax]);
	// 	}
	// 	double unmapped_perc = ((double)unmappedBaseNumTotal / ((double)readLengthMax * readNum)) * 100;
	// 	softClippedStats_ofs << endl;
	// 	softClippedStats_ofs << alignerNameVec[tmpAligner] << "\t" << unmappedBaseNumTotal << " --- " << unmapped_perc << "%" << endl;
	// }

	return 0;
}	