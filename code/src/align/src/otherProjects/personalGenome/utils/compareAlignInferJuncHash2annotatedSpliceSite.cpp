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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <InputIndexFolderPath> <SJfromAnnotationFile> <ToCompareJuncFile> <outputFolder> <inputSNPfile>" << endl;
		exit(1);
	}
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "SJfromAnnotationFile: " << argv[2] << endl;
	log_ofs << "ToCompareJuncFile: " << argv[3] << endl;
	log_ofs << "outputFolder: " << argv[4] << endl;

	cout << "defining output files ......" << endl;
	log_ofs << endl << "defining output files ......" << endl;
	string outputCompareFileStr_total = outputFolderStr + "total.compare";
	ofstream compare_ofs(outputCompareFileStr_total.c_str());
	string annotated_bothEnds_JuncFile = outputFolderStr + "annotated_bothEnds.junc";
	string annotated_oneEnd_JuncFile = outputFolderStr + "annotated_oneEnd.junc";
	string annotated_oneEnd_JuncFile_withOtherJuncInfo = outputFolderStr + "annotated_oneEnd_detailInfo.junc";
	string unannotated_JuncFile = outputFolderStr + "unannotated.junc";

	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	string indexStr = indexParameterFileStr;
	indexStr += "/";
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of loading indexes" << endl;

	cout << "generating 2 sam2alignInferJuncHash" << endl;
	log_ofs << "generating 2 sam2alignInferJuncHash" << endl;
	string juncHash_file_1_str = argv[2];
	cout << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	log_ofs << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	AlignInferJunctionHash_Info* juncHash_1 = new AlignInferJunctionHash_Info();
	juncHash_1->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to read 2 juncfiles ...." << endl;
	log_ofs << "start to read 2 juncfiles ...." << endl;
	juncHash_1->insertJuncFromJuncFile_chrNamePosOnly(juncHash_file_1_str, indexInfo);
	cout << "start to generate SJ " << endl;
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	juncHash_1->convert2SJhashInfo(SJ, indexInfo);
	cout << "start to generate spliceSiteSetVec" << endl;
	vector< set<int> > spliceSiteSetVec_doner;
	vector< set<int> > spliceSiteSetVec_acceptor;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpSpliceSiteSet_doner;
		set<int> tmpSpliceSiteSet_acceptor;
		spliceSiteSetVec_doner.push_back(tmpSpliceSiteSet_doner);
		spliceSiteSetVec_acceptor.push_back(tmpSpliceSiteSet_acceptor);
	}
	juncHash_1->generateSpliceSiteSetVec(spliceSiteSetVec_doner, spliceSiteSetVec_acceptor);

	ofstream bothEndsAnnotated_ofs(annotated_bothEnds_JuncFile.c_str());
	ofstream oneEndAnnotated_ofs(annotated_oneEnd_JuncFile.c_str());
	ofstream oneEndAnnotated_detail_ofs(annotated_oneEnd_JuncFile_withOtherJuncInfo.c_str());
	ofstream unannotated_ofs(unannotated_JuncFile.c_str());
	string juncHash_file_2_str = argv[3];
	ifstream junc_2_ifs(juncHash_file_2_str.c_str());
	while(!junc_2_ifs.eof())
	{
		string tmpStr;
		getline(junc_2_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpChrName = tmpStr.substr(0, tabLoc_1);
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		string tmpStartPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpEndPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		int tmpEndPos = atoi(tmpEndPosStr.c_str());
		bool startPosFoundBool = ((spliceSiteSetVec_doner[tmpChrNameInt]).find(tmpStartPos) 
			!= (spliceSiteSetVec_doner[tmpChrNameInt]).end());
		bool endPosFoundBool = ((spliceSiteSetVec_acceptor[tmpChrNameInt]).find(tmpEndPos) 
			!= (spliceSiteSetVec_acceptor[tmpChrNameInt]).end());		
		if(startPosFoundBool && endPosFoundBool)
			bothEndsAnnotated_ofs << tmpStr << endl;
		else if(startPosFoundBool || endPosFoundBool)
		{
			oneEndAnnotated_detail_ofs << "AlterSpliceJunc:";
			if(startPosFoundBool)
			{
				vector<int> tmpAlterAcceptorSpliceSiteVec;
				SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(
					tmpChrNameInt, tmpStartPos, tmpAlterAcceptorSpliceSiteVec);
				for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSiteVec.size(); tmp++)
					oneEndAnnotated_detail_ofs << "\t" << tmpStartPos 
						<< ":" << tmpAlterAcceptorSpliceSiteVec[tmp];
			}
			else
			{
				vector<int> tmpAlterDonerSpliceSiteVec;
				SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(
					tmpChrNameInt, tmpEndPos, tmpAlterDonerSpliceSiteVec);
				for(int tmp = 0; tmp < tmpAlterDonerSpliceSiteVec.size(); tmp++)
					oneEndAnnotated_detail_ofs << "\t" << tmpAlterDonerSpliceSiteVec[tmp]
						<< ":" << tmpEndPos;
			}
			oneEndAnnotated_detail_ofs << endl << tmpStr << endl;
			oneEndAnnotated_ofs << tmpStr << endl;;
		}
		else
			unannotated_ofs << tmpStr << endl;
	}
	junc_2_ifs.close();
	unannotated_ofs.close();
	oneEndAnnotated_ofs.close();
	oneEndAnnotated_detail_ofs.close();
	bothEndsAnnotated_ofs.close();
	cout << "All jobs done ! " << endl;
	log_ofs << "All jobs done ! " << endl;	
	delete juncHash_1;
	delete indexInfo;
	compare_ofs.close();
	parameter_ifs.close();
	log_ofs.close();
}	