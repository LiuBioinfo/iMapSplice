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
#include "../../../general/alignInferJunctionHash_info.h"
#include "./../general/SNPhash_info.h"

using namespace std;

int returnSJindex(string& tmpSJgroup_geneName, string& tmpSJgroup_strand, int tmpSJgroup_chrNameInt, 
	vector< vector<string> >& geneNameVecVec_ori, vector< vector<string> >& geneStrandVecVec_ori)
{
	//cout << "tmpSJgroup_geneName: " << tmpSJgroup_geneName << endl;
	//cout << "tmpSJgroup_strand: " << tmpSJgroup_strand << endl;
	//cout << "tmpSJgroup_chrNameInt: " << tmpSJgroup_chrNameInt << endl;
	int tmpChrGeneNameVecSize = geneNameVecVec_ori[tmpSJgroup_chrNameInt].size();
	//cout << "tmpChrGeneNameVecSize: " << tmpChrGeneNameVecSize << endl;
	for(int tmp = 0; tmp < tmpChrGeneNameVecSize; tmp++)
	{
		string tmpGeneName_existing = (geneNameVecVec_ori[tmpSJgroup_chrNameInt])[tmp];
		string tmpGeneStrand_existing = (geneStrandVecVec_ori[tmpSJgroup_chrNameInt])[tmp];
		//cout << "tmpGeneName_existing: " << tmpGeneName_existing << endl;
		//cout << "tmpGeneStrand_existing: " << tmpGeneStrand_existing << endl;
		if((tmpSJgroup_geneName == tmpGeneName_existing)
			&&(tmpSJgroup_strand == tmpGeneStrand_existing))
			return tmp;
	}
	return -1;
}

void generateRawGeneInfoVecVecFromGTFfile(string& gtfFile,
	vector< vector<string> >& geneNameVecVec_ori, vector< vector<string> >& geneStrandVecVec_ori,
	vector< vector< pair<int,int> > >& geneBoundaryPosPairVecVec_ori, Index_Info* indexInfo)
{
	ifstream gtf_ifs(gtfFile.c_str());
	while(!gtf_ifs.eof())
	{
		string tmpStr;
		getline(gtf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.substr(0,1) == "#")
			continue;
		vector<string> gtfFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpGtfField = tmpStr.substr(startLoc, tabLoc-startLoc);
			gtfFieldVec.push_back(tmpGtfField);
			startLoc = tabLoc + 1;
		}
		string tmpChrName = gtfFieldVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
			continue;	
		string tmpFeature = gtfFieldVec[2];
		if(tmpFeature != "gene")
			continue;
		int tmpStartPos = atoi((gtfFieldVec[3]).c_str());
		int tmpEndPos = atoi((gtfFieldVec[4]).c_str());
		string tmpStrand = gtfFieldVec[6];
		int tmpGeneNameLocInStr = tmpStr.find("gene_name");
		if(tmpGeneNameLocInStr == string::npos)
			continue;
		int firstQuota_geneName = tmpStr.find("\"", tmpGeneNameLocInStr + 1);
		int secondQuota_geneName = tmpStr.find("\"", firstQuota_geneName + 1);
		if((firstQuota_geneName == string::npos)||(secondQuota_geneName == string::npos))
			continue;
		string tmpGeneName = tmpStr.substr(firstQuota_geneName + 1,
			secondQuota_geneName - firstQuota_geneName - 1);
		geneNameVecVec_ori[tmpChrNameInt].push_back(tmpGeneName);
		geneStrandVecVec_ori[tmpChrNameInt].push_back(tmpStrand);
		geneBoundaryPosPairVecVec_ori[tmpChrNameInt].push_back(pair<int,int>(tmpStartPos, tmpEndPos));		
	}
	gtf_ifs.close();
}

void generateRawSJgroupInfoVecFromSJfile(string& SJfile, vector<string>& SJgroup_geneNameVec, 
	vector<string>& SJgroup_strandVec, vector<int>& SJgroup_chrNameIntVec, 
	vector< vector< pair<int,int> > >& SJgroup_posPairVecVec, Index_Info* indexInfo)
{
	//cout << "start to do generateRawSJgroupInfoVecFromSJfile" << endl;
	ifstream SJ_ifs(SJfile.c_str());
	while(!SJ_ifs.eof())
	{
		string tmpStr;
		getline(SJ_ifs, tmpStr);
		if(tmpStr == "")
			break;		
		vector<string> SJfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 5; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpSJfield = tmpStr.substr(startLoc, tabLoc-startLoc);
			SJfieldVec.push_back(tmpSJfield);
			startLoc = tabLoc + 1;
		}
		SJfieldVec.push_back(tmpStr.substr(startLoc));
		string tmpSJ_chrName = SJfieldVec[0];
		string tmpSJ_startPosStr = SJfieldVec[1];
		string tmpSJ_endPosStr = SJfieldVec[2];
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrName);
		int tmpSJ_startPos = atoi(tmpSJ_startPosStr.c_str());
		int tmpSJ_endPos = atoi(tmpSJ_endPosStr.c_str());
		string tmpSJ_geneName_withSemiColon = SJfieldVec[5];
		int tmpSJ_geneName_withSemiColon_length = tmpSJ_geneName_withSemiColon.length();
		string tmpSJ_geneName = tmpSJ_geneName_withSemiColon.substr(0, tmpSJ_geneName_withSemiColon_length - 1);
		string tmpSJ_strand = SJfieldVec[3];

		if(SJgroup_geneNameVec.size() == 0)
		{
			SJgroup_geneNameVec.push_back(tmpSJ_geneName);
			SJgroup_strandVec.push_back(tmpSJ_strand);
			SJgroup_chrNameIntVec.push_back(tmpSJ_chrNameInt);
			vector< pair<int,int> > tmpNewSJgroupPosPairVec;
			tmpNewSJgroupPosPairVec.push_back(pair<int,int>(tmpSJ_startPos, tmpSJ_endPos));
			SJgroup_posPairVecVec.push_back(tmpNewSJgroupPosPairVec);
			continue;			
		}

		int currentLastSJindex = SJgroup_geneNameVec.size() - 1;
		string currentLastSJ_geneName = SJgroup_geneNameVec[currentLastSJindex];
		string currentLastSJ_strand = SJgroup_strandVec[currentLastSJindex];
		int currentLastSJ_chrNameInt = SJgroup_chrNameIntVec[currentLastSJindex];
		if((tmpSJ_geneName == currentLastSJ_geneName)
			&&(tmpSJ_strand == currentLastSJ_strand)
			&&(tmpSJ_chrNameInt == currentLastSJ_chrNameInt))
			SJgroup_posPairVecVec[currentLastSJindex].push_back(pair<int,int>(tmpSJ_startPos, tmpSJ_endPos));
		else
		{
			SJgroup_geneNameVec.push_back(tmpSJ_geneName);
			SJgroup_strandVec.push_back(tmpSJ_strand);
			SJgroup_chrNameIntVec.push_back(tmpSJ_chrNameInt);
			vector< pair<int,int> > tmpNewSJgroupPosPairVec;
			tmpNewSJgroupPosPairVec.push_back(pair<int,int>(tmpSJ_startPos, tmpSJ_endPos));
			SJgroup_posPairVecVec.push_back(tmpNewSJgroupPosPairVec);
		}
	}
	SJ_ifs.close();
	//cout << "end of doing generateRawSJgroupInfoVecFromSJfile" << endl;
}

void generateGeneSpliceSiteVecFromGTF(string& gtfFile, string& SJfile,
		vector<string>& geneNameVec, vector<string>& geneStrandVec,
		vector<int>& geneChrNameIntVec, vector< pair<int,int> >& geneBoundaryPosPairVec,
		vector< vector<int> >& geneSpliceDonerSiteVecVec, vector< vector<int> >& geneSpliceAcceptorSiteVecVec, 
		Index_Info* indexInfo)
{
	//cout << "generateGeneSpliceSiteVecFromGTF starts ......" << endl; 
	// read GTF file to generate original gene info vec
	vector< vector<string> > geneNameVecVec_ori;
	vector< vector<string> > geneStrandVecVec_ori;
	vector< vector< pair<int,int> > > geneBoundaryPosPairVecVec_ori;
	int chromNum = indexInfo->returnChromNum();
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		vector<string> tmpGeneNameVec;
		vector<string> tmpGeneStrandVec;
		vector< pair<int,int> > tmpGeneBoundaryPosPairVec;
		geneNameVecVec_ori.push_back(tmpGeneNameVec);
		geneStrandVecVec_ori.push_back(tmpGeneStrandVec);
		geneBoundaryPosPairVecVec_ori.push_back(tmpGeneBoundaryPosPairVec);
	}
	//cout << "generateRawGeneInfoVecVecFromGTFfile starts ......" << endl;
	generateRawGeneInfoVecVecFromGTFfile(gtfFile, geneNameVecVec_ori, 
		geneStrandVecVec_ori, geneBoundaryPosPairVecVec_ori, indexInfo);
	//for(int tmp = 0; tmp < chromNum; tmp++)
	//	cout << "chrom: " << tmp << " geneNum: " << (geneNameVecVec_ori[tmp]).size() << endl;
	// read SJ file to generate original SJ info vec
	vector<string> SJgroup_geneNameVec;
	vector<string> SJgroup_strandVec;
	vector<int> SJgroup_chrNameIntVec;
	vector< vector< pair<int,int> > > SJgroup_posPairVecVec;
	//cout << "generateRawSJgroupInfoVecFromSJfile starts ......" << endl;
	generateRawSJgroupInfoVecFromSJfile(SJfile, SJgroup_geneNameVec,
		SJgroup_strandVec, SJgroup_chrNameIntVec, SJgroup_posPairVecVec, indexInfo);
	// based on original SJ info vec, generate original spliceSite vec
	// SJ_posPairVecVec -> SJ_donerSiteVecVec, SJ_acceptorSiteVecVec
	// remove duplicate doner / acceptor sites
	//cout << "start to generate SJgroup_donerSiteVecVec and SJgroup_acceptorSiteVecVec " << endl;
	//cout << "SJgroup_geneNameVec.size(): " << SJgroup_geneNameVec.size() << endl;
	//cout << "SJgroup_posPairVecVec.size(): " << SJgroup_posPairVecVec.size() << endl;
	vector< vector<int> > SJgroup_donerSiteVecVec;
	vector< vector<int> > SJgroup_acceptorSiteVecVec;
	for(int tmp = 0; tmp < SJgroup_posPairVecVec.size(); tmp++)
	{
		//cout << "tmp: " << tmp << endl;
		vector<int> tmpSJdonerSiteVec;
		vector<int> tmpSJacceptorSiteVec;
		int tmpSJposPairVecSize = SJgroup_posPairVecVec[tmp].size();
		//cout << "tmpSJposPairVecSize " << tmpSJposPairVecSize << endl; 
 		for(int tmp2 = 0; tmp2 < tmpSJposPairVecSize; tmp2 ++)
		{
			int tmpSJdonerSite_raw = (SJgroup_posPairVecVec[tmp])[tmp2].first;
			int currentTmpSJdonerSiteVecSize = tmpSJdonerSiteVec.size();
			bool theSameAsOneExistingSJdonerSiteBool = false;
			for(int tmp3 = 0; tmp3 < currentTmpSJdonerSiteVecSize; tmp3++)
			{
				int tmpSJdonerSite_existing = tmpSJdonerSiteVec[tmp3];
				if(tmpSJdonerSite_raw == tmpSJdonerSite_existing)
				{
					theSameAsOneExistingSJdonerSiteBool = true;
					break;
				}
			}
			if(!theSameAsOneExistingSJdonerSiteBool)
				tmpSJdonerSiteVec.push_back(tmpSJdonerSite_raw);
			
			int tmpSJacceptorSite_raw = (SJgroup_posPairVecVec[tmp])[tmp2].second;
			int currentTmpSJacceptorSiteVecSize = tmpSJacceptorSiteVec.size();
			bool theSameAsOneExistingSJacceptorSiteBool = false;
			for(int tmp3 = 0; tmp3 < currentTmpSJacceptorSiteVecSize; tmp3++)
			{
				int tmpSJacceptorSite_existing = tmpSJacceptorSiteVec[tmp3];
				if(tmpSJacceptorSite_raw == tmpSJacceptorSite_existing)
				{
					theSameAsOneExistingSJacceptorSiteBool = true;
					break;
				}
			}
			if(!theSameAsOneExistingSJacceptorSiteBool)
				tmpSJacceptorSiteVec.push_back(tmpSJacceptorSite_raw);
		}
		//cout << "tmpSJdonerSiteVec.size():" << tmpSJdonerSiteVec.size() << endl;
		//cout << "tmpSJacceptorSiteVec.size():" << tmpSJacceptorSiteVec.size() << endl;
 		SJgroup_donerSiteVecVec.push_back(tmpSJdonerSiteVec);
		SJgroup_acceptorSiteVecVec.push_back(tmpSJacceptorSiteVec);
	}

	// generate final vec ...., search for matched gene name
	//cout << "start to generate final vec ......" << endl;
	int SJgroupNum = SJgroup_geneNameVec.size();
	//cout << "SJgroupNum: " << SJgroupNum << endl;
	for(int tmp = 0; tmp < SJgroupNum; tmp++)
	{
		string tmpSJgroup_geneName = SJgroup_geneNameVec[tmp];
		string tmpSJgroup_strand = SJgroup_strandVec[tmp];
		int tmpSJgroup_chrNameInt = SJgroup_chrNameIntVec[tmp];
		// check with geneNameVecVec_ori, geneStrandVecVec_ori, geneBoundaryPosPairVecVec_ori
		int tmpIndex = returnSJindex(tmpSJgroup_geneName, tmpSJgroup_strand, tmpSJgroup_chrNameInt,
			geneNameVecVec_ori, geneStrandVecVec_ori);
		//cout << "tmpIndex: " << tmpIndex << endl;
		if(tmpIndex < 0)
			continue;
		//cout << "tmpIndex: " << tmpIndex << endl;
		int tmpGene_startPos = (geneBoundaryPosPairVecVec_ori[tmpSJgroup_chrNameInt])[tmpIndex].first;
		int tmpGene_endPos = (geneBoundaryPosPairVecVec_ori[tmpSJgroup_chrNameInt])[tmpIndex].second;
		
		geneNameVec.push_back(tmpSJgroup_geneName);
		geneStrandVec.push_back(tmpSJgroup_strand);
		geneChrNameIntVec.push_back(tmpSJgroup_chrNameInt);
		geneBoundaryPosPairVec.push_back(pair<int,int>(tmpGene_startPos, tmpGene_endPos));
		geneSpliceDonerSiteVecVec.push_back(SJgroup_donerSiteVecVec[tmp]);
		geneSpliceAcceptorSiteVecVec.push_back(SJgroup_acceptorSiteVecVec[tmp]);
	}
	//cout << "geneNameVec.size(): " << geneNameVec.size() << endl; 
}

int returnLeftMostPos(vector<int>& tmpIntVec)
{
	int tmpLeftMostPos = tmpIntVec[0];
	if(tmpIntVec.size() > 1)
	{
		for(int tmp = 1; tmp < tmpIntVec.size(); tmp++)
		{
			int tmpInt = tmpIntVec[tmp];
			if(tmpInt < tmpLeftMostPos)
				tmpLeftMostPos = tmpInt;
		}
	}
	return tmpLeftMostPos;
}

int returnRightMostPos(vector<int>& tmpIntVec)
{
	int tmpRightMostPos = tmpIntVec[0];
	if(tmpIntVec.size() > 1)
	{
		for(int tmp = 1; tmp < tmpIntVec.size(); tmp++)
		{
			int tmpInt = tmpIntVec[tmp];
			if(tmpInt > tmpRightMostPos)
				tmpRightMostPos = tmpInt;
		}
	}
	return tmpRightMostPos;
}

void generatePsSJsiteVec(vector<int>& tmpGenePsSpliceDonerSiteVec, vector<int>& tmpGenePsSpliceAcceptorSiteVec, 
	int tmpGeneChrNameInt, int tmpGeneStartPos, int tmpGeneEndPos,
	vector<int>& tmpGeneSpliceDonerSiteVec, vector<int>& tmpGeneSpliceAcceptorSiteVec, 
	vector<int>& tmpSNPposVecInGene, string& tmpGeneStrand, Index_Info* indexInfo_ps)
{
	string canonicalFlankString_doner = "GT";
	string canonicalFlankString_acceptor = "AG";
	if(tmpGeneStrand == "-")
	{
		canonicalFlankString_doner = "CT";
		canonicalFlankString_acceptor = "AC";
	}
	int tmpGeneSpliceDonerSite_leftMost = returnLeftMostPos(tmpGeneSpliceDonerSiteVec);
	int tmpGeneSpliceAcceptorSite_rightMost = returnRightMostPos(tmpGeneSpliceAcceptorSiteVec);

	int tmpSNPposVecInGeneSize = tmpSNPposVecInGene.size();
	// generate tmpCandiPsDonerSiteVec with "GT" or "CT" flank string.	
	// generate tmpCandiPsAcceptorSiteVec with "AG" or "AC" flank string.
	vector<int> tmpCandiPsDonerSiteVec;
	vector<int> tmpCandiPsAcceptorSiteVec;
	for(int tmp = 0; tmp < tmpSNPposVecInGeneSize; tmp++)
	{
		int tmpSNPpos = tmpSNPposVecInGene[tmp];
		
		int tmpCandiPsDonerSite_1 = tmpSNPpos - 2;
		string tmpCandiPsDonerSite_1_twoBaseString 
			= indexInfo_ps->returnTwoBasesString(tmpGeneChrNameInt, tmpCandiPsDonerSite_1 + 1);
		if((tmpCandiPsDonerSite_1_twoBaseString == canonicalFlankString_doner)
			&&(tmpCandiPsDonerSite_1 < tmpGeneSpliceAcceptorSite_rightMost - 1))
			tmpCandiPsDonerSiteVec.push_back(tmpCandiPsDonerSite_1);
		
		int tmpCandiPsDonerSite_2 = tmpSNPpos - 1;
		string tmpCandiPsDonerSite_2_twoBaseString 
			= indexInfo_ps->returnTwoBasesString(tmpGeneChrNameInt, tmpCandiPsDonerSite_2 + 1);		
		if((tmpCandiPsDonerSite_2_twoBaseString == canonicalFlankString_doner)
			&&(tmpCandiPsDonerSite_2 < tmpGeneSpliceAcceptorSite_rightMost - 1))
			tmpCandiPsDonerSiteVec.push_back(tmpCandiPsDonerSite_2);		
		
		int tmpCandiPsAcceptorSite_1 = tmpSNPpos + 2;
		string tmpCandiPsAcceptorSite_1_twoBaseString
			= indexInfo_ps->returnTwoBasesString(tmpGeneChrNameInt, tmpCandiPsAcceptorSite_1 - 2);
		if((tmpCandiPsAcceptorSite_1_twoBaseString == canonicalFlankString_acceptor)
			&&(tmpCandiPsAcceptorSite_1 > tmpGeneSpliceDonerSite_leftMost + 1))
			tmpCandiPsAcceptorSiteVec.push_back(tmpCandiPsAcceptorSite_1);
		
		int tmpCandiPsAcceptorSite_2 = tmpSNPpos + 1;
		string tmpCandiPsAcceptorSite_2_twoBaseString
			= indexInfo_ps->returnTwoBasesString(tmpGeneChrNameInt, tmpCandiPsAcceptorSite_2 - 2);
		if((tmpCandiPsAcceptorSite_2_twoBaseString == canonicalFlankString_acceptor)
			&&(tmpCandiPsAcceptorSite_2 > tmpGeneSpliceDonerSite_leftMost + 1))
			tmpCandiPsAcceptorSiteVec.push_back(tmpCandiPsAcceptorSite_2);
	}

	//filter those ps sites already in annotation
	vector<int> tmpCandiPsDonerSiteVec_novel;
	vector<int> tmpCandiPsAcceptorSiteVec_novel;	
	for(int tmp = 0; tmp < tmpCandiPsDonerSiteVec.size(); tmp++)
	{
		int tmpCandiPsDonerSite = tmpCandiPsDonerSiteVec[tmp];
		bool existingBool = false;
		for(int tmp2 = 0; tmp2 < tmpGeneSpliceDonerSiteVec.size(); tmp2++)
		{
			int tmpGeneSpliceDonerSite = tmpGeneSpliceDonerSiteVec[tmp2];
			if(tmpCandiPsDonerSite == tmpGeneSpliceDonerSite)
			{
				existingBool = true;
				break;
			}
		}
		if(!existingBool)
			tmpCandiPsDonerSiteVec_novel.push_back(tmpCandiPsDonerSite);
	}
	for(int tmp = 0; tmp < tmpCandiPsAcceptorSiteVec.size(); tmp++)
	{
		int tmpCandiPsAcceptorSite = tmpCandiPsAcceptorSiteVec[tmp];
		bool existingBool = false;
		for(int tmp2 = 0; tmp2 < tmpGeneSpliceAcceptorSiteVec.size(); tmp2++)
		{
			int tmpGeneSpliceAcceptorSite = tmpGeneSpliceAcceptorSiteVec[tmp2];
			if(tmpCandiPsAcceptorSite == tmpGeneSpliceAcceptorSite)
			{
				existingBool = true;
				break;
			}
		}
		if(!existingBool)
			tmpCandiPsAcceptorSiteVec_novel.push_back(tmpCandiPsAcceptorSite);
	}

	// filter those duplicate sites
	for(int tmp = 0; tmp < tmpCandiPsDonerSiteVec_novel.size(); tmp++)
	{
		int tmpCandiPsDonerSite_novel = tmpCandiPsDonerSiteVec_novel[tmp];
		int currentPsDonerSiteVecSize = tmpGenePsSpliceDonerSiteVec.size();
		bool alreadyAddedBool = false;
		for(int tmp2 = 0; tmp2 < currentPsDonerSiteVecSize; tmp2++)
		{
			int tmpCurrentPsDonerSite = tmpGenePsSpliceDonerSiteVec[tmp2];
			if(tmpCandiPsDonerSite_novel == tmpCurrentPsDonerSite)
			{
				alreadyAddedBool = true;
				break;
			}
		}
		if(!alreadyAddedBool)
			tmpGenePsSpliceDonerSiteVec.push_back(tmpCandiPsDonerSite_novel);
	}
	for(int tmp = 0; tmp < tmpCandiPsAcceptorSiteVec_novel.size(); tmp++)
	{
		int tmpCandiPsAcceptorSite_novel = tmpCandiPsAcceptorSiteVec_novel[tmp];
		int currentPsAcceptorSiteVecSize = tmpGenePsSpliceAcceptorSiteVec.size();
		bool alreadyAddedBool = false;
		for(int tmp2 = 0; tmp2 < currentPsAcceptorSiteVecSize; tmp2++)
		{
			int tmpCurrentPsAcceptorSite = tmpGenePsSpliceAcceptorSiteVec[tmp2];
			if(tmpCandiPsAcceptorSite_novel == tmpCurrentPsAcceptorSite)
			{
				alreadyAddedBool = true;
				break;				
			}
		}
		if(!alreadyAddedBool)
			tmpGenePsSpliceAcceptorSiteVec.push_back(tmpCandiPsAcceptorSite_novel);
	}
}

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputSdIndex inputPsIndex inputSNPfile inputGTFannFile inputSJfile outputPsSJinfoFilePrefix" << endl;
		exit(1);
	}
	cout << "start to initiate indexInfo for both sd and ps genome" << endl;
	string indexFolderPath_sd = argv[1];
	string indexFolderPath_ps = argv[2];
	indexFolderPath_sd += "/";
	indexFolderPath_ps += "/";
	string chrom_bit_file_sd = indexFolderPath_sd; chrom_bit_file_sd.append("_chrom");
	string chrom_bit_file_ps = indexFolderPath_ps; chrom_bit_file_ps.append("_chrom"); 
	ifstream chrom_bit_file_ifs_sd(chrom_bit_file_sd.c_str(),ios::binary);
	ifstream chrom_bit_file_ifs_ps(chrom_bit_file_ps.c_str(),ios::binary);	
	string indexParameterFileStr_sd = indexFolderPath_sd + "/_parameter";
	string indexParameterFileStr_ps = indexFolderPath_ps + "/_parameter";
	ifstream parameter_ifs_sd(indexParameterFileStr_sd.c_str());
	ifstream parameter_ifs_ps(indexParameterFileStr_ps.c_str());
	Index_Info* indexInfo_sd = new Index_Info(parameter_ifs_sd);
	Index_Info* indexInfo_ps = new Index_Info(parameter_ifs_ps);
	int chromNum_sd = indexInfo_sd->returnChromNum();
	int chromNum_ps = indexInfo_ps->returnChromNum();
	if(chromNum_sd != chromNum_ps)
	{
		cout << "different chrom # for sd and ps genome" << endl;
		exit(1);
	}
	char *chrom_sd; chrom_sd = (char*)malloc((indexInfo_sd->returnIndexSize()) * sizeof(char));
	char *chrom_ps; chrom_ps = (char*)malloc((indexInfo_ps->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs_sd.read((char*)chrom_sd, (indexInfo_sd->returnIndexSize()) * sizeof(char));
	chrom_bit_file_ifs_ps.read((char*)chrom_ps, (indexInfo_ps->returnIndexSize()) * sizeof(char));
	indexInfo_sd->readGenome(chrom_sd);
	indexInfo_ps->readGenome(chrom_ps);
	indexInfo_sd->initiate();
	indexInfo_ps->initiate();
	indexInfo_sd->initiateChrNameIndexArray(1000);
	indexInfo_ps->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	cout << "start to load SNP file " << endl;
	string inputFormattedSNPfile = argv[3];
	SNPhash_Info tmpSNPhashInfo;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum_sd);
	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	tmpSNPhashInfo.generateSNPhash_formattedSNPfile(inputFormattedSNPfile, indexInfo_sd);

	cout << "start to generate gene boundary and corresponding splice sites" << endl; 
	string gtfFile = argv[4];
	string SJfile = argv[5];
	vector<string> geneNameVec;
	vector<string> geneStrandVec;
	vector<int> geneChrNameIntVec;
	vector< pair<int,int> > geneBoundaryPosPairVec;
	vector< vector<int> > geneSpliceDonerSiteVecVec;
	vector< vector<int> > geneSpliceAcceptorSiteVecVec;
	generateGeneSpliceSiteVecFromGTF(gtfFile, SJfile, geneNameVec, geneStrandVec, geneChrNameIntVec, 
		geneBoundaryPosPairVec, geneSpliceDonerSiteVecVec, geneSpliceAcceptorSiteVecVec, indexInfo_sd);

	cout << "start to locate SNPs in each gene and derive candidiate personal SJs for each gene" << endl;
	int geneNum = geneNameVec.size();
	cout << "geneNum: " << geneNum << endl;
	cout << "start to locate SNPs in each gene and derive candidiate personal splice sites for each gene " << endl;
	vector< vector<int> > genePsSpliceDonerSiteVecVec;
	vector< vector<int> > genePsSpliceAcceptorSiteVecVec;
	for(int tmp = 0; tmp < geneNum; tmp++)
	{
		vector<int> tmpGenePsSpliceDonerSiteVec;
		vector<int> tmpGenePsSpliceAcceptorSiteVec;
		int existingSpliceSiteNum = geneSpliceDonerSiteVecVec[tmp].size();
		if(existingSpliceSiteNum == 0)
		{
			genePsSpliceDonerSiteVecVec.push_back(tmpGenePsSpliceDonerSiteVec);
			genePsSpliceAcceptorSiteVecVec.push_back(tmpGenePsSpliceAcceptorSiteVec);
			continue;
		}
		int tmpGeneChrNameInt = geneChrNameIntVec[tmp];
		int tmpGeneStartPos = geneBoundaryPosPairVec[tmp].first;
		int tmpGeneEndPos = geneBoundaryPosPairVec[tmp].second;
		string tmpGeneStrand = geneStrandVec[tmp];
		vector<int> tmpSNPposVecInGene;
		tmpSNPhashInfo.returnSNPposVecWithinRegion(tmpGeneChrNameInt, tmpGeneStartPos, 
			tmpGeneEndPos, tmpSNPposVecInGene);
		int tmpSNPnum = tmpSNPposVecInGene.size();
		//cout << "tmpSNPnum: " << tmpSNPnum << endl;
		if(tmpSNPnum == 0)
		{
			genePsSpliceDonerSiteVecVec.push_back(tmpGenePsSpliceDonerSiteVec);
			genePsSpliceAcceptorSiteVecVec.push_back(tmpGenePsSpliceAcceptorSiteVec);
			continue;
		}
		generatePsSJsiteVec(tmpGenePsSpliceDonerSiteVec, tmpGenePsSpliceAcceptorSiteVec,
			tmpGeneChrNameInt, tmpGeneStartPos, tmpGeneEndPos,
			geneSpliceDonerSiteVecVec[tmp], geneSpliceAcceptorSiteVecVec[tmp], 
			tmpSNPposVecInGene, tmpGeneStrand, indexInfo_ps);
		genePsSpliceDonerSiteVecVec.push_back(tmpGenePsSpliceDonerSiteVec);
		genePsSpliceAcceptorSiteVecVec.push_back(tmpGenePsSpliceAcceptorSiteVec);
	}

	cout << "start to output personal SJs " << endl;
	string outputFilePrefix = argv[6];
	/*string outputCandiPsSJ = outputFilePrefix + "_PsSJ.txt";
	ofstream candiPsSJ_ofs(outputCandiPsSJ.c_str());
	for(int tmp = 0; tmp < geneNum; tmp++)
	{
		int tmpGenePsSJnum = genePsSpliceSitePairVecVec[tmp].size();
		if(tmpGenePsSJnum == 0)
			continue;
		string tmpGeneName = geneNameVec[tmp];
		int tmpGeneChrNameInt = geneChrNameIntVec[tmp];
		string tmpGeneStrand = geneStrandVec[tmp];
		string tmpGeneChrName = indexInfo_ps->returnChrNameStr(tmpGeneChrNameInt);
		for(int tmpPsSJindex = 0; tmpPsSJindex < tmpGenePsSJnum; tmpPsSJindex ++)
		{
			int tmpPsSJ_donerSite = (genePsSpliceSitePairVecVec[tmp])[tmpPsSJindex].first;
			int tmpPsSJ_acceptorSite = (genePsSpliceSitePairVecVec[tmp])[tmpPsSJindex].second;
			candiPsSJ_ofs << tmpGeneChrName << "\t" << tmpPsSJ_donerSite << "\t"
				<< tmpPsSJ_acceptorSite << "\t" << tmpGeneStrand << "\t" << tmpGeneName << endl;
		}
	}
	candiPsSJ_ofs.close();*/
	string outputCandiPsSpliceSite = outputFilePrefix + "_PsSpliceSite.txt";
	ofstream candiPsSpliceSite_ofs(outputCandiPsSpliceSite.c_str());
	string outputCandiPsSpliceSite_doner = outputFilePrefix + "_PsSpliceSite_doner.txt";
	ofstream candiPsSpliceSite_doner_ofs(outputCandiPsSpliceSite_doner.c_str());
	string outputCandiPsSpliceSite_acceptor = outputFilePrefix + "_PsSpliceSite_acceptor.txt";
	ofstream candiPsSpliceSite_acceptor_ofs(outputCandiPsSpliceSite_acceptor.c_str());
	int psDonerSite_total = 0;
	int psAcceptorSite_total = 0;
	int annotatedDonerSite_total = 0;
	int annotatedAcceptorSite_total = 0;
	for(int tmp = 0; tmp < geneNum; tmp++)
	{	
		int tmpGenePsSpliceDonerSiteNum = genePsSpliceDonerSiteVecVec[tmp].size();
		int tmpGenePsSpliceAcceptorSiteNum = genePsSpliceAcceptorSiteVecVec[tmp].size();
		if(tmpGenePsSpliceDonerSiteNum + tmpGenePsSpliceAcceptorSiteNum == 0)
			continue;
		psDonerSite_total += tmpGenePsSpliceDonerSiteNum;
		psAcceptorSite_total += tmpGenePsSpliceAcceptorSiteNum;
		string tmpGeneName = geneNameVec[tmp];
		int tmpGeneChrNameInt = geneChrNameIntVec[tmp];
		string tmpGeneChrName = indexInfo_ps->returnChrNameStr(tmpGeneChrNameInt);
		int tmpGeneStartPos = geneBoundaryPosPairVec[tmp].first;
		int tmpGeneEndPos = geneBoundaryPosPairVec[tmp].second;
		string tmpGeneStrand = geneStrandVec[tmp];
		int tmpGeneSpliceDonerSiteNum = geneSpliceDonerSiteVecVec[tmp].size();
		int tmpGeneSpliceAcceptorSiteNum = geneSpliceAcceptorSiteVecVec[tmp].size();
		annotatedDonerSite_total += tmpGeneSpliceDonerSiteNum;
		annotatedAcceptorSite_total += tmpGeneSpliceAcceptorSiteNum;
		candiPsSpliceSite_ofs << tmpGeneChrName << "\t" << tmpGeneStartPos << "\t"
			<< tmpGeneEndPos << "\t" << tmpGeneStrand << "\t" << tmpGeneName << "\t"
			<< tmpGeneSpliceDonerSiteNum << "\t" << tmpGeneSpliceAcceptorSiteNum << "\t"
			<< tmpGenePsSpliceDonerSiteNum << "\t" << tmpGenePsSpliceAcceptorSiteNum << endl;
		for(int tmp2 = 0; tmp2 < tmpGeneSpliceDonerSiteNum; tmp2 ++)
			candiPsSpliceSite_ofs << (geneSpliceDonerSiteVecVec[tmp])[tmp2] << ",";
		candiPsSpliceSite_ofs << endl;
		for(int tmp2 = 0; tmp2 < tmpGeneSpliceAcceptorSiteNum; tmp2 ++)
			candiPsSpliceSite_ofs << (geneSpliceAcceptorSiteVecVec[tmp])[tmp2] << ",";
		candiPsSpliceSite_ofs << endl;
		if(tmpGenePsSpliceDonerSiteNum == 0)
			candiPsSpliceSite_ofs << "NULL";
		else
		{	
			for(int tmp2 = 0; tmp2 < tmpGenePsSpliceDonerSiteNum; tmp2 ++)
			{
				candiPsSpliceSite_ofs << (genePsSpliceDonerSiteVecVec[tmp])[tmp2] << ",";
				candiPsSpliceSite_doner_ofs << tmpGeneChrName << "\t" << (genePsSpliceDonerSiteVecVec[tmp])[tmp2] << endl;
			}
		}
		candiPsSpliceSite_ofs << endl;
		if(tmpGenePsSpliceAcceptorSiteNum == 0)
			candiPsSpliceSite_ofs << "NULL";
		else
		{	
			for(int tmp2 = 0; tmp2 < tmpGenePsSpliceAcceptorSiteNum; tmp2 ++)
			{
				candiPsSpliceSite_ofs << (genePsSpliceAcceptorSiteVecVec[tmp])[tmp2] << ",";
				candiPsSpliceSite_acceptor_ofs << tmpGeneChrName << "\t" << (genePsSpliceAcceptorSiteVecVec[tmp])[tmp2] << endl;
			}
		}
		candiPsSpliceSite_ofs << endl;
	}
	candiPsSpliceSite_ofs.close();
	cout << "all jobs done" << endl;
	cout << endl << "psDonerSite_total: " << psDonerSite_total << endl;
	cout << "psAcceptorSite_total: " << psAcceptorSite_total << endl;
	cout << endl << "annotatedDonerSite_total: " << annotatedDonerSite_total << endl;
	cout << "annotatedAcceptorSite_total: " << annotatedAcceptorSite_total << endl;

	delete indexInfo_sd;
	delete indexInfo_ps;
	free(chrom_sd);
	free(chrom_ps);
	chrom_bit_file_ifs_sd.close();
	chrom_bit_file_ifs_ps.close();
	parameter_ifs_sd.close();
	parameter_ifs_ps.close();
	candiPsSpliceSite_ofs.close();
	candiPsSpliceSite_doner_ofs.close();
	candiPsSpliceSite_acceptor_ofs.close();
	return 0;
}