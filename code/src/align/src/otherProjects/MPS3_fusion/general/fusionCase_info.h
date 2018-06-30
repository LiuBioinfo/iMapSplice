// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FUSIONCASE_INFO_H
#define FUSIONCASE_INFO_H

#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"

using namespace std;

class FusionCase_Info
{
private:
	vector<string> canonicalFusionJuncFlankStringVec;
	vector<string> fusionStrandVec;
	vector<string> fusionMapDirStrVec;
public:
	FusionCase_Info()
	{}

	void initiate()
	{
		canonicalFusionJuncFlankStringVec.push_back("GTAG"); // case 1
		canonicalFusionJuncFlankStringVec.push_back("GTAG"); // case 2
		canonicalFusionJuncFlankStringVec.push_back("NULL"); // case 3
		canonicalFusionJuncFlankStringVec.push_back("CTAC"); // case 4
		canonicalFusionJuncFlankStringVec.push_back("CTAC"); // case 5
		canonicalFusionJuncFlankStringVec.push_back("NULL"); // case 6
		canonicalFusionJuncFlankStringVec.push_back("CTGT"); // case 7
		canonicalFusionJuncFlankStringVec.push_back("GTCT"); // case 8
		canonicalFusionJuncFlankStringVec.push_back("NULL"); // case 9
		canonicalFusionJuncFlankStringVec.push_back("ACAG"); // case 10
		canonicalFusionJuncFlankStringVec.push_back("AGAC"); // case 11
		canonicalFusionJuncFlankStringVec.push_back("NULL"); // case 12

		fusionStrandVec.push_back("++"); // case 1
		fusionStrandVec.push_back("++"); // case 2
		fusionStrandVec.push_back("++"); // case 3
		fusionStrandVec.push_back("--"); // case 4
		fusionStrandVec.push_back("--"); // case 5
		fusionStrandVec.push_back("--"); // case 6
		fusionStrandVec.push_back("+-"); // case 7
		fusionStrandVec.push_back("+-"); // case 8
		fusionStrandVec.push_back("+-"); // case 9
		fusionStrandVec.push_back("-+"); // case 10
		fusionStrandVec.push_back("-+"); // case 11
		fusionStrandVec.push_back("-+"); // case 12

		fusionMapDirStrVec.push_back("FOR_1~REV_1+REV_2"); // case 1
		fusionMapDirStrVec.push_back("FOR_1+FOR_2~REV_2"); // case 2
		fusionMapDirStrVec.push_back("FOR_1~REV_2"); // case 3
		fusionMapDirStrVec.push_back("FOR_2~REV_2+REV_1"); // case 4
		fusionMapDirStrVec.push_back("FOR_2+FOR_1~REV_1"); // case 5
		fusionMapDirStrVec.push_back("FOR_2~REV_1"); // case 6
		fusionMapDirStrVec.push_back("FOR_2~REV_2+FOR_1"); // case 7
		fusionMapDirStrVec.push_back("FOR_1~REV_1+FOR_2"); // case 8
		fusionMapDirStrVec.push_back("FOR_1~FOR_2"); // case 9
		fusionMapDirStrVec.push_back("REV_1+FOR_2~REV_2"); // case 10
		fusionMapDirStrVec.push_back("REV_2+FOR_1~REV_1"); // case 11
		fusionMapDirStrVec.push_back("REV_1~REV_2"); // case 12
	}

	string returnCanonicalFlankString_fusionCase(int tmpCase)
	{
		int tmpCaseIndex = tmpCase-1;
		return canonicalFusionJuncFlankStringVec[tmpCaseIndex];
	}

	string returnFusionStrand_fusionCase(int tmpCase)
	{
		int tmpCaseIndex = tmpCase-1;
		return fusionStrandVec[tmpCaseIndex];
	}

	string returnFusionMapDir_fusioncase(int tmpCase)
	{
		int tmpCaseIndex = tmpCase-1;
		return fusionMapDirStrVec[tmpCaseIndex];
	}
};
#endif