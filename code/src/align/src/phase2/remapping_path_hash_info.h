// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef REMAPPING_PATH_HASH_INFO_H
#define REMAPPING_PATH_HASH_INFO_H

#include <string>
#include <string.h>
#include <map>
#include <set>

using namespace std;

typedef map<string, vector< Remapping_Path_Info > > Remapping_Path_Info_Map;
typedef map<int, Remapping_Path_Info_Map> Remapping_Path_Hash_Info_Map;

class Remapping_Path_Hash_Info
{
private:
	vector<Remapping_Path_Hash_Info_Map> mapPathInfoHashVec_forw;
	vector<Remapping_Path_Hash_Info_Map> mapPathInfoHashVec_reve;

	int anchorStringLength;

	int areaSize;

	vector<SJareaHash> SJstartPosAreaHash;
	vector<SJareaHash> SJendPosAreaHash;
public:
	Remapping_Path_Hash_Info()
	{
		areaSize = 1000;
		anchorStringLength = 3;		
	}


};
#endif