// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ENHANCED_SUFFIX_ARRAY_INFO_H
#define ENHANCED_SUFFIX_ARRAY_INFO_H

#include <string>
#include <vector>
#include <set>
#include <malloc.h>

using namespace std;


class Enhanced_Suffix_Array_Info
{
private:
	unsigned long long index_length;
	unsigned long long child_SA_size;
	unsigned long long verifyChild_lcpCompress_size;
	unsigned int* child_SA;
	BYTE* verifyChild_lcpCompress;

public:
	Enhanced_Suffix_Array_Info(unsigned int index_size)
	{
		index_length = (unsigned long long)index_size;
		cout << "index_length: " << index_length << endl;
		child_SA_size = (unsigned long long)(index_length * 2);
		verifyChild_lcpCompress_size = child_SA_size;
		cout << "child_SA_size: " << child_SA_size << endl; 
		cout << "verifyChild_lcpCompress_size: " << verifyChild_lcpCompress_size << endl;
		child_SA = (unsigned int*)malloc((child_SA_size) * sizeof(unsigned int));
		verifyChild_lcpCompress = (BYTE*)malloc(verifyChild_lcpCompress_size * sizeof(BYTE));
	}

	void freeMem()
	{
		free(child_SA);
		free(verifyChild_lcpCompress);
	}

	void loadAllIndexFiles(const string& indexStr, unsigned int index_size, ofstream& log_ofs)
	{
		string SA_file = indexStr; SA_file.append("_SA"); ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
		string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
		string childTab_file = indexStr; childTab_file.append("_childTab"); ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
		string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);
		string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
		string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);	

		cout << "index_size: " << index_size << endl;

		log_ofs << "start to load SA" << endl;
	    unsigned int *sa; sa = (unsigned int*)malloc(index_size * sizeof(unsigned int)); SA_file_ifs.read((char*)sa, index_size * sizeof(unsigned int));
		log_ofs << "start to load lcpCompress" << endl;
		BYTE *lcpCompress; lcpCompress = (BYTE*)malloc(index_size * sizeof(BYTE)); lcpCompress_file_ifs.read((char*)lcpCompress, index_size * sizeof(BYTE));	
		log_ofs << "start to load childTab " << endl;
		unsigned int *childTab; childTab = (unsigned int*)malloc(index_size * sizeof(unsigned int)); childTab_file_ifs.read((char*)childTab, index_size * sizeof(unsigned int));
		log_ofs << "start to load detChild" << endl;
		BYTE *verifyChild; verifyChild = (BYTE*)malloc(index_size * sizeof(BYTE)); verifyChild_file_ifs.read((char*)verifyChild, index_size * sizeof(BYTE));
		log_ofs << "All index files loaded" << endl;
	
		this->initiate_child_SA(childTab, sa);
		free(childTab);
		free(sa);

		this->initiate_verifyChild_lcpCompress(verifyChild, lcpCompress);

		free(lcpCompress);//free(child_up);free(child_down);free(child_next);		
		free(verifyChild);	
	}



	void initiate_child_SA(unsigned int* child, unsigned int* SA)
	{		
		cout << "start to initiate child_SA " << endl;
		for(unsigned long long tmp = 0; tmp < index_length; tmp++)
		{
			unsigned long long tmp_index = tmp * 2;
			//cout << "tmp: " << tmp << endl;
			child_SA[tmp_index] = child[tmp];
			child_SA[tmp_index + 1] = SA[tmp];
		}
		cout << "finish initiating child_SA" << endl;
	}

	void initiate_verifyChild_lcpCompress(BYTE* verifyChild, BYTE* lcpCompress)
	{
		cout << "start to initiate verifyChild_lcpCompress " << endl;
		for(unsigned long long tmp = 0; tmp < index_length; tmp++)
		{
			unsigned long long tmp_index = tmp * 2;
			//cout << "tmp: " << tmp << endl;
			verifyChild_lcpCompress[tmp_index] = verifyChild[tmp];
			verifyChild_lcpCompress[tmp_index + 1] = lcpCompress[tmp];
		}		
		cout << "finish initiating verifyChild_lcpCompress " << endl;
	}

	unsigned int returnChild(unsigned int index_child)
	{
		unsigned long long tmp_index = (unsigned long long)index_child * 2;
		return child_SA[tmp_index];
	}
	unsigned int returnSA(unsigned int index_SA)
	{
		unsigned long long tmp_index = (unsigned long long)index_SA * 2 + 1;
		return child_SA[tmp_index];
	}

	unsigned int returnVerifyChild(unsigned int index_verifyChild)
	{
		unsigned long long tmp_index = (unsigned long long)index_verifyChild * 2;
		return verifyChild_lcpCompress[tmp_index];
	}
	unsigned int returnLcpCompress(unsigned int index_lcpCompress)
	{
		unsigned long long tmp_index = (unsigned long long)index_lcpCompress * 2 + 1;
		return verifyChild_lcpCompress[tmp_index];
	}

};

/*
class Enhanced_Suffix_Array_Info
{
private:
	unsigned long long index_length;
	unsigned long long ESA_size;
	unsigned int* ESA;
	//BYTE* verifyChild_lcpCompress;

public:
	Enhanced_Suffix_Array_Info(unsigned int index_size)
	{
		index_length = (unsigned long long)index_size;
		cout << "index_length: " << index_length << endl;
		ESA_size = (unsigned long long)(index_length * 4);
		cout << "ESA_size: " << ESA_size << endl; 
		ESA = (unsigned int*)malloc((ESA_size) * sizeof(unsigned int));
	}

	void loadAllIndexFiles(const string& indexStr, unsigned int index_size, ofstream& log_ofs)
	{
		string SA_file = indexStr; SA_file.append("_SA"); ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
		string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
		string childTab_file = indexStr; childTab_file.append("_childTab"); ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
		string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);
		string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
		string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);	

		cout << "index_size: " << index_size << endl;

		log_ofs << "start to load SA" << endl;
	    unsigned int *sa; sa = (unsigned int*)malloc(index_size * sizeof(unsigned int)); SA_file_ifs.read((char*)sa, index_size * sizeof(unsigned int));
		log_ofs << "start to load lcpCompress" << endl;
		BYTE *lcpCompress; lcpCompress = (BYTE*)malloc(index_size * sizeof(BYTE)); lcpCompress_file_ifs.read((char*)lcpCompress, index_size * sizeof(BYTE));	
		log_ofs << "start to load childTab " << endl;
		unsigned int *childTab; childTab = (unsigned int*)malloc(index_size * sizeof(unsigned int)); childTab_file_ifs.read((char*)childTab, index_size * sizeof(unsigned int));
		log_ofs << "start to load detChild" << endl;
		BYTE *verifyChild; verifyChild = (BYTE*)malloc(index_size * sizeof(BYTE)); verifyChild_file_ifs.read((char*)verifyChild, index_size * sizeof(BYTE));
		log_ofs << "All index files loaded" << endl;
	
		this->initiate_ESA(childTab, sa, verifyChild, lcpCompress);

		free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
		free(childTab);
		free(verifyChild);	
	}



	void initiate_ESA(unsigned int* child, unsigned int* SA, BYTE* verifyChild, BYTE* lcpCompress)
	{		
		cout << "start to initiate ESA " << endl;
		for(unsigned long long tmp = 0; tmp < index_length; tmp++)
		{
			unsigned long long tmp_index = tmp * 4;
			//cout << "tmp: " << tmp << endl;
			ESA[tmp_index] = child[tmp];
			ESA[tmp_index + 1] = SA[tmp];
			ESA[tmp_index + 2] = verifyChild[tmp];
			ESA[tmp_index + 3] = lcpCompress[tmp];
		}
		cout << "finish initiating ESA" << endl;
	}

	unsigned int returnChild(unsigned int index_child)
	{
		unsigned long long tmp_index = (unsigned long long)index_child * 4;
		return ESA[tmp_index];
	}
	unsigned int returnSA(unsigned int index_SA)
	{
		unsigned long long tmp_index = (unsigned long long)index_SA * 4 + 1;
		return ESA[tmp_index];
	}

	unsigned int returnVerifyChild(unsigned int index_verifyChild)
	{
		unsigned long long tmp_index = (unsigned long long)index_verifyChild * 4 + 2;
		return ESA[tmp_index];
	}
	unsigned int returnLcpCompress(unsigned int index_lcpCompress)
	{
		unsigned long long tmp_index = (unsigned long long)index_lcpCompress * 4 + 3;
		return ESA[tmp_index];
	}

};*/

#endif