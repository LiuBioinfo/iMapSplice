// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef STATS_INFO_H
#define STATS_INFO_H

#include <string>
#include <string.h>
#include <vector>

using namespace std;

class Stats_Info
{
private:
	// *******************************  SE -- phase1 *************************************  //
	vector<unsigned int> SE_num_complete_unmapped_lowScore_unique_phase1_thread;
	vector<unsigned int> SE_num_complete_unmapped_lowScore_multi_phase1_thread;

	vector<unsigned int> SE_num_complete_unique_phase1_thread;
	vector<unsigned int> SE_num_complete_multi_phase1_thread;
	vector<unsigned int> SE_num_incomplete_unique_phase1_thread;
	vector<unsigned int> SE_num_incomplete_multi_phase1_thread;

	vector<unsigned int> SE_num_unmapped_phase1_thread;
	vector<unsigned int> SE_num_unmapped_mappedToRepeatRegion_phase1_thread;

	unsigned int SE_num_complete_unmapped_lowScore_unique_phase1;
	unsigned int SE_num_complete_unmapped_lowScore_multi_phase1;

	unsigned int SE_num_complete_unique_phase1;
	unsigned int SE_num_complete_multi_phase1;
	unsigned int SE_num_incomplete_unique_phase1;
	unsigned int SE_num_incomplete_multi_phase1;

	unsigned int SE_num_unmapped_phase1;
	unsigned int SE_num_unmapped_mappedToRepeatRegion_phase1;

	// *******************************  PE  -- phase1 *************************************  //
	vector<unsigned int> PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1_thread;
	unsigned int PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1;
	vector<unsigned int> PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1_thread;
	unsigned int PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1;

	vector<unsigned int> PE_pair_num_paired_complete_unique_phase1_thread;
	vector<unsigned int> PE_pair_num_paired_complete_multi_phase1_thread;	
	vector<unsigned int> PE_pair_num_paired_incomplete_unique_phase1_thread;
	vector<unsigned int> PE_pair_num_paired_incomplete_multi_phase1_thread;

	vector<unsigned int> PE_pair_num_unmapped_phase1_thread;
	vector<unsigned int> PE_pair_num_unmapped_mappedToRepeatRegion_phase1_thread;
	vector<unsigned int> PE_pair_num_unpaired_phase1_thread;

	unsigned int PE_pair_num_paired_complete_unique_phase1;
	unsigned int PE_pair_num_paired_complete_multi_phase1;
	unsigned int PE_pair_num_paired_incomplete_unique_phase1;
	unsigned int PE_pair_num_paired_incomplete_multi_phase1;

	unsigned int PE_pair_num_unmapped_phase1;
	unsigned int PE_pair_num_unmapped_mappedToRepeatRegion_phase1;
	unsigned int PE_pair_num_unpaired_phase1; 

	// ************  SE  -- unpaired -- fixUnpaired is not needed for SE reads *********** //
	// ********************************  PE  -- unpaired ****************************** //
	vector<unsigned int> PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair_thread;
	unsigned int PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair;
	vector<unsigned int> PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair_thread;
	unsigned int PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair;

	vector<unsigned int> PE_pair_num_paired_complete_unique_fixUnpaired_thread;
	vector<unsigned int> PE_pair_num_paired_complete_multi_fixUnpaired_thread;
	vector<unsigned int> PE_pair_num_paired_incomplete_unique_fixUnpaired_thread;
	vector<unsigned int> PE_pair_num_paired_incomplete_multi_fixUnpaired_thread;

	vector<unsigned int> PE_pair_num_unpaired_complete_fixUnpaired_thread;
	vector<unsigned int> PE_pair_num_unpaired_incomplete_fixUnpaired_thread;	

	unsigned int PE_pair_num_paired_complete_unique_fixUnpaired;
	unsigned int PE_pair_num_paired_complete_multi_fixUnpaired;
	unsigned int PE_pair_num_paired_incomplete_unique_fixUnpaired;
	unsigned int PE_pair_num_paired_incomplete_multi_fixUnpaired;

	unsigned int PE_pair_num_unpaired_complete_fixUnpaired;
	unsigned int PE_pair_num_unpaired_incomplete_fixUnpaired;	

	// *****************************  SE -- fixHeadTail  ********************************** //
	vector<unsigned int> SE_num_complete_unmapped_lowScore_unique_fixHeadTail_thread;
	vector<unsigned int> SE_num_complete_unmapped_lowScore_multi_fixHeadTail_thread; 
	vector<unsigned int> SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail_thread;
	vector<unsigned int> SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail_thread; 

	unsigned int SE_num_complete_unmapped_lowScore_unique_fixHeadTail;
	unsigned int SE_num_complete_unmapped_lowScore_multi_fixHeadTail;
	unsigned int SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail;
	unsigned int SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail;

	vector<unsigned int> SE_num_complete_unique_fixHeadTail_thread;
	vector<unsigned int> SE_num_complete_multi_fixHeadTail_thread;
	vector<unsigned int> SE_num_incomplete_unique_fixHeadTail_thread;
	vector<unsigned int> SE_num_incomplete_multi_fixHeadTail_thread;

	unsigned int SE_num_complete_unique_fixHeadTail;
	unsigned int SE_num_complete_multi_fixHeadTail;
	unsigned int SE_num_incomplete_unique_fixHeadTail;
	unsigned int SE_num_incomplete_multi_fixHeadTail;

	// *****************************  PE -- fixHeadTail  ********************************** //
	vector<unsigned int> PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail_thread;
	unsigned int PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail;
	vector<unsigned int> PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail_thread;
	unsigned int PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail;

	vector<unsigned int> PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail_thread;
	unsigned int PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail;
	vector<unsigned int> PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail_thread;
	unsigned int PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail;

	vector<unsigned int> PE_pair_num_paired_complete_unique_fixHeadTail_thread;
	vector<unsigned int> PE_pair_num_paired_complete_multi_fixHeadTail_thread;
	vector<unsigned int> PE_pair_num_paired_incomplete_unique_fixHeadTail_thread;
	vector<unsigned int> PE_pair_num_paired_incomplete_multi_fixHeadTail_thread;

	vector<unsigned int> PE_pair_num_unpaired_complete_fixHeadTail_thread;
	vector<unsigned int> PE_pair_num_unpaired_incomplete_fixHeadTail_thread;

	unsigned int PE_pair_num_paired_complete_unique_fixHeadTail;
	unsigned int PE_pair_num_paired_complete_multi_fixHeadTail;
	unsigned int PE_pair_num_paired_incomplete_unique_fixHeadTail;
	unsigned int PE_pair_num_paired_incomplete_multi_fixHeadTail;

	unsigned int PE_pair_num_unpaired_complete_fixHeadTail;
	unsigned int PE_pair_num_unpaired_incomplete_fixHeadTail;

public:
	Stats_Info()
	{
		// *************************   SE   ************************  //
		SE_num_complete_unmapped_lowScore_unique_phase1 = 0;
		SE_num_complete_unmapped_lowScore_multi_phase1 = 0;

		SE_num_complete_unique_phase1 = 0;
		SE_num_complete_multi_phase1 = 0;
		SE_num_incomplete_unique_phase1 = 0;
		SE_num_incomplete_multi_phase1 = 0;

		SE_num_unmapped_phase1 = 0;
		SE_num_unmapped_mappedToRepeatRegion_phase1 = 0;		

		SE_num_complete_unmapped_lowScore_unique_fixHeadTail = 0;
		SE_num_complete_unmapped_lowScore_multi_fixHeadTail = 0;
		SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail = 0;
		SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail = 0;

		SE_num_complete_unique_fixHeadTail = 0;
		SE_num_complete_multi_fixHeadTail = 0;
		SE_num_incomplete_unique_fixHeadTail = 0;
		SE_num_incomplete_multi_fixHeadTail = 0;
			
		// *************************   PE   ************************  //
		PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1 = 0;
		PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1 = 0;

		PE_pair_num_paired_complete_unique_phase1 = 0;
		PE_pair_num_paired_complete_multi_phase1 = 0;	
		PE_pair_num_paired_incomplete_unique_phase1 = 0;
		PE_pair_num_paired_incomplete_multi_phase1 = 0;

		PE_pair_num_unmapped_phase1 = 0;	
		PE_pair_num_unmapped_mappedToRepeatRegion_phase1 = 0;
		PE_pair_num_unpaired_phase1 = 0; 

		PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair = 0;
		PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair = 0;

		PE_pair_num_paired_complete_unique_fixUnpaired = 0;
		PE_pair_num_paired_complete_multi_fixUnpaired = 0;
		PE_pair_num_paired_incomplete_unique_fixUnpaired = 0;
		PE_pair_num_paired_incomplete_multi_fixUnpaired = 0;

		PE_pair_num_unpaired_complete_fixUnpaired = 0;
		PE_pair_num_unpaired_incomplete_fixUnpaired = 0;	

		PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail = 0;
		PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail = 0;
		PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail = 0;
		PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail = 0;

		PE_pair_num_paired_complete_unique_fixHeadTail = 0;
		PE_pair_num_paired_complete_multi_fixHeadTail = 0;
		PE_pair_num_paired_incomplete_unique_fixHeadTail = 0;
		PE_pair_num_paired_incomplete_multi_fixHeadTail = 0;

		PE_pair_num_unpaired_complete_fixHeadTail = 0;
		PE_pair_num_unpaired_incomplete_fixHeadTail = 0;
	}

	void outputAllStats(ofstream& output, bool Do_phase1_only_or_not_bool)
	{
		unsigned int read_num;
	
		unsigned int bothEndsUnmapped_repeatRegion_num = PE_pair_num_unmapped_mappedToRepeatRegion_phase1;
		unsigned int bothEndsUnmapped_num = PE_pair_num_unmapped_phase1;
		unsigned int unmapped_dueToLowScore_num 
			= PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1
				+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair
				+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail
				+ PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail
				+ PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail;

		if(Do_phase1_only_or_not_bool)
		{
			read_num = PE_pair_num_paired_complete_unique_phase1 + PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_multi_phase1 + PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_unpaired_phase1 + PE_pair_num_unmapped_phase1
				+ bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		}
		else
		{
			read_num = 	PE_pair_num_paired_complete_unique_phase1 //+ PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_unique_fixUnpaired //+ PE_pair_num_paired_incomplete_unique_fixUnpaired
				+ PE_pair_num_paired_complete_unique_fixHeadTail + PE_pair_num_paired_incomplete_unique_fixHeadTail;
				+ PE_pair_num_paired_complete_multi_phase1 //+ PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_paired_complete_multi_fixUnpaired //+ PE_pair_num_paired_incomplete_multi_fixUnpaired
				+ PE_pair_num_paired_complete_multi_fixHeadTail + PE_pair_num_paired_incomplete_multi_fixHeadTail;
				+ PE_pair_num_unpaired_complete_fixUnpaired + PE_pair_num_unpaired_incomplete_fixUnpaired
				+ PE_pair_num_unpaired_complete_fixHeadTail + PE_pair_num_unpaired_incomplete_fixHeadTail
				+ bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		}

		if(Do_phase1_only_or_not_bool)
		{
			read_num = PE_pair_num_paired_complete_unique_phase1 + PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_multi_phase1 + PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_unpaired_phase1 + PE_pair_num_unmapped_phase1
				+ bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		}
		else
		{
			read_num = 	PE_pair_num_paired_complete_unique_phase1 //+ PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_unique_fixUnpaired //+ PE_pair_num_paired_incomplete_unique_fixUnpaired
				+ PE_pair_num_paired_complete_unique_fixHeadTail + PE_pair_num_paired_incomplete_unique_fixHeadTail;
				+ PE_pair_num_paired_complete_multi_phase1 //+ PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_paired_complete_multi_fixUnpaired //+ PE_pair_num_paired_incomplete_multi_fixUnpaired
				+ PE_pair_num_paired_complete_multi_fixHeadTail + PE_pair_num_paired_incomplete_multi_fixHeadTail;
				+ PE_pair_num_unpaired_complete_fixUnpaired + PE_pair_num_unpaired_incomplete_fixUnpaired
				+ PE_pair_num_unpaired_complete_fixHeadTail + PE_pair_num_unpaired_incomplete_fixHeadTail
				+ bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		}

		output << endl << "total reads #: " << read_num << endl << endl;
		output << endl << "phase1 statics: " << endl << endl 
			<< "PE_pair_num_paired_complete_unique_phase1: " //<< endl
			<< PE_pair_num_paired_complete_unique_phase1 << endl;
		output << "PE_pair_num_paired_complete_multi_phase1: " //<< endl
			<< PE_pair_num_paired_complete_multi_phase1 << endl;
		output << "PE_pair_num_paired_incomplete_unique_phase1: " //<< endl
			<< PE_pair_num_paired_incomplete_unique_phase1 << endl;
		output << "PE_pair_num_paired_incomplete_multi_phase1: " //<< endl
			<< PE_pair_num_paired_incomplete_multi_phase1 << endl;	

		output << "PE_pair_num_unpaired_phase1: " //<< endl
			<< PE_pair_num_unpaired_phase1 << endl;
		output << "PE_pair_num_unmapped_mappedToRepeatRegion_phase1: " //<< endl
			<< PE_pair_num_unmapped_mappedToRepeatRegion_phase1 << endl;
		output << "PE_pair_num_unmapped_phase1: " //<< endl
			<< PE_pair_num_unmapped_phase1 << endl << endl;

		output << "PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1 << endl;
		output << "PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1 << endl;	


		output << endl << endl << "fix unpaired reads output: " << endl << endl;
		output << "PE_pair_num_paired_complete_unique_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_complete_unique_fixUnpaired << endl;
		output << "PE_pair_num_paired_complete_multi_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_complete_multi_fixUnpaired << endl;
		output << "PE_pair_num_paired_incomplete_unique_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_incomplete_unique_fixUnpaired << endl;
		output << "PE_pair_num_paired_incomplete_multi_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_incomplete_multi_fixUnpaired << endl;			

		output << "PE_pair_num_unpaired_complete_fixUnpaired: " //<< endl
			<< PE_pair_num_unpaired_complete_fixUnpaired << endl; 
		output << "PE_pair_num_unpaired_incomplete_fixUnpaired: " //<< endl
			<< PE_pair_num_unpaired_incomplete_fixUnpaired << endl << endl; 

		output << "PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair: " 
			<< PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair << endl;
		output << "PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair: " 
			<< PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair << endl;


		output << endl << endl << "fix headTail output: " << endl << endl;
		output << "PE_pair_num_paired_complete_unique_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_complete_unique_fixHeadTail << endl;
		output << "PE_pair_num_paired_complete_multi_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_complete_multi_fixHeadTail << endl;
		output << "PE_pair_num_paired_incomplete_unique_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_incomplete_unique_fixHeadTail << endl;
		output << "PE_pair_num_paired_incomplete_multi_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_incomplete_multi_fixHeadTail << endl;	

		output << "PE_pair_num_unpaired_complete_fixHeadTail: " //<< endl
			<< PE_pair_num_unpaired_complete_fixHeadTail << endl; 
		output << "PE_pair_num_unpaired_incomplete_fixHeadTail: " //<< endl
			<< PE_pair_num_unpaired_incomplete_fixHeadTail << endl << endl; 

		output << "PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail << endl;
		output << "PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail << endl;			
		output << "PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail: "
			<< PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail << endl;
		output << "PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail: "
			<< PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail << endl;	
	}

	void outputAllStats(ofstream& output, bool Do_phase1_only_or_not_bool, unsigned int read_num)
	{
		double tmpPerc;
		output << endl << "total reads #: " << read_num << endl << endl;
		
		tmpPerc = ((double)PE_pair_num_paired_complete_unique_phase1/read_num)*100;
		output << endl << "phase1 statics: " << endl << endl 
			<< "PE_pair_num_paired_complete_unique_phase1: " //<< endl
			<< PE_pair_num_paired_complete_unique_phase1 << " -- " << tmpPerc << "%" << endl;

		tmpPerc = ((double)PE_pair_num_paired_complete_multi_phase1/read_num)*100;
		output << "PE_pair_num_paired_complete_multi_phase1: " //<< endl
			<< PE_pair_num_paired_complete_multi_phase1 << " -- " << tmpPerc << "%" << endl;
		
		tmpPerc = ((double)PE_pair_num_paired_incomplete_unique_phase1/read_num)*100;
		output << "PE_pair_num_paired_incomplete_unique_phase1: " //<< endl
			<< PE_pair_num_paired_incomplete_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		
		tmpPerc = ((double)PE_pair_num_paired_incomplete_multi_phase1/read_num)*100;
		output << "PE_pair_num_paired_incomplete_multi_phase1: " //<< endl
			<< PE_pair_num_paired_incomplete_multi_phase1 << " -- " << tmpPerc << "%" << endl;	

		tmpPerc = ((double)PE_pair_num_unpaired_phase1/read_num)*100;
		output << "PE_pair_num_unpaired_phase1: " //<< endl
			<< PE_pair_num_unpaired_phase1 << " -- " << tmpPerc << "%" << endl;
		
		tmpPerc = ((double)PE_pair_num_unmapped_mappedToRepeatRegion_phase1/read_num)*100;
		output << "PE_pair_num_unmapped_mappedToRepeatRegion_phase1: " //<< endl
			<< PE_pair_num_unmapped_mappedToRepeatRegion_phase1 << " -- " << tmpPerc << "%" << endl;
		
		tmpPerc = ((double)PE_pair_num_unmapped_phase1/read_num)*100;
		output << "PE_pair_num_unmapped_phase1: " //<< endl
			<< PE_pair_num_unmapped_phase1 << " -- " << tmpPerc << "%" <<  endl << endl;

		tmpPerc = ((double)PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1/read_num)*100;
		output << "PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		
		tmpPerc = ((double)PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1/read_num)*100;
		output << "PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1 << " -- " << tmpPerc << "%" << endl << endl;	


		output << endl << endl << "fix unpaired reads output: " << endl << endl;
		
		tmpPerc = ((double)PE_pair_num_paired_complete_unique_fixUnpaired/read_num)*100;
		output << "PE_pair_num_paired_complete_unique_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_complete_unique_fixUnpaired << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_complete_multi_fixUnpaired/read_num)*100;
		output << "PE_pair_num_paired_complete_multi_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_complete_multi_fixUnpaired << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_incomplete_unique_fixUnpaired/read_num)*100;
		output << "PE_pair_num_paired_incomplete_unique_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_incomplete_unique_fixUnpaired << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_incomplete_multi_fixUnpaired/read_num)*100;
		output << "PE_pair_num_paired_incomplete_multi_fixUnpaired: " //<< endl
			<< PE_pair_num_paired_incomplete_multi_fixUnpaired << " -- " << tmpPerc << "%" << endl << endl;			

		tmpPerc = ((double)PE_pair_num_unpaired_complete_fixUnpaired/read_num)*100;
		output << "PE_pair_num_unpaired_complete_fixUnpaired: " //<< endl
			<< PE_pair_num_unpaired_complete_fixUnpaired << " -- " << tmpPerc << "%" << endl; 
		tmpPerc = ((double)PE_pair_num_unpaired_incomplete_fixUnpaired/read_num)*100;
		output << "PE_pair_num_unpaired_incomplete_fixUnpaired: " //<< endl
			<< PE_pair_num_unpaired_incomplete_fixUnpaired << " -- " << tmpPerc << "%" << endl << endl; 
		tmpPerc = ((double)PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair/read_num)*100;	
		output << "PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair: " 
			<< PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair/read_num)*100;
		output << "PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair: " 
			<< PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair << " -- " << tmpPerc << "%" << endl;


		output << endl << endl << "fix headTail output: " << endl << endl;
		tmpPerc = ((double)PE_pair_num_paired_complete_unique_fixHeadTail/read_num)*100;
		output << "PE_pair_num_paired_complete_unique_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_complete_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_complete_multi_fixHeadTail/read_num)*100;
		output << "PE_pair_num_paired_complete_multi_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_complete_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_incomplete_unique_fixHeadTail/read_num)*100;
		output << "PE_pair_num_paired_incomplete_unique_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_incomplete_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_incomplete_multi_fixHeadTail/read_num)*100;
		output << "PE_pair_num_paired_incomplete_multi_fixHeadTail: " //<< endl
			<< PE_pair_num_paired_incomplete_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;	

		tmpPerc = ((double)PE_pair_num_unpaired_complete_fixHeadTail/read_num)*100;	
		output << "PE_pair_num_unpaired_complete_fixHeadTail: " //<< endl
			<< PE_pair_num_unpaired_complete_fixHeadTail << " -- " << tmpPerc << "%" << endl; 
		tmpPerc = ((double)PE_pair_num_unpaired_incomplete_fixHeadTail/read_num)*100;
		output << "PE_pair_num_unpaired_incomplete_fixHeadTail: " //<< endl
			<< PE_pair_num_unpaired_incomplete_fixHeadTail << " -- " << tmpPerc << "%" << endl << endl; 

		tmpPerc = ((double)PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail/read_num)*100;	
		output << "PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail/read_num)*100;
		output << "PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail: "
			<< PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;			
		tmpPerc = ((double)PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail/read_num)*100;
		output << "PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail: "
			<< PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail/read_num)*100;
		output << "PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail: "
			<< PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;	
	}

	void outputFinalStats(ofstream& output, bool Do_phase1_only_or_not_bool//, unsigned int read_num
		)
	{
		//unsigned int complete_num = 
		//unsigned int incomplete_num = 
		unsigned int read_num;
	
		unsigned int bothEndsUnmapped_repeatRegion_num = PE_pair_num_unmapped_mappedToRepeatRegion_phase1;
		double bothEndsUnmapped_repeatRegion_perc = ((double)bothEndsUnmapped_repeatRegion_num/read_num)*100;
		unsigned int bothEndsUnmapped_num = PE_pair_num_unmapped_phase1;
		double bothEndsUnmapped_perc = ((double)bothEndsUnmapped_num/read_num)*100;
		unsigned int unmapped_dueToLowScore_num 
			= PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1
				+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair
				+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail
				+ PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail
				+ PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail;
		double unmapped_dueToLowScore_perc = ((double)unmapped_dueToLowScore_num/read_num)*100; 

		if(Do_phase1_only_or_not_bool)
		{
			read_num = PE_pair_num_paired_complete_unique_phase1 + PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_multi_phase1 + PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_unpaired_phase1 + PE_pair_num_unmapped_phase1
				+ bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		}
		else
		{
			read_num = 	PE_pair_num_paired_complete_unique_phase1 //+ PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_unique_fixUnpaired //+ PE_pair_num_paired_incomplete_unique_fixUnpaired
				+ PE_pair_num_paired_complete_unique_fixHeadTail + PE_pair_num_paired_incomplete_unique_fixHeadTail;
				+ PE_pair_num_paired_complete_multi_phase1 //+ PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_paired_complete_multi_fixUnpaired //+ PE_pair_num_paired_incomplete_multi_fixUnpaired
				+ PE_pair_num_paired_complete_multi_fixHeadTail + PE_pair_num_paired_incomplete_multi_fixHeadTail;
				+ PE_pair_num_unpaired_complete_fixUnpaired + PE_pair_num_unpaired_incomplete_fixUnpaired
				+ PE_pair_num_unpaired_complete_fixHeadTail + PE_pair_num_unpaired_incomplete_fixHeadTail
				+ bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		}

		unsigned int pair_unique_num;// = 
		double pair_unique_perc;
		unsigned int pair_multi_num;// = 
		double pair_multi_perc;
		unsigned int oneEndUnmapped_num;// = 
		double oneEndUnmapped_perc;
		// unsigned int bothEndsUnmapped_repeatRegion_num = PE_pair_num_unmapped_mappedToRepeatRegion_phase1;
		// double bothEndsUnmapped_repeatRegion_perc = ((double)bothEndsUnmapped_repeatRegion_num/read_num)*100;
		// unsigned int bothEndsUnmapped_num = PE_pair_num_unmapped_phase1;
		// double bothEndsUnmapped_perc = ((double)bothEndsUnmapped_num/read_num)*100;
		// unsigned int unmapped_dueToLowScore_num 
		// 	= PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1
		// 		+ PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1
		// 		+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair
		// 		+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair
		// 		+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail
		// 		+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail
		// 		+ PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail
		// 		+ PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail;
		// double unmapped_dueToLowScore_perc = ((double)unmapped_dueToLowScore_num/read_num)*100; 

		unsigned int unmapped_num = bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		double unmapped_perc = ((double)unmapped_num/read_num)*100;
		
		if(Do_phase1_only_or_not_bool)
		{
			pair_unique_num = PE_pair_num_paired_complete_unique_phase1 + PE_pair_num_paired_incomplete_unique_phase1;
			pair_multi_num = PE_pair_num_paired_complete_multi_phase1 + PE_pair_num_paired_incomplete_multi_phase1;
			oneEndUnmapped_num = PE_pair_num_unpaired_phase1;
		}
		else
		{
			pair_unique_num = PE_pair_num_paired_complete_unique_phase1 //+ PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_unique_fixUnpaired //+ PE_pair_num_paired_incomplete_unique_fixUnpaired
				+ PE_pair_num_paired_complete_unique_fixHeadTail + PE_pair_num_paired_incomplete_unique_fixHeadTail;
			pair_multi_num = PE_pair_num_paired_complete_multi_phase1 //+ PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_paired_complete_multi_fixUnpaired //+ PE_pair_num_paired_incomplete_multi_fixUnpaired
				+ PE_pair_num_paired_complete_multi_fixHeadTail + PE_pair_num_paired_incomplete_multi_fixHeadTail;
			oneEndUnmapped_num = PE_pair_num_unpaired_complete_fixUnpaired + PE_pair_num_unpaired_incomplete_fixUnpaired
				+ PE_pair_num_unpaired_complete_fixHeadTail + PE_pair_num_unpaired_incomplete_fixHeadTail;
		}

		pair_unique_perc = ((double)pair_unique_num/read_num)*100;
		pair_multi_perc = ((double)pair_multi_num/read_num)*100;
		double pair_perc = pair_unique_perc + pair_multi_perc;
		oneEndUnmapped_perc = ((double)oneEndUnmapped_num/read_num)*100;

		output << endl << "Total # of read pairs: " << read_num << endl;
		output << endl << endl << "unique-paired alignment read pairs: " << pair_unique_num << " -- " << pair_unique_perc << "%" << endl << endl;
		output << "multi-paired alignment read pairs: " << pair_multi_num << " -- " << pair_multi_perc << "%" << endl << endl;
		output << "paired alignment read pairs: " << pair_unique_num + pair_multi_num << " -- " << pair_perc << "%" << endl << endl;
		output << "unpaired read pairs: " << oneEndUnmapped_num << " --  " << oneEndUnmapped_perc << "%" << endl << endl;
		output << "unmapped read pairs: " << unmapped_num << " -- " << unmapped_perc << "%" << endl << endl;
		output << "           unmapped reads pairs (reason -- low mapping score): " 
			<< unmapped_dueToLowScore_num << " -- " << unmapped_dueToLowScore_perc << "%" << endl << endl;		
		output << "           unmapped reads pairs (reason -- perfectly mapped to repeat regions): " 
			<< bothEndsUnmapped_repeatRegion_num << " -- " << bothEndsUnmapped_repeatRegion_perc << "%" << endl << endl;		
		output << "           unmapped reads pairs (other reasons): " 
			<< bothEndsUnmapped_num << " -- " << bothEndsUnmapped_perc << "%" << endl << endl;
	}

	void outputFinalStats(ofstream& output, bool Do_phase1_only_or_not_bool, unsigned int read_num)
	{
		//unsigned int complete_num = 
		//unsigned int incomplete_num = 
		unsigned int pair_unique_num;// = 
		double pair_unique_perc;
		unsigned int pair_multi_num;// = 
		double pair_multi_perc;
		unsigned int oneEndUnmapped_num;// = 
		double oneEndUnmapped_perc;
		unsigned int bothEndsUnmapped_repeatRegion_num = PE_pair_num_unmapped_mappedToRepeatRegion_phase1;
		double bothEndsUnmapped_repeatRegion_perc = ((double)bothEndsUnmapped_repeatRegion_num/read_num)*100;
		unsigned int bothEndsUnmapped_num = PE_pair_num_unmapped_phase1;
		double bothEndsUnmapped_perc = ((double)bothEndsUnmapped_num/read_num)*100;
		unsigned int unmapped_dueToLowScore_num 
			= PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1
				+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair
				+ PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail
				+ PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail
				+ PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail
				+ PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail;
		double unmapped_dueToLowScore_perc = ((double)unmapped_dueToLowScore_num/read_num)*100; 

		unsigned int unmapped_num = bothEndsUnmapped_repeatRegion_num + bothEndsUnmapped_num + unmapped_dueToLowScore_num;
		double unmapped_perc = ((double)unmapped_num/read_num)*100;
		
		if(Do_phase1_only_or_not_bool)
		{
			pair_unique_num = PE_pair_num_paired_complete_unique_phase1 + PE_pair_num_paired_incomplete_unique_phase1;
			pair_multi_num = PE_pair_num_paired_complete_multi_phase1 + PE_pair_num_paired_incomplete_multi_phase1;
			oneEndUnmapped_num = PE_pair_num_unpaired_phase1;
		}
		else
		{
			pair_unique_num = PE_pair_num_paired_complete_unique_phase1 //+ PE_pair_num_paired_incomplete_unique_phase1
				+ PE_pair_num_paired_complete_unique_fixUnpaired //+ PE_pair_num_paired_incomplete_unique_fixUnpaired
				+ PE_pair_num_paired_complete_unique_fixHeadTail + PE_pair_num_paired_incomplete_unique_fixHeadTail;
			pair_multi_num = PE_pair_num_paired_complete_multi_phase1 //+ PE_pair_num_paired_incomplete_multi_phase1
				+ PE_pair_num_paired_complete_multi_fixUnpaired //+ PE_pair_num_paired_incomplete_multi_fixUnpaired
				+ PE_pair_num_paired_complete_multi_fixHeadTail + PE_pair_num_paired_incomplete_multi_fixHeadTail;
			oneEndUnmapped_num = PE_pair_num_unpaired_complete_fixUnpaired + PE_pair_num_unpaired_incomplete_fixUnpaired
				+ PE_pair_num_unpaired_complete_fixHeadTail + PE_pair_num_unpaired_incomplete_fixHeadTail;
		}

		pair_unique_perc = ((double)pair_unique_num/read_num)*100;
		pair_multi_perc = ((double)pair_multi_num/read_num)*100;
		double pair_perc = pair_unique_perc + pair_multi_perc;
		oneEndUnmapped_perc = ((double)oneEndUnmapped_num/read_num)*100;

		output << endl << endl << "unique-paired alignment reads: " << pair_unique_num << " -- " << pair_unique_perc << "%" << endl << endl;
		output << "multi-paired alignment reads: " << pair_multi_num << " -- " << pair_multi_perc << "%" << endl << endl;
		output << "paired alignment read: " << pair_unique_num + pair_multi_num << " -- " << pair_perc << "%" << endl << endl;
		output << "unpaired reads: " << oneEndUnmapped_num << " --  " << oneEndUnmapped_perc << "%" << endl << endl;
		output << "unmapped reads: " << unmapped_num << " -- " << unmapped_perc << "%" << endl << endl;
		output << "           unmapped reads (reason -- low mapping score): " 
			<< unmapped_dueToLowScore_num << " -- " << unmapped_dueToLowScore_perc << "%" << endl << endl;		
		output << "           unmapped reads (reason -- perfectly mapped to repeat regions): " 
			<< bothEndsUnmapped_repeatRegion_num << " -- " << bothEndsUnmapped_repeatRegion_perc << "%" << endl << endl;		
		output << "           unmapped reads (other reasons): " 
			<< bothEndsUnmapped_num << " -- " << bothEndsUnmapped_perc << "%" << endl << endl;
	}

	void getPhase1Stats()
	{
		int size = PE_pair_num_paired_complete_unique_phase1_thread.size();
		for(int tmp = 0; tmp < size; tmp++)
		{
			PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1 += PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1_thread[tmp];
			PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1 += PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1_thread[tmp];

			PE_pair_num_paired_complete_unique_phase1 += PE_pair_num_paired_complete_unique_phase1_thread[tmp];
			PE_pair_num_paired_complete_multi_phase1 += PE_pair_num_paired_complete_multi_phase1_thread[tmp];	
			PE_pair_num_paired_incomplete_unique_phase1 += PE_pair_num_paired_incomplete_unique_phase1_thread[tmp];
			PE_pair_num_paired_incomplete_multi_phase1 += PE_pair_num_paired_incomplete_multi_phase1_thread[tmp];

			PE_pair_num_unmapped_phase1 += PE_pair_num_unmapped_phase1_thread[tmp];	
			PE_pair_num_unmapped_mappedToRepeatRegion_phase1 += PE_pair_num_unmapped_mappedToRepeatRegion_phase1_thread[tmp];
			PE_pair_num_unpaired_phase1 += PE_pair_num_unpaired_phase1_thread[tmp]; 
		}

		// after adding the feature of filtering out low score alignments
		PE_pair_num_paired_complete_unique_phase1 
			= PE_pair_num_paired_complete_unique_phase1 - PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1;
		PE_pair_num_paired_complete_multi_phase1 
			= PE_pair_num_paired_complete_multi_phase1 - PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1;
	}

	void getFixUnpairedStats()
	{
		int size = PE_pair_num_paired_complete_unique_fixUnpaired_thread.size();
		for(int tmp = 0; tmp < size; tmp++)
		{
			PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair += PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair_thread[tmp];
			PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair += PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair_thread[tmp];

			PE_pair_num_paired_complete_unique_fixUnpaired += PE_pair_num_paired_complete_unique_fixUnpaired_thread[tmp];
			PE_pair_num_paired_complete_multi_fixUnpaired += PE_pair_num_paired_complete_multi_fixUnpaired_thread[tmp];
			PE_pair_num_paired_incomplete_unique_fixUnpaired += PE_pair_num_paired_incomplete_unique_fixUnpaired_thread[tmp];
			PE_pair_num_paired_incomplete_multi_fixUnpaired += PE_pair_num_paired_incomplete_multi_fixUnpaired_thread[tmp];
			
			PE_pair_num_unpaired_complete_fixUnpaired += PE_pair_num_unpaired_complete_fixUnpaired_thread[tmp];
			PE_pair_num_unpaired_incomplete_fixUnpaired += PE_pair_num_unpaired_incomplete_fixUnpaired_thread[tmp];
		}

		// after adding the feature of filtering out low score alignments
		PE_pair_num_paired_complete_unique_fixUnpaired 
			= PE_pair_num_paired_complete_unique_fixUnpaired - PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair;
		PE_pair_num_paired_complete_multi_fixUnpaired
			= PE_pair_num_paired_complete_multi_fixUnpaired - PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair; 
	}

	void getFixHeadTailStats()
	{
		int size = PE_pair_num_paired_complete_unique_fixHeadTail_thread.size();
		for(int tmp = 0; tmp < size; tmp++)
		{
			PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail += PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail_thread[tmp];
			PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail += PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail_thread[tmp];
			PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail += PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail_thread[tmp];				
			PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail += PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail_thread[tmp];			

			PE_pair_num_paired_complete_unique_fixHeadTail += PE_pair_num_paired_complete_unique_fixHeadTail_thread[tmp];
			PE_pair_num_paired_complete_multi_fixHeadTail += PE_pair_num_paired_complete_multi_fixHeadTail_thread[tmp];
			PE_pair_num_paired_incomplete_unique_fixHeadTail += PE_pair_num_paired_incomplete_unique_fixHeadTail_thread[tmp];
			PE_pair_num_paired_incomplete_multi_fixHeadTail += PE_pair_num_paired_incomplete_multi_fixHeadTail_thread[tmp];

			PE_pair_num_unpaired_complete_fixHeadTail += PE_pair_num_unpaired_complete_fixHeadTail_thread[tmp];
			PE_pair_num_unpaired_incomplete_fixHeadTail += PE_pair_num_unpaired_incomplete_fixHeadTail_thread[tmp];
		}	

		// after adding the feature of filtering out low score alignments
		PE_pair_num_paired_complete_unique_fixHeadTail
			= PE_pair_num_paired_complete_unique_fixHeadTail - PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail;
		PE_pair_num_paired_complete_multi_fixHeadTail
			= PE_pair_num_paired_complete_multi_fixHeadTail - PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail;
		PE_pair_num_paired_incomplete_unique_fixHeadTail
			= PE_pair_num_paired_incomplete_unique_fixHeadTail - PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail;
		PE_pair_num_paired_incomplete_multi_fixHeadTail
			= PE_pair_num_paired_incomplete_multi_fixHeadTail - PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail;

	}

	void getPhase1Stats_SE()
	{
		int size = SE_num_complete_unmapped_lowScore_unique_phase1_thread.size();
		for(int tmp = 0; tmp < size; tmp++)
		{
			SE_num_complete_unmapped_lowScore_unique_phase1 += SE_num_complete_unmapped_lowScore_unique_phase1_thread[tmp];
			SE_num_complete_unmapped_lowScore_multi_phase1 += SE_num_complete_unmapped_lowScore_multi_phase1_thread[tmp];	
			SE_num_complete_unique_phase1 += SE_num_complete_unique_phase1_thread[tmp];
			SE_num_complete_multi_phase1 += SE_num_complete_multi_phase1_thread[tmp];
			SE_num_incomplete_unique_phase1 += SE_num_incomplete_unique_phase1_thread[tmp];
			SE_num_incomplete_multi_phase1 += SE_num_incomplete_multi_phase1_thread[tmp];
			SE_num_unmapped_phase1 += SE_num_unmapped_phase1_thread[tmp];
			SE_num_unmapped_mappedToRepeatRegion_phase1 += SE_num_unmapped_mappedToRepeatRegion_phase1_thread[tmp];
		}
	}

	void getFixHeadTailStats_SE()
	{
		int size = SE_num_complete_unmapped_lowScore_unique_fixHeadTail_thread.size();
		for(int tmp = 0; tmp < size; tmp++)
		{
			SE_num_complete_unmapped_lowScore_unique_fixHeadTail += SE_num_complete_unmapped_lowScore_unique_fixHeadTail_thread[tmp];
			SE_num_complete_unmapped_lowScore_multi_fixHeadTail += SE_num_complete_unmapped_lowScore_multi_fixHeadTail_thread[tmp];
			SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail += SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail_thread[tmp];
			SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail += SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail_thread[tmp];

			SE_num_complete_unique_fixHeadTail += SE_num_complete_unique_fixHeadTail_thread[tmp];
			SE_num_complete_multi_fixHeadTail += SE_num_complete_multi_fixHeadTail_thread[tmp];
			SE_num_incomplete_unique_fixHeadTail += SE_num_incomplete_unique_fixHeadTail_thread[tmp];
			SE_num_incomplete_multi_fixHeadTail += SE_num_incomplete_multi_fixHeadTail_thread[tmp];
		}
	}

	void outputAllStats_SE_phase1(ofstream& output, unsigned int read_num)
	{
		double tmpPerc;
		output << endl << "total reads #: " << read_num << endl << endl;
		output << endl << "phase1 statics: " << endl << endl;
		tmpPerc = ((double)SE_num_complete_unique_phase1/read_num)*100;
		output << "SE_complete_unique_phase1: " << SE_num_complete_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_complete_multi_phase1/read_num)*100;
		output << "SE_complete_multi_phase1: " << SE_num_complete_multi_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unique_phase1/read_num)*100;
		output << "SE_incomplete_unique_phase1: " << SE_num_incomplete_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_multi_phase1/read_num)*100;
		output << "SE_incomplete_multi_phase1: " << SE_num_incomplete_multi_phase1 << " -- " << tmpPerc << "%" << endl;

		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_unique_phase1/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_unique_phase1: " << SE_num_complete_unmapped_lowScore_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_multi_phase1/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_multi_phase1: " << SE_num_complete_unmapped_lowScore_multi_phase1 << " -- " << tmpPerc << "%" << endl;

		tmpPerc = ((double)SE_num_unmapped_phase1/read_num)*100;
		output << "SE_num_unmapped_phase1: " << SE_num_unmapped_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_unmapped_mappedToRepeatRegion_phase1/read_num)*100;
		output << "SE_num_unmapped_mappedToRepeatRegion_phase1: " << SE_num_unmapped_mappedToRepeatRegion_phase1 << " -- " << tmpPerc << "%" << endl;		
	}

	void outputAllStats_SE_fixHeadTail(ofstream& output, unsigned int read_num)
	{
		double tmpPerc;
		output << endl << endl << "fix headTail output: " << endl << endl;
		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_unique_fixHeadTail/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_unique_fixHeadTail: " << SE_num_complete_unmapped_lowScore_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;				
		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_multi_fixHeadTail/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_multi_fixHeadTail: " << SE_num_complete_unmapped_lowScore_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail: " << SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail: " << SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;

		tmpPerc = ((double)SE_num_complete_unique_fixHeadTail/read_num)*100;
		output << "SE_num_complete_unique_fixHeadTail: " << SE_num_complete_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_complete_multi_fixHeadTail/read_num)*100;
		output << "SE_num_complete_multi_fixHeadTail: " << SE_num_complete_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unique_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_unique_fixHeadTail: " << SE_num_incomplete_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;										
		tmpPerc = ((double)SE_num_incomplete_multi_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_multi_fixHeadTail: " << SE_num_incomplete_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;		
	}

	void outputAllStats_SE(ofstream& output, unsigned int read_num)
	{
		double tmpPerc;
		output << endl << "total reads #: " << read_num << endl << endl;
		output << endl << "phase1 statics: " << endl << endl;
		tmpPerc = ((double)SE_num_complete_unique_phase1/read_num)*100;
		output << "SE_complete_unique_phase1: " << SE_num_complete_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_complete_multi_phase1/read_num)*100;
		output << "SE_complete_multi_phase1: " << SE_num_complete_multi_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unique_phase1/read_num)*100;
		output << "SE_incomplete_unique_phase1: " << SE_num_incomplete_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_multi_phase1/read_num)*100;
		output << "SE_incomplete_multi_phase1: " << SE_num_incomplete_multi_phase1 << " -- " << tmpPerc << "%" << endl;

		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_unique_phase1/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_unique_phase1: " << SE_num_complete_unmapped_lowScore_unique_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_multi_phase1/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_multi_phase1: " << SE_num_complete_unmapped_lowScore_multi_phase1 << " -- " << tmpPerc << "%" << endl;

		tmpPerc = ((double)SE_num_unmapped_phase1/read_num)*100;
		output << "SE_num_unmapped_phase1: " << SE_num_unmapped_phase1 << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_unmapped_mappedToRepeatRegion_phase1/read_num)*100;
		output << "SE_num_unmapped_mappedToRepeatRegion_phase1: " << SE_num_unmapped_mappedToRepeatRegion_phase1 << " -- " << tmpPerc << "%" << endl;

		output << endl << endl << "fix headTail output: " << endl << endl;
		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_unique_fixHeadTail/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_unique_fixHeadTail: " << SE_num_complete_unmapped_lowScore_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;				
		tmpPerc = ((double)SE_num_complete_unmapped_lowScore_multi_fixHeadTail/read_num)*100;
		output << "SE_num_complete_unmapped_lowScore_multi_fixHeadTail: " << SE_num_complete_unmapped_lowScore_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail: " << SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail: " << SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;

		tmpPerc = ((double)SE_num_complete_unique_fixHeadTail/read_num)*100;
		output << "SE_num_complete_unique_fixHeadTail: " << SE_num_complete_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_complete_multi_fixHeadTail/read_num)*100;
		output << "SE_num_complete_multi_fixHeadTail: " << SE_num_complete_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;
		tmpPerc = ((double)SE_num_incomplete_unique_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_unique_fixHeadTail: " << SE_num_incomplete_unique_fixHeadTail << " -- " << tmpPerc << "%" << endl;										
		tmpPerc = ((double)SE_num_incomplete_multi_fixHeadTail/read_num)*100;
		output << "SE_num_incomplete_multi_fixHeadTail: " << SE_num_incomplete_multi_fixHeadTail << " -- " << tmpPerc << "%" << endl;
	}

	void outputFinalStats_SE(ofstream& output, unsigned int read_num)
	{
		unsigned int unique_num = SE_num_complete_unique_phase1 + SE_num_complete_unique_fixHeadTail + SE_num_incomplete_unique_fixHeadTail;
		double unique_perc = ((double)unique_num/read_num)*100;
		unsigned int multi_num = SE_num_complete_multi_phase1 + SE_num_complete_multi_fixHeadTail + SE_num_incomplete_multi_fixHeadTail;
		double multi_perc = ((double)multi_num/read_num)*100;
		unsigned int mapped_num = unique_num + multi_num;
		double mapped_perc = ((double)mapped_num/read_num)*100;
		unsigned int unmap_num = SE_num_complete_unmapped_lowScore_unique_phase1 + SE_num_complete_unmapped_lowScore_multi_phase1
			+ SE_num_unmapped_phase1 + SE_num_unmapped_mappedToRepeatRegion_phase1
			+ SE_num_complete_unmapped_lowScore_unique_fixHeadTail + SE_num_complete_unmapped_lowScore_multi_fixHeadTail
			+ SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail + SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail;
		double unmap_perc = ((double)unmap_num/read_num)*100;
		unsigned int unmap_lowScore_num = SE_num_complete_unmapped_lowScore_unique_phase1 + SE_num_complete_unmapped_lowScore_multi_phase1
			+ SE_num_complete_unmapped_lowScore_unique_fixHeadTail + SE_num_complete_unmapped_lowScore_multi_fixHeadTail
			+ SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail + SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail;
		double unmap_lowScore_perc = ((double)unmap_lowScore_num/read_num)*100;
		double unmap_repeatRegion_perc = ((double)SE_num_unmapped_mappedToRepeatRegion_phase1/read_num)*100;
		double unmap_otherReason_perc = ((double)SE_num_unmapped_phase1/read_num)*100;
		
		output << endl << endl << "unique alignment reads: " << unique_num << " -- " << unique_perc << "%" << endl << endl;
		output << "multi alignment reads: " << multi_num << " -- " << multi_perc << "%" << endl << endl;
		output << "mapped reads: " << mapped_num << " -- " << mapped_perc << "%" << endl << endl;
		output << "unmapped reads: " << unmap_num << " -- " << unmap_perc << "%" << endl << endl;
		output << "           unmapped reads (reason -- low mapping score): " 
			<< unmap_lowScore_num << " -- " << unmap_lowScore_perc << "%" << endl << endl;		
		output << "           unmapped reads (reason -- perfectly mapped to repeat regions): " 
			<< SE_num_unmapped_mappedToRepeatRegion_phase1 << " -- " << unmap_repeatRegion_perc << "%" << endl << endl;		
		output << "           unmapped reads (other reasons): " 
			<< SE_num_unmapped_phase1 << " -- " << unmap_otherReason_perc << "%" << endl << endl;		
	}

	void initiate_stats_info_SE(int thread_num)
	{
		for(int tmp = 0; tmp < thread_num; tmp++)
		{
			SE_num_complete_unmapped_lowScore_unique_phase1_thread.push_back(0);
			SE_num_complete_unmapped_lowScore_multi_phase1_thread.push_back(0);

			SE_num_complete_unique_phase1_thread.push_back(0);
			SE_num_complete_multi_phase1_thread.push_back(0);
			SE_num_incomplete_unique_phase1_thread.push_back(0);
			SE_num_incomplete_multi_phase1_thread.push_back(0);

			SE_num_unmapped_phase1_thread.push_back(0);
			SE_num_unmapped_mappedToRepeatRegion_phase1_thread.push_back(0);			

			SE_num_complete_unmapped_lowScore_unique_fixHeadTail_thread.push_back(0);
			SE_num_complete_unmapped_lowScore_multi_fixHeadTail_thread.push_back(0); 
			SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail_thread.push_back(0);
			SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail_thread.push_back(0);

			SE_num_complete_unique_fixHeadTail_thread.push_back(0);
			SE_num_complete_multi_fixHeadTail_thread.push_back(0);
			SE_num_incomplete_unique_fixHeadTail_thread.push_back(0);
			SE_num_incomplete_multi_fixHeadTail_thread.push_back(0);
		}
	}

	void initiate_stats_info_PE(int thread_num)
	{
		for(int tmp = 0; tmp < thread_num; tmp++)
		{
			PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1_thread.push_back(0);
			PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1_thread.push_back(0);
			PE_pair_num_paired_complete_unique_phase1_thread.push_back(0);
			PE_pair_num_paired_complete_multi_phase1_thread.push_back(0);	
			PE_pair_num_paired_incomplete_unique_phase1_thread.push_back(0);
			PE_pair_num_paired_incomplete_multi_phase1_thread.push_back(0);

			PE_pair_num_unmapped_phase1_thread.push_back(0);
			PE_pair_num_unmapped_mappedToRepeatRegion_phase1_thread.push_back(0);
			PE_pair_num_unpaired_phase1_thread.push_back(0);

			PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair_thread.push_back(0);
			PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair_thread.push_back(0);
			PE_pair_num_paired_complete_unique_fixUnpaired_thread.push_back(0);
			PE_pair_num_paired_complete_multi_fixUnpaired_thread.push_back(0);
			PE_pair_num_paired_incomplete_unique_fixUnpaired_thread.push_back(0);
	 		PE_pair_num_paired_incomplete_multi_fixUnpaired_thread.push_back(0);

			PE_pair_num_unpaired_complete_fixUnpaired_thread.push_back(0);
			PE_pair_num_unpaired_incomplete_fixUnpaired_thread.push_back(0);

			PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail_thread.push_back(0);
			PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail_thread.push_back(0);
			PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail_thread.push_back(0);
			PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail_thread.push_back(0);
			PE_pair_num_paired_complete_unique_fixHeadTail_thread.push_back(0);
			PE_pair_num_paired_complete_multi_fixHeadTail_thread.push_back(0);
			PE_pair_num_paired_incomplete_unique_fixHeadTail_thread.push_back(0);
			PE_pair_num_paired_incomplete_multi_fixHeadTail_thread.push_back(0);

			PE_pair_num_unpaired_complete_fixHeadTail_thread.push_back(0);
			PE_pair_num_unpaired_incomplete_fixHeadTail_thread.push_back(0);
		}
	}

	void increLowScoreComplete_phase1(int thread, bool unique_bool)
	{
		if(unique_bool)
			PE_pair_num_paired_complete_unmapped_lowScore_unique_phase1_thread[thread] ++;
		else
			PE_pair_num_paired_complete_unmapped_lowScore_multi_phase1_thread[thread] ++;
	}

	void increLowScoreComplete_fixUnpair(int thread, bool unique_bool)
	{
		if(unique_bool)
			PE_pair_num_paired_complete_unmapped_lowScore_unique_fixUnpair_thread[thread] ++;
		else
			PE_pair_num_paired_complete_unmapped_lowScore_multi_fixUnpair_thread[thread] ++;
	}

	void increLowScoreComplete_fixHeadTail(int thread, bool complete_bool, bool unique_bool)
	{
		if(complete_bool)
		{	
			if(unique_bool)
				PE_pair_num_paired_complete_unmapped_lowScore_unique_fixHeadTail_thread[thread] ++;
			else
				PE_pair_num_paired_complete_unmapped_lowScore_multi_fixHeadTail_thread[thread] ++;				
		}
		else
		{
			if(unique_bool)
				PE_pair_num_paired_incomplete_unmapped_lowScore_unique_fixHeadTail_thread[thread] ++;
			else
				PE_pair_num_paired_incomplete_unmapped_lowScore_multi_fixHeadTail_thread[thread] ++;
		}
	}

	void increPairedNum_fixHeadTail( int thread, bool allFinalPairAlignmentCompleteBool, 
		bool unique_bool)
	{
		if(allFinalPairAlignmentCompleteBool)
		{
			if(unique_bool)
				PE_pair_num_paired_complete_unique_fixHeadTail_thread[thread] ++;
			else
				PE_pair_num_paired_complete_multi_fixHeadTail_thread[thread] ++;
		}
		else
		{
			if(unique_bool)
				PE_pair_num_paired_incomplete_unique_fixHeadTail_thread[thread] ++;
			else
				PE_pair_num_paired_incomplete_multi_fixHeadTail_thread[thread] ++;
		}
	}

	void increUnpairedNum_fixHeadTail(int thread, bool allUnpairedAlignmentCompleteBool)
	{
		if(allUnpairedAlignmentCompleteBool)
			PE_pair_num_unpaired_complete_fixHeadTail_thread[thread] ++;
		else
			PE_pair_num_unpaired_incomplete_fixHeadTail_thread[thread] ++;
	}
	
	void increNum_fixUnpaired(int thread, bool pairExistsBool, 
		bool allAlignmentCompleteBool, bool allUnpairedAlignmentCompleteBool, 
		bool unique_bool)
	{
		if(pairExistsBool && allAlignmentCompleteBool) // some pair exists, all completed, print out paired SAM info
		{
			if(unique_bool)
				PE_pair_num_paired_complete_unique_fixUnpaired_thread[thread] ++;
			else
				PE_pair_num_paired_complete_multi_fixUnpaired_thread[thread] ++;
		}
		else if(pairExistsBool && (!allAlignmentCompleteBool)) // pair exists, incomplete
		{
			if(unique_bool)
				PE_pair_num_paired_incomplete_unique_fixUnpaired_thread[thread] ++;
			else
				PE_pair_num_paired_incomplete_multi_fixUnpaired_thread[thread] ++;
		}
		else if((!pairExistsBool) && (allUnpairedAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
		{
			PE_pair_num_unpaired_complete_fixUnpaired_thread[thread] ++;
		}
		else // no pair exists, incomplete, print out alignInfo
		{
			PE_pair_num_unpaired_incomplete_fixUnpaired_thread[thread] ++;
		}
	}

	void increPairedNum_phase1(int thread, bool complete_bool, bool unique_bool)
	{
		//read_pair_num_paired_phase1_thread[thread] ++;
		if(complete_bool)
		{
				if(unique_bool)
					PE_pair_num_paired_complete_unique_phase1_thread[thread] ++;
				else
					PE_pair_num_paired_complete_multi_phase1_thread[thread]++;
		}
		else
		{
			if(unique_bool)
				PE_pair_num_paired_incomplete_unique_phase1_thread[thread] ++;
			else
				PE_pair_num_paired_incomplete_multi_phase1_thread[thread]++;
		}
	}

	void increUnpairedNum_phase1(int thread, bool exist_bool, bool repeat_bool)
	{
		if(exist_bool)
		{
			PE_pair_num_unpaired_phase1_thread[thread] ++;
		}
		else if(repeat_bool)
		{
			PE_pair_num_unmapped_mappedToRepeatRegion_phase1_thread[thread] ++;
		}
		else
		{
			PE_pair_num_unmapped_phase1_thread[thread] ++;
		}
	}

	void increCompleteOrTooLowScore_phase1_SE(int threadNO, bool completeAlignmentScore_tooLow, bool uniqueBool)
	{
		if(completeAlignmentScore_tooLow)
		{
			if(uniqueBool)
			{
				SE_num_complete_unmapped_lowScore_unique_phase1_thread[threadNO] ++;
			}
			else
			{
				SE_num_complete_unmapped_lowScore_multi_phase1_thread[threadNO] ++;
			}
		}
		else
		{
			if(uniqueBool)
			{
				SE_num_complete_unique_phase1_thread[threadNO] ++;
			}
			else
			{
				SE_num_complete_multi_phase1_thread[threadNO] ++;
			}
		}
	}

	void increIncomplete_phase1_SE(int threadNO, bool uniqueBool)
	{
		if(uniqueBool)
		{
			SE_num_incomplete_unique_phase1_thread[threadNO] ++;
		}
		else
		{
			SE_num_incomplete_multi_phase1_thread[threadNO] ++;
		}
	}

	void increUnmap_phase1_SE(int threadNO, bool mappedToRepeatRegionBool)
	{
		if(mappedToRepeatRegionBool)
		{
			SE_num_unmapped_mappedToRepeatRegion_phase1_thread[threadNO] ++;
		}
		else
		{
			SE_num_unmapped_phase1_thread[threadNO] ++;
		}
	}

	void increNum_fixHeadTail_SE(int threadNO, bool completeOrNotBool, bool uniqueOrNotBool, bool lowScoreOrNotBool)
	{
		if(completeOrNotBool)
		{
			if(uniqueOrNotBool)
			{
				if(lowScoreOrNotBool)
				{
					SE_num_complete_unmapped_lowScore_unique_fixHeadTail_thread[threadNO] ++;
				}
				else
				{
					SE_num_complete_unique_fixHeadTail_thread[threadNO] ++;
				}
			}
			else
			{
				if(lowScoreOrNotBool)
				{
					SE_num_complete_unmapped_lowScore_multi_fixHeadTail_thread[threadNO] ++;
				}
				else
				{
					SE_num_complete_multi_fixHeadTail_thread[threadNO] ++;
				}
			}
		}
		else
		{
			if(uniqueOrNotBool)
			{
				if(lowScoreOrNotBool)
				{
					SE_num_incomplete_unmapped_lowScore_unique_fixHeadTail_thread[threadNO] ++;
				}
				else
				{
					SE_num_incomplete_unique_fixHeadTail_thread[threadNO] ++;
				}
			}
			else
			{
				if(lowScoreOrNotBool)
				{
					SE_num_incomplete_unmapped_lowScore_multi_fixHeadTail_thread[threadNO] ++;
				}
				else
				{
					SE_num_incomplete_multi_fixHeadTail_thread[threadNO] ++;
				}
			}
		}
	}
};

#endif