// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef READ_INFO_H
#define READ_INFO_H

#include <string>
#include <string.h>

using namespace std;

/*
inline char getCharRevComp(char ch)
{
	int chInt = ch - 'A';
	static const char alphatChar[26] = {'T', 'N', 'G', 'N', 'N', 'N', 'C',
		'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'A',
		'N', 'N', 'N', 'N', 'N', 'N'};
	return alphatChar[chInt];
}

string getRcmSeq(const string& readSeq)
{
	int readSeqLength = readSeq.length();

	char readRcmSeqChar[readSeqLength];

	readRcmSeqChar[0] = getCharRevComp(readSeq.at(readSeqLength-1));

	for(int tmp = 1; tmp < readSeqLength; tmp ++)
	{
		readRcmSeqChar[tmp] = getCharRevComp((readSeq.at(readSeqLength - tmp - 1)));
	}
	string rcmSeq = readRcmSeqChar;
	return rcmSeq.substr(0, readSeqLength);
}

string getRevSeq(const string& qualitySeq)
{
	int qualitySeqLength = qualitySeq.length();

	char revQualitySeqChar[qualitySeqLength];

	revQualitySeqChar[0] = qualitySeq.at(qualitySeqLength-1);
	for(int tmp = 1; tmp < qualitySeqLength; tmp++)
	{
		revQualitySeqChar[tmp] = qualitySeq.at(qualitySeqLength-1-tmp);
	}
	string revSeq = revQualitySeqChar;
	return revSeq.substr(0, qualitySeqLength);
}*/

class Ori_Read_Info
{
//public:
private:
	string readName;
	int readSeqLength;
	string readSeq;
	string rcmReadSeq;
	string readQual;
	string rcmReadQual;
	//string readPeInfo;
public:
	int returnReadNameSize()
	{
		return readName.length();
	}
	void assignReadName(const string& tmpReadName)
	{
		readName = tmpReadName;
	}
	void assignReadSeqLength(int tmpReadSeqLength)
	{
		readSeqLength = tmpReadSeqLength;
	}
	void assignReadSeq(const string& tmpReadSeq)
	{
		readSeq = tmpReadSeq;
	}
	void assignRcmReadSeq(const string& tmpRcmReadSeq)
	{
		rcmReadSeq = tmpRcmReadSeq;
	}
	void assignReadQual(const string& tmpReadQual)
	{
		readQual = tmpReadQual;
	}
	void assignRcmReadQual(const string& tmpRcmReadQual)
	{
		rcmReadQual = tmpRcmReadQual;
	}
	string returnReadName()
	{
		return readName;
	}
	string returnReadSeq()
	{
		return readSeq;
	}
	string returnRcmReadSeq()
	{
		return rcmReadSeq;
	}
	string returnReadQual()
	{
		return readQual;
	}
	string returnRcmReadQual()
	{
		return rcmReadQual;
	}
	bool getSeqLength()
	{
		//return readSeq.length();
		readSeqLength = readSeq.length();
		if(readSeqLength <= 0)
		{
			return false;
		}	
		else
		{
			return true;
		}
	}
	int returnReadSeqLength()
	{
		return readSeqLength;
	}
	int returnReadLength()
	{
		return readSeqLength;
	}
};

class PE_Read_Info
{
//public:
private:
	Ori_Read_Info readInfo_pe1;
	Ori_Read_Info readInfo_pe2; 
public:
	bool checkPEreadsOverlapped()
	{
		// fix me
		//string read2prefix = readInfo_pe2.

		return false;
	}

	bool notContainMarkedMismatchQualSeq(const string& tmpQualSeq, const string& markedMismatchQualStr)
	{
		int foundMarkedBase = tmpQualSeq.find(markedMismatchQualStr);
		if(foundMarkedBase == tmpQualSeq.npos)
			return true;
		else
			return false;
	}
	bool checkSegTrustableOrNot(int tmpAlignInfoType, int segStartInRead, int segEndInRead)
	{
		string tmpQualSeq = this->returnSegQualSeq(tmpAlignInfoType, segStartInRead, segEndInRead);
		//cout << "tmpQualSeq: " << endl << tmpQualSeq << endl;
		// fix me, useing 'notContainMarkedMismatchQualSeq' for evaluation
		bool segSeqConfident = this->notContainMarkedMismatchQualSeq(tmpQualSeq, "a");
		return segSeqConfident;
	}
	bool checkConfidenceInShortAnchorHeadSeq(int tmpAlignInfoType, int headSeqLength)
	{
		string headQualSeq = this->returnHeadQualSeq(tmpAlignInfoType, headSeqLength);
		
		// fix me, useing 'notContainMarkedMismatchQualSeq' for evaluation
		bool headSeqConfident = this->notContainMarkedMismatchQualSeq(headQualSeq, "a");
		return headSeqConfident;
	}

	bool checkConfidenceInShortAnchorTailSeq(int tmpAlignInfoType, int tailSeqLength)
	{
		string tailQualSeq = this->returnTailQualSeq(tmpAlignInfoType, tailSeqLength);
	
		// fix me, useing 'notContainMarkedMismatchQualSeq' for evaluation 
		bool tailSeqConfident = this->notContainMarkedMismatchQualSeq(tailQualSeq, "a");
		return tailSeqConfident;
	}

	void initiateReadInfo(const string& readName_1, const string& readName_2,
		const string& readSeq_1, const string& readSeq_2,
		const string& readQualitySeq_1, const string& readQualitySeq_2, bool fasta_or_fastq, bool SE_or_PE_bool)
	{
		// cout << "readName_1: " << endl << readName_1 << endl;
		// cout << "readSeq_1: " << endl << readSeq_1 << endl;
		// cout << "rcmReadSeq_1: " << endl << getRcmSeq(readSeq_1) << endl;
		readInfo_pe1.assignReadName(readName_1);
		readInfo_pe1.assignReadSeq(readSeq_1);
		readInfo_pe1.assignRcmReadSeq(getRcmSeq(readSeq_1));
		readInfo_pe1.assignReadSeqLength((readInfo_pe1.returnReadSeq()).length());
		// cout << "readLength_1: " << endl << (readInfo_pe1.returnReadSeq()).length() << endl;
		// cout << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
		if(!SE_or_PE_bool)
		{	
			readInfo_pe2.assignReadName(readName_2);
			readInfo_pe2.assignReadSeq(readSeq_2);		
			readInfo_pe2.assignRcmReadSeq(getRcmSeq(readSeq_2));
			readInfo_pe2.assignReadSeqLength((readInfo_pe2.returnReadSeq()).length());	
		}
			
		if(fasta_or_fastq)
		{
			readInfo_pe1.assignReadQual("*");
			readInfo_pe1.assignRcmReadQual("*");
			if(!SE_or_PE_bool)
			{	
				readInfo_pe2.assignReadQual("*");
				readInfo_pe2.assignRcmReadQual("*");
			}
		}
		else
		{
			// cout << "readQual: " << endl << readQualitySeq_1 << endl;
			// cout << "rcmReadQual: " << endl << getRevSeq(readQualitySeq_1) << endl;
			readInfo_pe1.assignReadQual(readQualitySeq_1);
			readInfo_pe1.assignRcmReadQual(getRevSeq(readQualitySeq_1));
			if(!SE_or_PE_bool)
			{
				readInfo_pe2.assignReadQual(readQualitySeq_2);
				readInfo_pe2.assignRcmReadQual(getRevSeq(readQualitySeq_2));
			}
		}		
	}

	void initiateReadInfo(const string& readName_1, const string& readName_2,
		const string& readSeq_1, const string& readSeq_2,
		const string& readQualitySeq_1, const string& readQualitySeq_2, bool fasta_or_fastq)
	{
		initiateReadInfo_PE(readName_1, readName_2, readSeq_1, readSeq_2, 
			readQualitySeq_1, readQualitySeq_2, fasta_or_fastq);
	}	

	void initiateReadInfo_PE(const string& readName_1, const string& readName_2,
		const string& readSeq_1, const string& readSeq_2,
		const string& readQualitySeq_1, const string& readQualitySeq_2, bool fasta_or_fastq)
	{
		readInfo_pe1.assignReadName(readName_1);
		readInfo_pe2.assignReadName(readName_2);
		readInfo_pe1.assignReadSeq(readSeq_1);
		readInfo_pe2.assignReadSeq(readSeq_2);
		readInfo_pe1.assignRcmReadSeq(getRcmSeq(readSeq_1));
		readInfo_pe2.assignRcmReadSeq(getRcmSeq(readSeq_2));

		if(fasta_or_fastq)
		{
			readInfo_pe1.assignReadQual("*");
			readInfo_pe2.assignReadQual("*");
			readInfo_pe1.assignRcmReadQual("*");
			readInfo_pe2.assignRcmReadQual("*");
		}
		else
		{
			readInfo_pe1.assignReadQual(readQualitySeq_1);
			readInfo_pe2.assignReadQual(readQualitySeq_2);
			readInfo_pe1.assignRcmReadQual(getRevSeq(readQualitySeq_1));
			readInfo_pe2.assignRcmReadQual(getRevSeq(readQualitySeq_2));
		}
		readInfo_pe1.assignReadSeqLength((readInfo_pe1.returnReadSeq()).length());
		readInfo_pe2.assignReadSeqLength((readInfo_pe2.returnReadSeq()).length());			
	}

	void getFastaFormatReadInfo_new(const string& readName_1, const string& readName_2, 
		const string& readSeq_1, const string& readSeq_2)
	{
		readInfo_pe1.assignReadName(readName_1);
		readInfo_pe2.assignReadName(readName_2);
		readInfo_pe1.assignReadSeq(readSeq_1);
		readInfo_pe2.assignReadSeq(readSeq_2);
		readInfo_pe1.assignRcmReadSeq(getRcmSeq(readSeq_1));
		readInfo_pe2.assignRcmReadSeq(getRcmSeq(readSeq_2));

		readInfo_pe1.assignReadSeqLength((readInfo_pe1.returnReadSeq()).length());
		readInfo_pe2.assignReadSeqLength((readInfo_pe2.returnReadSeq()).length());	
	}

	int returnAlignmentScoreMinOutput_withComplement_perHundredBases(int complementPerHundredBases)
	{
		int tmpPairReadLength = this->returnReadLength_end1() + this->returnReadLength_end2();
		int complement = (tmpPairReadLength/100) * complementPerHundredBases;
		return (tmpPairReadLength - complement);
	}

	int returnAlignmentScoreMinOutput_withComplement_perHundredBases_SE(int complementPerHundredBases)
	{
		int tmpReadLength = this->returnReadLength_SE();// + this->returnReadLength_end2();
		int complement = (tmpReadLength/100) * complementPerHundredBases;
		return (tmpReadLength - complement);
	}

	// int returnAlignmentScoreMinOutput_withComplement(int complement)
	// {
	// 	int tmpPairReadLength = this->returnReadLength_end1() + this->returnReadLength_end2();
	// 	return (tmpPairReadLength - complement);
	// }
	string returnReadQual_SE()
	{
		return readInfo_pe1.returnReadQual();
	}
	string returnQualitySeq_SE()
	{
		return readInfo_pe1.returnReadQual();
	}
	string returnRcmQualitySeq_SE()
	{
		return readInfo_pe1.returnRcmReadQual();
	}		
	string returnReadName_SE()
	{
		return readInfo_pe1.returnReadName();
	}
	string returnReadSeq_SE()
	{
		return readInfo_pe1.returnReadSeq();
	}
	string returnRcmReadSeq_SE()
	{
		return readInfo_pe1.returnRcmReadSeq();
	}
	int returnReadLength_SE()
	{
		return readInfo_pe1.returnReadLength();
	}
	int returnReadSeqLength_SE()
	{
		return readInfo_pe1.returnReadLength();
	}



	string returnReadQual_1()
	{
		return readInfo_pe1.returnReadQual();
	}
	string returnQualitySeq_1()
	{
		return readInfo_pe1.returnReadQual();
	}
	string returnRcmQualitySeq_1()
	{
		return readInfo_pe1.returnRcmReadQual();
	}		
	string returnReadQual_2()
	{
		return readInfo_pe2.returnReadQual();
	}
	string returnQualitySeq_2()
	{
		return readInfo_pe2.returnReadQual();
	}	
	string returnRcmQualitySeq_2()
	{
		return readInfo_pe2.returnRcmReadQual();
	}	
	string returnReadName_1()
	{
		return readInfo_pe1.returnReadName();
	}
	string returnReadName_2()
	{
		return readInfo_pe2.returnReadName();
	}	
	string returnReadSeq_1()
	{
		return readInfo_pe1.returnReadSeq();
	}
	string returnReadSeq_2()
	{
		return readInfo_pe2.returnReadSeq();
	}
	string returnRcmReadSeq_1()
	{
		return readInfo_pe1.returnRcmReadSeq();
	}
	string returnRcmReadSeq_2()
	{
		return readInfo_pe2.returnRcmReadSeq();
	}
	string returnReadName_beforeSlash_1()
	{
		int readNameSize_1 = ((readInfo_pe1).returnReadName()).length();
		return (readInfo_pe1.returnReadName()).substr(0, readNameSize_1 - 2);
	}
	string returnReadName_beforeSlash_2()
	{
		int readNameSize_2 = ((readInfo_pe2).returnReadName()).length();
		return (readInfo_pe2.returnReadName()).substr(0, readNameSize_2 - 2);
	}
	int returnReadNameSize_1()
	{
		return ((readInfo_pe1).returnReadName()).length();
	}
	int returnReadNameSize_2()
	{
		return ((readInfo_pe2).returnReadName()).length();
	}
	int returnReadLength_end1()
	{
		return readInfo_pe1.returnReadLength();
	}
	int returnReadLength_end2()
	{
		return readInfo_pe2.returnReadLength();
	}
	int returnReadSeqLength_1()
	{
		return readInfo_pe1.returnReadLength();
	}
	int returnReadSeqLength_2()
	{
		return readInfo_pe2.returnReadLength();
	}
	string returnSegQualSeq(int alignInfoType, int segStartInRead, int segEndInRead)
	{
		if(alignInfoType == 1)
		{
			return ((readInfo_pe1).returnReadQual()).substr(segStartInRead - 1, segEndInRead - segStartInRead + 1);
		}
		else if(alignInfoType == 2)
		{
			return ((readInfo_pe1).returnRcmReadQual()).substr(segStartInRead - 1, segEndInRead - segStartInRead + 1);
		}
		else if(alignInfoType == 3)
		{
			return ((readInfo_pe2).returnReadQual()).substr(segStartInRead - 1, segEndInRead - segStartInRead + 1);
		}
		else
		{
			return ((readInfo_pe2).returnRcmReadQual()).substr(segStartInRead - 1, segEndInRead - segStartInRead + 1);
		}		
	}
	string returnIncompleteLongHeadSeq(int alignInfoType, int incompleteHeadLength)
	{
		if(alignInfoType == 1)
		{
			return ((readInfo_pe1).returnReadSeq()).substr(0, incompleteHeadLength);
		}
		else if(alignInfoType == 2)
		{
			return ((readInfo_pe1).returnRcmReadSeq()).substr(0, incompleteHeadLength);
		}
		else if(alignInfoType == 3)
		{
			return ((readInfo_pe2).returnReadSeq()).substr(0, incompleteHeadLength);
		}
		else
		{
			return ((readInfo_pe2).returnRcmReadSeq()).substr(0, incompleteHeadLength);
		}		
	}
	string returnHeadQualSeq(int alignInfoType, int headLength)
	{
		if(alignInfoType == 1)
		{
			return ((readInfo_pe1).returnReadQual()).substr(0, headLength);
		}
		else if(alignInfoType == 2)
		{
			return ((readInfo_pe1).returnRcmReadQual()).substr(0, headLength);
		}
		else if(alignInfoType == 3)
		{
			return ((readInfo_pe2).returnReadQual()).substr(0, headLength);
		}
		else
		{
			return ((readInfo_pe2).returnRcmReadQual()).substr(0, headLength);
		}		
	}
	string returnIncompleteLongTailSeq(int alignInfoType, int incompleteTailLength)
	{
		if(alignInfoType == 1)
		{
			int readLength = readInfo_pe1.returnReadSeqLength();
			int substrStartLocInRead = readLength - incompleteTailLength + 1;
			return ((readInfo_pe1).returnReadSeq()).substr(substrStartLocInRead - 1, incompleteTailLength);
		}
		else if(alignInfoType == 2)
		{
			int readLength = readInfo_pe1.returnReadSeqLength();
			int substrStartLocInRead = readLength - incompleteTailLength + 1;
			return ((readInfo_pe1).returnRcmReadSeq()).substr(substrStartLocInRead - 1, incompleteTailLength);
		}
		else if(alignInfoType == 3)
		{
			int readLength = readInfo_pe2.returnReadSeqLength();
			int substrStartLocInRead = readLength - incompleteTailLength + 1;
			return ((readInfo_pe2).returnReadSeq()).substr(substrStartLocInRead - 1, incompleteTailLength);
		}
		else
		{
			int readLength = readInfo_pe2.returnReadSeqLength();
			int substrStartLocInRead = readLength - incompleteTailLength + 1;
			return ((readInfo_pe2).returnRcmReadSeq()).substr(substrStartLocInRead - 1, incompleteTailLength);
		}		
	}
	string returnTailQualSeq(int alignInfoType, int tailLength)
	{
		if(alignInfoType == 1)
		{
			int readLength = readInfo_pe1.returnReadSeqLength();
			int substrStartLocInRead = readLength - tailLength + 1;
			return ((readInfo_pe1).returnReadQual()).substr(substrStartLocInRead - 1, tailLength);
		}
		else if(alignInfoType == 2)
		{
			int readLength = readInfo_pe1.returnReadSeqLength();
			int substrStartLocInRead = readLength - tailLength + 1;
			return ((readInfo_pe1).returnRcmReadQual()).substr(substrStartLocInRead - 1, tailLength);
		}
		else if(alignInfoType == 3)
		{
			int readLength = readInfo_pe2.returnReadSeqLength();
			int substrStartLocInRead = readLength - tailLength + 1;
			return ((readInfo_pe2).returnReadQual()).substr(substrStartLocInRead - 1, tailLength);
		}
		else
		{
			int readLength = readInfo_pe2.returnReadSeqLength();
			int substrStartLocInRead = readLength - tailLength + 1;
			return ((readInfo_pe2).returnRcmReadQual()).substr(substrStartLocInRead - 1, tailLength);
		}		
	}


	string returnReadName_alignInfoType(int alignInfoType)
	{
		if((alignInfoType == 1) || (alignInfoType == 2))
		{
			return (readInfo_pe1).returnReadName();
			//readSeqOriginal = (readInfo->readInfo_pe1).readSeq;
		}
		else
		{
			return (readInfo_pe2).returnReadName();
			//readSeqOriginal = (readInfo->readInfo_pe2).readSeq;
		}
		//return readName;
	}

	string returnOriReadSeq_alignInfoType(int alignInfoType)
	{
		if((alignInfoType == 1) || (alignInfoType == 2))
		{
			//readName = (readInfo->readInfo_pe1).readName;
			return (readInfo_pe1).returnReadSeq();
		}
		else
		{
			//readName = (readInfo->readInfo_pe2).readName;
			return (readInfo_pe2).returnReadSeq();
		}		
	}

	string returnReadSeqInDirection_alignInfoType(int alignInfoType)
	{
		if(alignInfoType == 1)
		{
			return (readInfo_pe1).returnReadSeq();
		}
		else if(alignInfoType == 2)
		{
			return (readInfo_pe1).returnRcmReadSeq();
		}
		else if(alignInfoType == 3)
		{
			return (readInfo_pe2).returnReadSeq();
		}
		else
		{
			return (readInfo_pe2).returnRcmReadSeq();
		}
	}


	int returnReadLength(bool End1OrEnd2_bool)
	{
		if(End1OrEnd2_bool)
		{
			return readInfo_pe1.returnReadSeqLength();
		}
		else
		{
			return readInfo_pe2.returnReadSeqLength();
		}	
	}

	int returnReadLength(int tmpAlignInfoType)
	{
		if(tmpAlignInfoType <= 2)
		{
			return (readInfo_pe1).returnReadSeqLength();
		}
		else
		{
			return (readInfo_pe2).returnReadSeqLength();
		}
	}

	PE_Read_Info()
	{}

	bool checkEnd1OrEnd2WithAlignInfoTypeNo(int alignInfoType)
	{
		if(alignInfoType <= 2)
			return true;
		else
			return false;
	}

	bool checkNorOrRcmWithAlignInfoTypeNo(int alignInfoType)
	{
		if((alignInfoType == 1)||(alignInfoType == 3))
			return true;
		else
			return false;
	}

	int checkReadLengthWithAlignInfoTypeNo(int alignInfoType)
	{
		if(alignInfoType <= 2)
			return (readInfo_pe1.returnReadSeqLength());
		else
			return (readInfo_pe2.returnReadSeqLength());
	}

	void get_PE_Read_Info(const string& readName1, const string& readName2,
		const string& readSeq1, const string& readSeq2)
	{
		readInfo_pe1.assignReadName(readName1);
		//readInfo_pe1.readSeqLength = readSeq1.length();
		readInfo_pe1.assignReadSeq(readSeq1);
		readInfo_pe1.assignRcmReadSeq(getRcmSeq(readSeq1));
		//readInfo_pe1.rcmReadSeq = rcmReadSeq1;
		//readInfo_pe1.readPeInfo = "1";

		readInfo_pe2.assignReadName(readName2);
		//readInfo_pe2.readSeqLength = readSeq2.length();
		readInfo_pe2.assignReadSeq(readSeq2);
		readInfo_pe2.assignRcmReadSeq(getRcmSeq(readSeq2));

		readInfo_pe1.assignReadSeqLength((readInfo_pe1.returnReadSeq()).length());
		readInfo_pe2.assignReadSeqLength((readInfo_pe2.returnReadSeq()).length());	

		//readInfo_pe2.rcmReadSeq = rcmReadSeq2;
		//readInfo_pe2.readPeInfo = "2";
	}

	PE_Read_Info(const string& readName1, const string& readName2,
		const string& readSeq1, const string& readSeq2)
	{
		//(readInfo_pe1.readName) = readName1;
		//readInfo_pe1.readSeqLength = readSeq1.length();
		//(readInfo_pe1.readSeq) = readSeq1;
		//readInfo_pe1.rcmReadSeq = rcmReadSeq1;
		//readInfo_pe1.readPeInfo = "1";
		readInfo_pe1.assignReadName(readName1);
		//readInfo_pe1.readSeqLength = readSeq1.length();
		readInfo_pe1.assignReadSeq(readSeq1);

		//(readInfo_pe2.readName) = readName2;
		//readInfo_pe2.readSeqLength = readSeq2.length();
		//(readInfo_pe2.readSeq) = readSeq2;
		//readInfo_pe2.rcmReadSeq = rcmReadSeq2;
		//readInfo_pe2.readPeInfo = "2";
		readInfo_pe2.assignReadName(readName2);
		//readInfo_pe2.readSeqLength = readSeq2.length();
		readInfo_pe2.assignReadSeq(readSeq2);
	}

	void getBothEndRcmReadSeq()
	{
		readInfo_pe1.assignRcmReadSeq(covertStringToReverseComplement(readInfo_pe1.returnReadSeq()));//, (readInfo_pe1.readSeq).length());
		readInfo_pe2.assignRcmReadSeq(covertStringToReverseComplement(readInfo_pe2.returnReadSeq()));//, (readInfo_pe2.readSeq).length());
	}

	void getReverseComplementReadSeq(char* readChar, char* readChar_PE)
	{
		//readInfo_pe1.rcmReadSeq = convertStringToReverseComplement(readInfo_pe1.readSeq);
		//readInfo_pe2.rcmReadSeq = convertStringToReverseComplement(readInfo_pe2.readSeq);		
		readInfo_pe1.assignRcmReadSeq( 
			convertCharArrayToReverseCompletmentStr(readChar, readInfo_pe1.returnReadSeqLength()));
		readInfo_pe2.assignRcmReadSeq( 
			convertCharArrayToReverseCompletmentStr(readChar_PE, readInfo_pe2.returnReadSeqLength()));
	}

	void getRcmReadSeq()
	{
		readInfo_pe1.assignRcmReadSeq(revcomp(readInfo_pe1.returnReadSeq()));
		readInfo_pe2.assignRcmReadSeq(revcomp(readInfo_pe2.returnReadSeq()));
	}

	PE_Read_Info(const string& readName1, const string& readName2,
		const string& readSeq1, const string& rcmReadSeq1, 
		const string& readSeq2, const string& rcmReadSeq2)
	{
		readInfo_pe1.assignReadName(readName1);
		readInfo_pe1.assignReadSeqLength(readSeq1.length());
		readInfo_pe1.assignReadSeq(readSeq1);
		readInfo_pe1.assignRcmReadSeq(rcmReadSeq1);
		//readInfo_pe1.readPeInfo = "1";

		readInfo_pe2.assignReadName(readName2);
		readInfo_pe2.assignReadSeqLength(readSeq2.length());
		readInfo_pe2.assignReadSeq(readSeq2);
		readInfo_pe2.assignRcmReadSeq(rcmReadSeq2);
		//readInfo_pe2.readPeInfo = "2";

	}

	string getIncompleteEndReadSeq(bool End1OrEnd2, bool NorOrRcm)
	{
		if(End1OrEnd2 && NorOrRcm)
		{
			return readInfo_pe2.returnRcmReadSeq();
		}
		else if(End1OrEnd2 && (!NorOrRcm))
		{
			return readInfo_pe2.returnReadSeq();
		}
		else if((!End1OrEnd2) && NorOrRcm)
		{
			return readInfo_pe1.returnRcmReadSeq();
		}
		else
		{
			return readInfo_pe1.returnReadSeq();
		}
	}

	string getReadSeq(bool End1OrEnd2, bool NorOrRcm)
	{
		if(End1OrEnd2 && NorOrRcm)
		{
			return readInfo_pe1.returnReadSeq();
		}
		else if(End1OrEnd2 && (!NorOrRcm))
		{
			return readInfo_pe1.returnRcmReadSeq();
		}
		else if((!End1OrEnd2) && NorOrRcm)
		{
			return readInfo_pe2.returnReadSeq();
		}
		else
		{
			return readInfo_pe2.returnRcmReadSeq();
		}
	}
	void printPEreadInfo()
	{
		cout << "end 1 read: " <<  readInfo_pe1.returnReadName() << endl;
		cout << "length: " << readInfo_pe1.returnReadSeqLength() << endl;
		cout << "oriReadSeq: " << endl << readInfo_pe1.returnReadSeq() << endl;
		cout << "rcmReadSeq: " << endl << readInfo_pe1.returnRcmReadSeq() << endl;

		cout << "end 2 read: " <<  readInfo_pe2.returnReadName() << endl;
		cout << "length: " << readInfo_pe2.returnReadSeqLength() << endl;
		cout << "oriReadSeq: " << endl << readInfo_pe2.returnReadSeq() << endl;		
		cout << "rcmReadSeq: " << endl << readInfo_pe2.returnRcmReadSeq() << endl;
	}

};

#endif

