// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef READ_SEQ_PREPROCESSING_H
#define READ_SEQ_PREPROCESSING_H

#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>

using namespace std;

class InputReadPreProcess
{
public:
	string validBaseCharStr;

	InputReadPreProcess()
	{
		validBaseCharStr = "ATGC";
	}

	string readPreProcess_upperCase_trim(const string& inputReadSeq)
	{
		//cout << "upperCaseSeq: " << endl << this->upperCaseReadSeq(inputReadSeq) << endl;
		//cout << "trimReadSeq: " << endl << this->trimReadSeq(this->upperCaseReadSeq(inputReadSeq)) << endl;
		return ( this->trimReadSeq(this->upperCaseReadSeq(inputReadSeq)) );
	}

	string trimReadSeq(const string& inputReadSeq)
	{
		int validReadSeqStart;
		int validReadSeqEnd;
		
		int inputReadSeqLength = inputReadSeq.length();
		
		validReadSeqStart = inputReadSeq.find_first_of(validBaseCharStr);
		if((validReadSeqStart >= inputReadSeqLength)||(validReadSeqStart < 0))
		{
			return "N";
		}

		validReadSeqEnd = inputReadSeq.find_last_of(validBaseCharStr);

		//cout << "validReadSeqStart: " << validReadSeqStart << endl;
		//cout << "validReadSeqEnd: " << validReadSeqEnd << endl;
		return inputReadSeq.substr(validReadSeqStart, validReadSeqEnd - validReadSeqStart + 1);
	}

	string upperCaseReadSeq(const string& inputReadSeq)
	{
		int stringLength = inputReadSeq.length();
		//cout << "stringLength: " << stringLength << endl;
		string upperCaseSeq;
		if(stringLength < 1)
			return "";
		else
		{

			char ch = inputReadSeq.at(0);
			//cout << "tmp:0 " << " ch: " << ch << endl;
			if((ch == 'A')||(ch == 'a'))
			{
				upperCaseSeq = "A";
			}
			else if((ch == 'C')||(ch == 'c'))
			{
				upperCaseSeq = "C";
			}
			else if((ch == 'G')||(ch == 'g'))
			{
				upperCaseSeq = "G";
			}
			else if((ch == 'T')||(ch == 't'))
			{
				upperCaseSeq = "T";
			}
			else if((ch == 'N')||(ch == 'n'))
			{
				upperCaseSeq = "N";
			}
			else
			{
				cout << "ch: " << ch << endl;
				return "N";
			}
		}

		for(int tmp = 1; tmp < stringLength; tmp ++)
		{
			char ch = inputReadSeq.at(tmp);
			//cout << "tmp: " << tmp << " ch: " << ch << endl;
			if((ch == 'A')||(ch == 'a'))
			{
				upperCaseSeq += "A";
			}
			else if((ch == 'C')||(ch == 'c'))
			{
				upperCaseSeq += "C";
			}
			else if((ch == 'G')||(ch == 'g'))
			{
				upperCaseSeq += "G";
			}
			else if((ch == 'T')||(ch == 't'))
			{
				upperCaseSeq += "T";
			}
			else if((ch == 'N')||(ch == 'n'))
			{
				upperCaseSeq += "N";
			}
			else
			{
				return "N";
			}
		}
		return upperCaseSeq;
	}
};

#endif