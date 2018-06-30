// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>




using namespace std;
#pragma warning(disable:4996)

class Exon
{
public:
	int start;
	int end;

};

class Transcript
{
public:
	string id;
	string chrom;
	string strand;
	vector<Exon> exon;
};

class Junction
{
public:
	string chrom;
	string strand;
	int start;
	int end;

	friend bool operator < (const Junction &ls, const Junction &rs);
};

inline bool operator < (const Junction &ls, const Junction &rs)
{
	if(ls.chrom != rs.chrom)
		return ls.chrom < rs.chrom;
	if(ls.start != rs.start)
		return ls.start < rs.start;
	if(ls.end != rs.end)
		return ls.end < rs.end;
	else
		return ls.strand < rs.strand;
}

class parse_annotation
{
public:
	ifstream input_fs;
	ofstream output_fs;
	map<int, int> donor_set;
	map<int, int> acceptor_set;
	vector<Transcript> transcript_set;
	vector<Junction> junction_set;
	string cur_gene;
	string cur_transcript;
	string cur_chrom;
	string cur_strand;
	int junc_num;

	void parse(char* input_file)
	{
		input_fs.open(input_file);
		if( !input_fs ) 
		{
			fprintf(stderr,"error: open input file error\n");exit(1);
		} 			
		string line;
		cur_gene = "";
		cur_transcript = "";
		cur_chrom = "";
	    cur_strand = "";
		junc_num = 0;
		while(getline(input_fs,line))
		{
			char chrom_tmp[1000], gene_tmp[1000], transcript_tmp[1000], type_tmp[1000], strand_tmp[1000];
			int start, end;
			sscanf(line.c_str(), "%s\t%*s\t%s\t%d\t%d\t%*s\t%s\t%*s\t%[^;];%[^;]", chrom_tmp, type_tmp, &start, &end, strand_tmp, gene_tmp, transcript_tmp);
			string chrom = chrom_tmp;
			string gene = gene_tmp;
			string transcript = transcript_tmp;
			string type = type_tmp;
			string strand = strand_tmp;
			if(type != "exon")
				continue;
			Exon new_exon;
			new_exon.start = start;
			new_exon.end = end;
			if(gene != cur_gene)
			{
				generate_splice_site();
				generate_junction();
				donor_set.clear();
				acceptor_set.clear();
				transcript_set.clear();
				cur_gene = gene;
				cur_chrom = chrom;
				cur_strand = strand;
			}
			if(transcript != cur_transcript)
			{
				cur_transcript = transcript;
				Transcript new_transcript;
				new_transcript.id = transcript;
				new_transcript.chrom = chrom;
				new_transcript.strand = strand;
				transcript_set.push_back(new_transcript);
			}
			transcript_set[transcript_set.size() - 1].exon.push_back(new_exon);
		}
		generate_splice_site();
		generate_junction();
		donor_set.clear();
		acceptor_set.clear();
		transcript_set.clear();
		input_fs.close();
	}

	void generate_splice_site()
	{
		for(size_t i = 0; i < transcript_set.size(); i++)
		{
			if(transcript_set[i].strand == "+")
			{
				for(size_t j = 0; j < transcript_set[i].exon.size(); j++)
				{
					int new_acceptor = transcript_set[i].exon[j].start;
					int new_donor = transcript_set[i].exon[j].end;
					if(j != 0 && acceptor_set.find(new_acceptor) == acceptor_set.end())
						acceptor_set.insert( make_pair(new_acceptor, 1) );
					if(j + 1 < transcript_set[i].exon.size() && donor_set.find(new_donor) == donor_set.end())
						donor_set.insert( make_pair(new_donor, 1) );
				}
			}
			else
			{
				for(size_t j = 0; j < transcript_set[i].exon.size(); j++)
				{
					int new_acceptor = transcript_set[i].exon[j].start;
					int new_donor = transcript_set[i].exon[j].end;
					if(j + 1 < transcript_set[i].exon.size() && acceptor_set.find(new_acceptor) == acceptor_set.end())
						acceptor_set.insert( make_pair(new_acceptor, 1) );
					if(j != 0 && donor_set.find(new_donor) == donor_set.end())
						donor_set.insert( make_pair(new_donor, 1) );
				}
			}
		}
	}

	void generate_junction()
	{
		for(map<int, int>::iterator it1 = donor_set.begin(); it1 != donor_set.end(); it1++)
		{
			for(map<int, int>::iterator it2 = acceptor_set.begin(); it2 != acceptor_set.end(); it2++)
			{
				if(it2->first - it1->first > 1)
				{
					Junction new_junc;
					new_junc.chrom = cur_chrom;
					new_junc.strand = cur_strand;
					new_junc.start = it1->first;
					new_junc.end = it2->first;
					junction_set.push_back(new_junc);
				}
			}
		}
	}

	void output_junction(char* output_file)
	{
		sort(junction_set.begin(), junction_set.end());
		output_fs.open(output_file);
		if( !output_fs ) 
		{
			fprintf(stderr,"error: write sequence file error\n");exit(1);
		}	
		output_fs << "track name=junctions description=\"Mapsplice junctions\"" <<endl;
		for(size_t i = 0; i< junction_set.size(); i++)
		{
			junc_num ++;
			output_fs << junction_set[i].chrom << "\t" << junction_set[i].start << "\t" <<  junction_set[i].end << "\t" << "JUNC_" << junc_num << "\t" << "1" << "\t" << junction_set[i].strand;
			output_fs << "\t" << junction_set[i].start << "\t" << junction_set[i].end << "\t255,0,0\t0\t75,75,\t0,21,\t0\t0\tGTAG\t1\t0.992798\t2\t2\t2"<< endl;
		}
		output_fs.close();
	}
};



void print_usage()
{
	fprintf(stderr,"parse_annotation [input_annotation_file] [output_junction_file]\n");
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc == 1)
		print_usage();
	else if (argc < 3)
	{
		fprintf(stderr,"error: too few arguments\n");
		print_usage();
		exit(1);
	} 
	else if (argc > 3)
	{
		fprintf(stderr,"error: too many arguments\n");
		print_usage();
		exit(1);
	}
	parse_annotation pa;
	pa.parse(argv[1]);
	pa.output_junction(argv[2]);
}
