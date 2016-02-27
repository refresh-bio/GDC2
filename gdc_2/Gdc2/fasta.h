/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef _FASTA_H
#define _FASTA_H

#include "defs.h"
#include <iostream>
#include <vector>
#include "dmutex.h"
using namespace std;

typedef enum {mode_fa_none, mode_fa_read, mode_fa_write, mode_fa_read_ra} t_fa_mode;

// ********************************************************************************************
class CFastaSequence {
public:
	uchar *raw_data;
	int32 *char_stat;
	uint32 size;
	string seq_name;
	uint32 line_len;
	uint32 start_pos;
	uchar eol_type;

public:
	CFastaSequence();
	CFastaSequence(uchar *_raw_data, string &_seq_name, uint32 _line_len, uint32 _start_pos, uint32 _size, uchar _eol_type);
	CFastaSequence(const CFastaSequence &x);
	~CFastaSequence();
};

// ********************************************************************************************
class CFastaFile {
	FILE *file;
	char *io_buffer;

public:
	DiskMutex* diskMtx;
	int64 file_size;
	t_fa_mode mode;

	uchar *raw_data;
	int64 raw_data_pos;
	int64 trans_data_pos;
	int64 data_size;
	uchar *raw_data_ptr;
	uchar *trans_data_ptr;
	uchar *eof_data_ptr;
    char new_line_at_end_of_file;

	vector<CFastaSequence> sequences;
	int32 char_stat[256];

	bool ReadSequence();
	bool WriteSequenceToMemory(CFastaSequence &seq);
	bool WriteSequenceToFile(CFastaSequence &seq, uchar * raw_data);

public:
	CFastaFile() {
		file     = NULL;
		mode     = mode_fa_none;
		raw_data = NULL;
		io_buffer = NULL;
		diskMtx = NULL;
	};

	~CFastaFile() {
		if(file)
			fclose(file);
		if(raw_data)
			delete[] raw_data;
		if(io_buffer)
			delete[] io_buffer;
	}

	bool Open(string file_name);
	bool Create(string file_name);
	bool Close();
	bool Read(vector<CFastaSequence> &_sequences);
	bool Write(vector<CFastaSequence> &_sequences);
	int64 GetFileSize() {return file_size;}
};

#endif
