/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#include "stdafx.h"


#include "defs.h"
#include "fasta.h"
#include <algorithm>

using namespace std;

const uint32 io_buffer_size = 1 << 16;

// ********************************************************************************************
// CFastaSequence
// ********************************************************************************************
CFastaSequence::CFastaSequence()
{
	start_pos = 0;
	size = 0;
	line_len = 0;
	raw_data = NULL;
}

// ********************************************************************************************
CFastaSequence::CFastaSequence(uchar *_raw_data, string &_seq_name, uint32 _line_len, uint32 _start_pos, uint32 _size, uchar _eol_type)
{
	raw_data  = _raw_data;
	seq_name  = _seq_name;
	line_len  = _line_len;
	size      = _size;
	start_pos = _start_pos;
	eol_type  = _eol_type;
}

// ********************************************************************************************
CFastaSequence::CFastaSequence(const CFastaSequence &x)
{
	raw_data  = x.raw_data;
	char_stat = x.char_stat;
	seq_name  = x.seq_name;

	start_pos = x.start_pos;
	size      = x.size;
	line_len  = x.line_len;
	eol_type  = x.eol_type;
}

// ********************************************************************************************
CFastaSequence::~CFastaSequence()
{
}


// ********************************************************************************************
// CFastaFile
// ********************************************************************************************
bool CFastaFile::Open(string file_name)
{
	if(mode != mode_fa_none)
		return false;

	if((file = my_fopen(file_name.c_str(), "rb")) != NULL)
		mode = mode_fa_read;

	if(mode != mode_fa_read)
		return false;

	my_fseek(file, 0, SEEK_END);
	file_size = my_ftell(file);
	my_fseek(file, 0, SEEK_SET);

	fill_n(char_stat, 256, 0);

	return true;
}

// ********************************************************************************************
bool CFastaFile::Create(string file_name)
{
	if((file = my_fopen(file_name.c_str(), "wb")) != NULL)
		mode = mode_fa_write;
	file_size = 0;

	io_buffer = new char[io_buffer_size];
	setvbuf(file, io_buffer, _IOFBF, io_buffer_size);

	return mode == mode_fa_write;
}

// ********************************************************************************************
bool CFastaFile::Close()
{
	if(!file)
		return false;

	fclose(file);
	file = NULL;
	mode = mode_fa_none;

	if(raw_data)
		delete[] raw_data;
	raw_data = NULL;

	return true;
}

// ********************************************************************************************
bool CFastaFile::ReadSequence()
{
	string seq_name;
	int32 line_len;
	int32 i;
	uchar c;
	uchar *start_ptr = trans_data_ptr;
	uchar eol_type;

	if(mode != mode_fa_read)
		return false;
	
	line_len = 0;

	c = *raw_data_ptr++;

	if(c != '>' && c != '@')
	{
		return false;
	}
	if(raw_data_ptr >= eof_data_ptr)
		return false;

	seq_name = "";
	while(true)
	{
		c = *raw_data_ptr++;
		if(c == '\n' || c == '\r')
			break;
		seq_name.push_back(c);
	}

	eol_type = c;
	if(*raw_data_ptr == '\n' || *raw_data_ptr == '\r')
	{
		++raw_data_ptr;
		eol_type = 255;				// 0x0d + 0x0a
	}

	uchar *tmp_data_ptr = raw_data_ptr;
	for(i = 0; ; ++i)
	{
		c = *tmp_data_ptr++;
		if(c == '>' || c == '\n' || c == '\r' || c == '@')
			break;
	}
	line_len = i;

	while(true)
	{
		c = *raw_data_ptr++;
		if(c >= 'A')
		{
			*trans_data_ptr++ = c;		
//			char_stat[c]++;
		}
		else if(c == '>' || c == '@')
		{
			raw_data_ptr--;
			break;
		}
	}

	sequences.push_back(CFastaSequence(NULL, seq_name, line_len, (uint32) (start_ptr - raw_data), (uint32) (trans_data_ptr - start_ptr), eol_type));

	return true;
}

// ********************************************************************************************
bool CFastaFile::Read(vector<CFastaSequence> &_sequences)
{
	raw_data = new uchar[file_size+1];
	raw_data_pos   = 0;
	trans_data_pos = 0;

	raw_data_ptr = raw_data;
	trans_data_ptr = raw_data;
	eof_data_ptr = raw_data+file_size;

	if(diskMtx)
		diskMtx->lock();
	fread(raw_data, 1, file_size, file);
    if(raw_data[file_size-1] == '\n' || raw_data[file_size-1]  == '\r')
        new_line_at_end_of_file = 1;
    else
        new_line_at_end_of_file = 0;
        
	if(diskMtx)
		diskMtx->unlock();
	raw_data[file_size] = '>';

	sequences.clear();
	while(ReadSequence())
		;

	for(size_t i = 0; i < sequences.size(); ++i)
	{
		sequences[i].raw_data  = raw_data;
		sequences[i].char_stat = char_stat;
	}

	_sequences = sequences;

	data_size = sequences.back().start_pos + sequences.back().size;

	return true;
}

// ********************************************************************************************
bool CFastaFile::WriteSequenceToMemory(CFastaSequence &seq)
{
	if(mode != mode_fa_write)
		return false;
	uchar eol = seq.eol_type;

	*raw_data_ptr++ = '>';
	copy_n(seq.seq_name.c_str(), seq.seq_name.length(), raw_data_ptr);
	raw_data_ptr += seq.seq_name.length();

	if(eol == 255)
	{
		*raw_data_ptr++ = 0x0d;
		*raw_data_ptr++ = 0x0a;
	}
	else
		*raw_data_ptr++ = eol;
	
	uchar *seq_start = seq.raw_data + seq.start_pos;

	int len = seq.line_len;
	int n_lines = ((int32) seq.size) / len;
	uint32 s_pos = 0;

	if(eol == 0xff)
	{
		for(int32 i = 0; i < n_lines; ++i, s_pos += len)
		{
			memcpy(raw_data_ptr+s_pos+2*i, seq_start+s_pos, len);
			raw_data_ptr[s_pos+2*i+len] = 0x0d;
			raw_data_ptr[s_pos+2*i+len+1] = 0x0a;
		}
	}
	else
	{
		for(int32 i = 0; i < n_lines; ++i, s_pos += len)
		{
			memcpy(raw_data_ptr+s_pos+i, seq_start+s_pos, len);
			raw_data_ptr[s_pos+i+len] = eol;
		}
	}

	raw_data_ptr += s_pos + n_lines;
	if(eol == 0xff)
		raw_data_ptr += n_lines;
	if(s_pos < seq.size)
	{
		// ???
		memcpy(raw_data_ptr, seq_start+s_pos, seq.size-s_pos);
		if(eol == 0xff)
		{
			raw_data_ptr[seq.size-s_pos] = 0x0d;
			raw_data_ptr[seq.size-s_pos+1] = 0x0a;
		}
		else
			raw_data_ptr[seq.size-s_pos] = eol;
		raw_data_ptr += seq.size-s_pos+1 + (eol == 0xff);
	}

	return true;
}

// ********************************************************************************************
bool CFastaFile::WriteSequenceToFile(CFastaSequence &seq, uchar * raw_data)
{
	if(mode != mode_fa_write)
		return false;
	uchar eol = seq.eol_type;

	putc('>', file);
	fwrite(seq.seq_name.c_str(), 1, seq.seq_name.length(), file);

	if(eol == 255)
	{
		putc(0x0d, file);
		putc(0x0a, file);
	}
	else
		putc(eol, file);
	
	uchar *seq_start = raw_data + seq.start_pos;

	int len = seq.line_len;
	int n_lines = ((int32) seq.size) / len;
	uint32 s_pos = 0;

	if(eol == 0xff)
	{
		for(int32 i = 0; i < n_lines; ++i, s_pos += len)
		{
			fwrite(seq_start+s_pos, 1, len, file);
			putc(0x0d, file);
			putc(0x0a, file);
		}
	}
	else
	{
		for(int32 i = 0; i < n_lines; ++i, s_pos += len)
		{
			fwrite(seq_start+s_pos, 1, len, file);
			putc(eol, file);
		}
	}

	if(s_pos < seq.size)
	{
		// ???
		fwrite(seq_start+s_pos, 1, seq.size-s_pos, file);
		if(eol == 0xff)
		{
			putc(0x0d, file);
			putc(0x0a, file);
		}
		else
			putc(eol, file);
	}

	return true;
}

// ********************************************************************************************
bool CFastaFile::Write(vector<CFastaSequence> &_sequences)
{
	if(mode != mode_fa_write)
		return false;

	uint32 tot_size = 0;
	for(vector<CFastaSequence>::iterator p = _sequences.begin(); p != _sequences.end(); ++p)
	{
		tot_size += p->size;															// data
		tot_size += (2) * (p->size / p->line_len + 1);			// EOLs
		tot_size += (uint32) p->seq_name.length() + 3;								// seq. name
	}

    data_size = 0;
    uchar * raw_data = _sequences[0].raw_data;

	for(vector<CFastaSequence>::iterator p = _sequences.begin(); p != _sequences.end(); ++p)
		WriteSequenceToFile(*p, raw_data);

	raw_data = NULL;

	return true;
}
