/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef COMPRESSOR_2STAGE_H
#define COMPRESSOR_2STAGE_H

#include "defs.h"
#include "queue.h"
#include "qelem.h"
#include "mmgr.h"
#include "qsmodel.h"
#include "rangecod.h"
#include "../libs/zlib.h"
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

typedef pair<int, int> ht_value_t;

struct Compressor2StageParams
{
	CRegisteringQueue<QElem>* elementQueue;
	CMemoryMonitor* memon;
	string archive_name;
	int input_size;
	FILE **rc_file;
	uint32 compression_level;
    uint32 compression_level_threshold;
    uint32 no_seq;
    int32 ref_id_in_list;
    vector<string>* inputFilenameVector;
};

class Compressor2Stage
{
	CRegisteringQueue<QElem>* elementQueue;
	CMemoryMonitor* memon;
	string archive_name;
	int input_size;
    uint32 no_seq;
    int32 ref_id_in_list;
    vector<string>* inputFilenameVector;
    
    uint32 compression_level;
	uint32 compression_level_threshold;
    
	// pliki
#ifdef COMPRESS_DESC
	gzFile_s *out_desc;
#else
	FILE  *out_desc;
#endif
	FILE **rc_file;
	
    
	// hash array
	ht_value_t *ht;
	unsigned long long ht_size, ht_filled, ht_size_mask, ht_threshold;
	
	// for hash array
	int hashing_len;
	int hashing_len_bytes;
	double ht_fill_factor;
	const int ref_match_len = 7;
	const int ref_literal_len = 1;
    
	vector<pair<unsigned char*, int>> file_data;
    
	// const for arithmetic coding
	static const int pos_short_limits = 64;
	static const int RC_CTX_FLAGS_LEN  = 2;
	static const int RC_NO_FLAGS_MODEL     = 3 * 3;
	static const int RC_NO_FLAGS_MODEL_MOD = 3;
	static const int RC_NO_FLAGS		= 3;
	static const int RC_SHIFT_FLAGS	= 15;
	static const int RC_SCALE_FLAGS	= 1 << 9;
#define UPD_FLAG_CTX(ctx, flag)	(((ctx) % RC_NO_FLAGS_MODEL_MOD) * (RC_NO_FLAGS_MODEL / RC_NO_FLAGS_MODEL_MOD) + (flag))
	static const int RC_NO_SYMBOLS_MODEL = 5;
	static const int RC_NO_SYMBOLS     = 64;
	static const int RC_SHIFT_SYMBOLS  = 15;
	static const int RC_SCALE_SYMBOLS  = 1 << 10;
	static const int RC_SHIFT_LEN      = 15;
	static const int RC_SCALE_LEN      = 1 << 12;
	static const int FLAG_MATCH        	= 0;
	static const int FLAG_LITERAL_LITERAL	= 1;
	static const int FLAG_LITERAL_MATCH    = 2;
	static const int FLAG_MATCH_LITERAL	= 0;		// for context calculation only
	static const int FLAG_MATCH_MATCH		= 3;		// for context calculation only
	static const int RC_SHIFT_2M_LEN   = 15;
	static const int RC_SCALE_2M_LEN   = 1 << 9;
	static const int RC_SHIFT_1M_POS   = 15;
	static const int RC_SCALE_1M_POS   = 1 << 10;
	static const int RC_SHIFT_SEQ_POS  = 15;
	static const int RC_SCALE_SEQ_POS  = 1 << 11;
	static const int RC_SHIFT_SEQ_ID   = 15;
	static const int RC_SCALE_SEQ_ID   = 1 << 9;
    
	// variables for arithmetic coding
	qsmodel
    qsm_flags[RC_NO_FLAGS_MODEL],
    qsm_symbols[5],
    qsm_len_m1_pref,
    qsm_len_m1_suf[6],
    qsm_len_m2_pref,
    qsm_len_m2_suf[7],
    qsm_pos_m1_pref,
    qsm_pos_m1_suf[5],
    qsm_seq_id_pref,
    *qsm_seq_id_suf,
    qsm_seq_pos_pref,
    qsm_seq_pos_suf[8];
	static int len_2m_model_sizes[7];
	int ctx_flag;
	int n_seq_ids_suf;
	rangecoder rc;
    
	// const for compression
	static const unsigned int PSEUDO_MATCH_N = 0x7ffffffful;
	static const unsigned int PSEUDO_MATCH_n = 0x7ffffffeul;
    
	// variables for compression
	int expected_match_pos;
	int ctx_symbol;
	vector<tuple<int, int, int>> v_last_seq_pos;		// tuple_pos, seq_pos, byte_pos
    
	// hash array
	void prepare_ht();
	void restruct_ht();
	void ht_insert(int seq_id, int seq_pos);
	unsigned long long hash_fun(unsigned char *seq, int seq_size, int seq_pos);
	void add_to_ht(int seq_id, unsigned char *data, int size);
    
	bool open_files();
    
	// arithmetic coding
	void init_rc_models(int mode);
	void delete_rc_models(void);
    
	// compression
	void store_1st_seq(unsigned char *data, int size);
	void store_literal_1m(unsigned char *s, int &encoded_symbols);
	void decode_match_1l_par(unsigned char *p, int &ref_pos, int &len);
	void store_literal_1l(unsigned char *s, int &encoded_symbols);
	int ctx_from_symbol(unsigned char c);
	void process_data(unsigned char *data, int size);
	void store_match(int seq_id, int seq_pos, int match_len, int &encoded_symbols);
	int symbols_in_match(unsigned char *p, int match_len);
	int tuple_distance(int seq_id, int pos1, int pos2);
	void advance_in_sequence(int seq_id, int tuple_pos, int seq_pos, int symbols_to_advance, pair<int, int> &expected_seq_pos);
	int reduce_match_len(unsigned char *s, int len);
	void determine_expected_pos(unsigned char *p, int match_len);
    
	// decompression
	vector<pair<string, int>> output_names;
    vector<uchar> eol_types;
    vector<char> new_line_at_end_of_files;
    vector<uint32> line_widths;
    vector<uint32*> seq_sizes;
    vector<vector<string>> seq_names;
	void unstore_1st_seq(unsigned char* data, int size, bool ht);
    uint32 cur_file_id;
	uint32 cur_file_pos;
    uint32 cur_file_size;
	void unstore_literal_1m(unsigned char *data, int &encoded_symbols, int &encoded_tuples, bool ht);
	void unstore_literal_1l(unsigned char *data, int &encoded_symbols, int &encoded_tuples, bool ht);
	void add_tuple_seq_pos_mapping(unsigned char *p, int size, int encoded_tuples);
	vector<vector<int>> v_tuple_pos_map; // mapping of tuple positions onto seq_pos positions
	void unprocess_data(unsigned char*&data, int &data_size, bool ht);
	void unstore_match(unsigned char *data, int &encoded_symbols, int &encoded_tuples, bool ht);
    int increase_match_len(unsigned char *s, int len);
	void putDataToQueue(uchar* d, int d_size, string filename, uchar eol_type, char new_line_at_end_of_file, uint32 line_width, uint32 * seq_sizes, vector<string> &_seq_names);
    
public:
	Compressor2Stage(Compressor2StageParams params) : expected_match_pos(0), ctx_symbol(0), cur_file_id(0), cur_file_pos(0), cur_file_size(0)
	{
		hashing_len = 3;
		hashing_len_bytes = 11;
		ht_fill_factor = 0.55;
        
		ctx_symbol = 0;
        
		elementQueue = params.elementQueue;
		memon = params.memon;
		archive_name = params.archive_name;
		input_size = params.input_size;
		rc_file = params.rc_file;
        compression_level_threshold = params.compression_level_threshold;
		compression_level =  params.compression_level;
        
        no_seq = params.no_seq;
        ref_id_in_list = params.ref_id_in_list;
        inputFilenameVector = params.inputFilenameVector;
	}
    
	void start_compression();
	void start_decompression();
};

#endif