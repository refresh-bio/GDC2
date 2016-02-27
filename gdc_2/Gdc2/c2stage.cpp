/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#include "c2stage.h"

int Compressor2Stage::len_2m_model_sizes[7] = {16, 32, 128, 256, 256, 256, 256};

void Compressor2Stage::start_compression()
{
	prepare_ht();
	open_files();
    //store compression level at the beginning of archive description in 1 byte
    char compr_lev = (char)compression_level;
    my_fwrite(&compr_lev, 1, 1, out_desc);
    
	init_rc_models(1);
	start_encoding(&rc, 0, 0);
	uint32 qelem_no = 0;
	while (!elementQueue->IsCompleted())
	{
		unsigned char *data;
		int size;
		QElem e;
		if (elementQueue->Pop(e))
		{
			data = e.data;
			size = e.dataSize;
            
			cout << "Processing " << e.filename << " (" << qelem_no+1 << " of " << input_size << ")\n";
            
			if(qelem_no > 0)		// not first file
				process_data(data, size);
			else
				store_1st_seq(data, size);		// store 1st sequence without 2nd level compression
            
			if(++qelem_no < compression_level_threshold) // memory limiters
			{
                file_data.push_back(make_pair(data, size));
                add_to_ht((int) file_data.size()-1, data, size);
			}
			else
			{
				delete data;
				data = NULL;
			}
			// Store files description
            my_fwrite(e.filename.c_str(), 1, e.filename.length(), out_desc);
			my_putc(0, out_desc);
			my_fwrite(&size, 4, 1, out_desc);
            my_fwrite(&e.eol_type, 1, 1, out_desc);
            my_fwrite(&e.new_line_at_end_of_file, 1, 1, out_desc);
            my_fwrite(&e.line_width, 4, 1, out_desc);
            
            for(uint32 i = 0; i < no_seq; i++)
            {
                my_fwrite(e.seq_names[i].c_str() ,  1, e.seq_names[i].length(), out_desc);
                my_putc(0, out_desc);
                my_fwrite(&e.seq_sizes[i], 4, 1, out_desc);
            }
            my_putc(1, out_desc);  // end indicator
            
            
            delete e.seq_sizes;
            
		}
	}
	done_encoding(&rc);
	delete_rc_models();
    fclose(*rc_file);
	my_fclose(out_desc);

}

bool Compressor2Stage::open_files()
{
#ifdef WIN32
	errno_t err;
	err = fopen_s(&*rc_file, (archive_name+".gdc2_rc").c_str(), "wb");
	if (err != 0)
		return false;
 #ifdef COMPRESS_DESC
	out_desc = gzopen((archive_name + ".gdc2_desc").c_str(), "wb");
 #else
	err = fopen_s(&out_desc, (archive_name+".gdc2_desc").c_str(), "wb");
	if (err != 0)
		return false;
 #endif

#else
	*rc_file = fopen((archive_name+".gdc2_rc").c_str(), "wb");
	if (*rc_file == NULL)
		return false;
 #ifdef COMPRESS_DESC
	out_desc = gzopen((archive_name + ".gdc2_desc").c_str(), "wb");
 #else
	out_desc = fopen((archive_name+".gdc2_desc").c_str(), "wb");
 #endif
	if (out_desc == NULL)
		return false;

#endif
	return true;
}

void Compressor2Stage::prepare_ht()
{
	ht_size = 0;
	ht_filled = 0;
	ht_threshold = 0;
	ht = nullptr;
	restruct_ht();
}

void Compressor2Stage::restruct_ht()
{
	unsigned long long ht_size_old = ht_size;
	if(ht_size_old == 0)
		ht_size = 1 << 25;
	else
	{
		// Estimate the final size of hash table
		ht_size = ht_size_old * ((double) input_size / file_data.size()) * 1.1;
        ht_size = ((double) compression_level/10) * ht_size;
       // Round to the nearest larger power of 2
		while(ht_size & (ht_size - 1))
			ht_size &= ht_size - 1;
		ht_size <<= 1;
	}
    
	ht_filled = 0;
	ht_size_mask = ht_size - 1;

	ht_threshold = (unsigned long long) (ht_size * ht_fill_factor);
    
	ht_value_t *ht_old = ht;
	ht = new ht_value_t[ht_size];
    
	for(uint64 i = 0; i < ht_size; ++i)
		ht[i].first = -1;
    
	for(uint64 i = 0; i < ht_size_old; ++i)
		if(ht_old[i].first >= 0)
			ht_insert(ht_old[i].first, ht_old[i].second);
    
	delete[] ht_old;
}

void Compressor2Stage::ht_insert(int seq_id, int seq_pos)
{
	if(ht_filled > ht_threshold)
		restruct_ht();
    
    unsigned char *seq = file_data[seq_id].first;
	int seq_size = file_data[seq_id].second;
    unsigned long long pos = hash_fun(seq, seq_size, seq_pos);
    
    while(true)
	{
		if(ht[pos].first < 0)
		{
			ht[pos].first  = seq_id;
			ht[pos].second = seq_pos;
			break;
		}
		pos = (pos == ht_size_mask) ? 0 : pos+1;
	}
	ht_filled++;
}

unsigned long long Compressor2Stage::hash_fun(unsigned char *seq, int seq_size, int seq_pos)
{
	unsigned long long val = 0;
    
	for(int i = 0; i < hashing_len_bytes; )
	{
		if(seq_pos >= seq_size)
			break;
        
		if(seq[seq_pos] < 128)		// match
		{
			for(int j = 0; j < ref_match_len; ++j)
				val = val * 127 + seq[seq_pos++];
			i += ref_match_len;
		}
		else						// literal
		{
			for(int j = 0; j < ref_literal_len; ++j)
				val = val * 63 + seq[seq_pos++];
			i += ref_literal_len;
		}
	}
    
	return val & ht_size_mask;
}

void Compressor2Stage::init_rc_models(int mode)
{
	// Flags
	for(int i = 0; i < RC_NO_FLAGS_MODEL; ++i)
		initqsmodel(&qsm_flags[i], RC_NO_FLAGS, RC_SHIFT_FLAGS, RC_SCALE_FLAGS, NULL, mode);
	ctx_flag = 0;
    
	// Literals
	for(int i = 0; i < RC_NO_SYMBOLS_MODEL; ++i)
		initqsmodel(&qsm_symbols[i], RC_NO_SYMBOLS, RC_SHIFT_SYMBOLS, RC_SCALE_SYMBOLS, NULL, mode);
    
	// Lens in 1st level matches
	initqsmodel(&qsm_len_m1_pref, 3, RC_SHIFT_LEN, RC_SCALE_LEN, NULL, mode);
	for(int i = 0; i < 6; ++i)
		initqsmodel(&qsm_len_m1_suf[i], 256, RC_SHIFT_LEN, RC_SCALE_LEN, NULL, mode);
    
	// Positions in 1st level matches
	initqsmodel(&qsm_pos_m1_pref, 3, RC_SHIFT_1M_POS, RC_SCALE_1M_POS, NULL, mode);
	initqsmodel(&qsm_pos_m1_suf[0], 2*pos_short_limits-2 + 2, RC_SHIFT_1M_POS, RC_SCALE_1M_POS, NULL, mode); // + 2 for pseudomatches
	for(int i = 1; i < 5; ++i)
		initqsmodel(&qsm_pos_m1_suf[i], 256, RC_SHIFT_1M_POS, RC_SCALE_1M_POS, NULL, mode);
    
	// Seq. ids in 2nd level matches
	n_seq_ids_suf = ((int) input_size - 1 + 256) / 256 + 1;
	initqsmodel(&qsm_seq_id_pref, n_seq_ids_suf, RC_SHIFT_SEQ_ID, RC_SCALE_SEQ_ID, NULL, mode);
	qsm_seq_id_suf = new qsmodel[n_seq_ids_suf];
	// !!! To do wyciecia, bo nie jest uzywane
	initqsmodel(&qsm_seq_id_suf[0], 1, RC_SHIFT_SEQ_ID, RC_SCALE_SEQ_ID, NULL, mode);
	for(int i = 1; i < n_seq_ids_suf; ++i)
		initqsmodel(&qsm_seq_id_suf[i], 256, RC_SHIFT_SEQ_ID, RC_SCALE_SEQ_ID, NULL, mode);
    
	// Seq. pos in 2nd level matches
	initqsmodel(&qsm_seq_pos_pref, 7, RC_SHIFT_SEQ_POS, RC_SCALE_SEQ_POS, NULL, mode);
	for(int i = 0; i < 2; ++i)
		initqsmodel(&qsm_seq_pos_suf[i], 15, RC_SHIFT_SEQ_POS, RC_SCALE_SEQ_POS, NULL, mode);
	for(int i = 2; i < 4; ++i)
		initqsmodel(&qsm_seq_pos_suf[i], 256-16, RC_SHIFT_SEQ_POS, RC_SCALE_SEQ_POS, NULL, mode);
	for(int i = 4; i < 8; ++i)
		initqsmodel(&qsm_seq_pos_suf[i], 256, RC_SHIFT_SEQ_POS, RC_SCALE_SEQ_POS, NULL, mode);
    
	// Lens in 2nd level matches
	initqsmodel(&qsm_len_m2_pref, 5, RC_SHIFT_LEN, RC_SCALE_2M_LEN, NULL, mode);
	for(int i = 0; i < 7; ++i)
		initqsmodel(&qsm_len_m2_suf[i], len_2m_model_sizes[i], RC_SHIFT_2M_LEN, RC_SCALE_2M_LEN, NULL, mode);
}

void Compressor2Stage::store_1st_seq(unsigned char *data, int size)
{
	int tmp = 0;
    
	for(int pos = 0; pos < size; )
	{
		if(data[pos] < 128)			// match
		{
			store_literal_1m(data+pos, tmp);
			pos += ref_match_len;
		}
		else						// literal
		{
			store_literal_1l(data+pos, tmp);
			pos += ref_literal_len;
		}
	}
}

void Compressor2Stage::store_literal_1m(unsigned char *s, int &encoded_symbols)
{
	int syfreq, ltfreq;
    
	// Flag
	qsgetfreq(&qsm_flags[ctx_flag], FLAG_LITERAL_MATCH, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_FLAGS);
	qsupdate(&qsm_flags[ctx_flag], FLAG_LITERAL_MATCH);
	ctx_flag = UPD_FLAG_CTX(ctx_flag, FLAG_LITERAL_MATCH);
    
	int len = s[4];
	len = (len << 8) + s[5];
	len = (len << 8) + s[6];
	int len_to_code;
    
	int pref, ctx_suf;
	if(len < 0x100)
	{
		pref = 0;
		ctx_suf = 0;
		len_to_code = len;
	}
	else if(len < 0x10000 + 0x100)
	{
		len_to_code = len - 0x100;
		pref = 1;
		ctx_suf = 1;
	}
	else
	{
		len_to_code = len - 0x10000 - 0x100;
		pref = 2;
		ctx_suf = 3;
	}
    
	// Len prefix
	qsgetfreq(&qsm_len_m1_pref, pref, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_LEN);
	qsupdate(&qsm_len_m1_pref, pref);
    
	// Len suffix
	for(int i = 0; i <= pref; ++i)
	{
		int b = (len_to_code >> (i * 8)) & 0xff;
		qsgetfreq(&qsm_len_m1_suf[ctx_suf+i], b , &syfreq, &ltfreq);
		encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_LEN);
		qsupdate(&qsm_len_m1_suf[ctx_suf+i], b);
	}
	
	int ref_pos;
	decode_match_1l_par(s, ref_pos, len);
    
	int pos_to_code;
	if(ref_pos != PSEUDO_MATCH_N && ref_pos != PSEUDO_MATCH_n)
		pos_to_code = expected_match_pos - ref_pos;
	else
		pos_to_code = ref_pos;
    
	if(pos_to_code == 0)
		pref = 0;
	else if(pos_to_code > -pos_short_limits && pos_to_code < pos_short_limits)
	{
		pref = 1;
		if(pos_to_code > 0)
			pos_to_code--;
		else
			pos_to_code += 2 * pos_short_limits - 2;
	}
	else
	{
		if(ref_pos == PSEUDO_MATCH_N)
		{
			pos_to_code = 2 * pos_short_limits - 2;
			pref = 1;
		}
		else if(ref_pos == PSEUDO_MATCH_n)
		{
			pos_to_code = 2 * pos_short_limits - 1;
			pref = 1;
		}
		else
			pref = 2;
	}
    
	qsgetfreq(&qsm_pos_m1_pref, pref, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_1M_POS);
	qsupdate(&qsm_pos_m1_pref, pref);
	
	if(pref == 1)
	{
		qsgetfreq(&qsm_pos_m1_suf[0], pos_to_code, &syfreq, &ltfreq);
		encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_1M_POS);
		qsupdate(&qsm_pos_m1_suf[0], pos_to_code);
	}
	else if(pref == 2)
		for(int i = 0; i < 4; ++i)
		{
			int b = (((unsigned int) pos_to_code) >> (8*i)) & 0xff;
			qsgetfreq(&qsm_pos_m1_suf[i+1], b, &syfreq, &ltfreq);
			encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_1M_POS);
			qsupdate(&qsm_pos_m1_suf[i+1], b);
		}
	
	if(ref_pos != PSEUDO_MATCH_N && ref_pos != PSEUDO_MATCH_n)
		expected_match_pos = ref_pos + len;
	else
		expected_match_pos += len;
    
	encoded_symbols += len;
}

void Compressor2Stage::decode_match_1l_par(unsigned char *p, int &ref_pos, int &len)
{
	ref_pos = (((int) p[0]) << 24) + (((int) p[1]) << 16) + (((int) p[2]) << 8) + ((int) p[3]);
	len     = (((int) p[4]) << 16) + (((int) p[5]) << 8) + ((int) p[6]);
}

void Compressor2Stage::store_literal_1l(unsigned char *s, int &encoded_symbols)
{
	int syfreq, ltfreq;
    
	// Flag
	qsgetfreq(&qsm_flags[ctx_flag], FLAG_LITERAL_LITERAL, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_FLAGS);
	qsupdate(&qsm_flags[ctx_flag], FLAG_LITERAL_LITERAL);
	ctx_flag = UPD_FLAG_CTX(ctx_flag, FLAG_LITERAL_LITERAL);
    
	// Symbol
	qsgetfreq(&qsm_symbols[ctx_symbol], s[0] - 0x80 - 0x40, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SYMBOLS);
	qsupdate(&qsm_symbols[ctx_symbol], s[0] - 0x80 - 0x40);
	
	ctx_symbol = ctx_from_symbol(s[0] - 0x80);
	expected_match_pos++;
	encoded_symbols++;
}

int Compressor2Stage::ctx_from_symbol(unsigned char c)
{
	c |= 32;
	if(c == 'a')
		return 1;
	if(c == 'c')
		return 2;
	if(c == 'g')
		return 3;
	if(c == 't')
		return 4;
	return 0;
}

void Compressor2Stage::process_data(unsigned char *data, int size)
{
	ctx_flag = 0;
	expected_match_pos = 0;
    
	v_last_seq_pos.clear();
	v_last_seq_pos.resize(input_size);
    
	int encoded_symbols = 0;
    
	for(int i = 0; i < size; )
	{
		unsigned long long ht_pos = hash_fun(data, size, i);
        //		unsigned long long ht_step = 0;
		int best_len = 0;
		int best_byte_len = 0;
		unsigned long long best_ht_pos = 0;

		// Find longest match
		while(true)
		{
			if(ht[ht_pos].first == -1)
				break;
            
			int ref_seq_id		   = ht[ht_pos].first;
			int ref_seq_pos		   = ht[ht_pos].second;
			int ref_seq_size	   = file_data[ref_seq_id].second;
			unsigned char *ref_seq = file_data[ref_seq_id].first;
            
			int len, byte_len;
			unsigned char *p1 = &ref_seq[ref_seq_pos];
			unsigned char *p2 = &data[i];

			for(len = byte_len = 0; ref_seq_pos+byte_len < ref_seq_size && i+byte_len < size; )
			{
				if(*p1 != *p2)
					break;
             
                if(*p2 < 128)				// match
				{
					p1++;	
					p2++;
					if((*p1++ != *p2++) || (*p1++ != *p2++) || (*p1++ != *p2++) || (*p1++ != *p2++) || (*p1++ != *p2++) || (*p1++ != *p2++))
						break;
					byte_len += ref_match_len;
					len++;
				}
				else									// literal
				{
					p1++;	
					p2++;
					byte_len += ref_literal_len;
					len++;
				}
			}
            
			if(byte_len > best_byte_len)
			{
				best_len      = len;
				best_byte_len = byte_len;
				best_ht_pos   = ht_pos;
			}
            
            ht_pos = (ht_pos == ht_size_mask) ? 0: ht_pos+1;
        }
        
		// Encode what was found
		if(best_byte_len >= hashing_len_bytes)	// 2nd level match found
		{
			store_match(ht[best_ht_pos].first, ht[best_ht_pos].second, best_len, encoded_symbols);
			i += best_byte_len;
		}
		else								// 2nd level literal
		{
			if(data[i] < 128)				// 1st level match
			{
				store_literal_1m(&data[i], encoded_symbols);
				i += ref_match_len;
			}
			else							// 1st level literal
			{
				store_literal_1l(&data[i], encoded_symbols);
				i += ref_literal_len;
			}
		}
	}
}

void Compressor2Stage::store_match(int seq_id, int seq_pos, int match_len, int &encoded_symbols)
{
	int syfreq, ltfreq;
    
	// Flag
	qsgetfreq(&qsm_flags[ctx_flag], FLAG_MATCH, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_FLAGS);
	qsupdate(&qsm_flags[ctx_flag], FLAG_MATCH);
	ctx_flag = UPD_FLAG_CTX(ctx_flag, FLAG_MATCH);
    
	int seq_id_pref, seq_id_suf;
    
	// Direct
	seq_id_pref = (seq_id >> 8) + 1;
	seq_id_suf  = seq_id & 0xff;
    
	// seq id prefix & suffix
	qsgetfreq(&qsm_seq_id_pref, seq_id_pref, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SEQ_ID);
	qsupdate(&qsm_seq_id_pref, seq_id_pref);
    
	qsgetfreq(&qsm_seq_id_suf[seq_id_pref], seq_id_suf, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SEQ_ID);
	qsupdate(&qsm_seq_id_suf[seq_id_pref], seq_id_suf);
	
	pair<int, int> expected_seq_pos;
	advance_in_sequence(seq_id, get<0>(v_last_seq_pos[seq_id]), get<1>(v_last_seq_pos[seq_id]), encoded_symbols-get<2>(v_last_seq_pos[seq_id]), expected_seq_pos);
    
	int to_code_seq_pos = tuple_distance(seq_id, expected_seq_pos.second, seq_pos);
    
	// Seq pos
	int sp_pref, sp_ctx_suf;
	int sp;
    
	if(to_code_seq_pos == 0)
	{
		sp_pref = 0;
	}
	else if(to_code_seq_pos > 0 && to_code_seq_pos < 16)
	{
		sp = to_code_seq_pos-1;		sp_pref = 1;		sp_ctx_suf = 0;
	}
	else if(to_code_seq_pos > -16 && to_code_seq_pos < 0)
	{
		sp = -to_code_seq_pos - 1;		sp_pref = 2;		sp_ctx_suf = 1;
	}
	else if(to_code_seq_pos >= 16 && to_code_seq_pos < 256)
	{
		sp = to_code_seq_pos-16;		sp_pref = 3;		sp_ctx_suf = 2;
	}
	else if(to_code_seq_pos > -256 && to_code_seq_pos <= -16)
	{
		sp = -to_code_seq_pos-16;		sp_pref = 4;		sp_ctx_suf = 3;
	}
	else if(to_code_seq_pos >= 256)
	{
		sp = to_code_seq_pos;		sp_pref = 5;		sp_ctx_suf = 4;
	}
	else
	{
		sp = -to_code_seq_pos;		sp_pref = 6;		sp_ctx_suf = 4;
	}
    
	// seq. pos. prefix
	qsgetfreq(&qsm_seq_pos_pref, sp_pref, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SEQ_POS);
	qsupdate(&qsm_seq_pos_pref, sp_pref);
    
	// Seq. pos. suffix
	if(sp_pref > 0 && sp_pref < 5)
	{
		qsgetfreq(&qsm_seq_pos_suf[sp_ctx_suf], sp, &syfreq, &ltfreq);
		encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SEQ_POS);
		qsupdate(&qsm_seq_pos_suf[sp_ctx_suf], sp);
	}
	else if(sp_pref >= 5)
		for(int i = 0; i < 4; ++i)
		{
			int b = (sp >> (i * 8)) & 0xff;
			qsgetfreq(&qsm_seq_pos_suf[sp_ctx_suf+i], b, &syfreq, &ltfreq);
			encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SEQ_POS);
			qsupdate(&qsm_seq_pos_suf[sp_ctx_suf+i], b);
		}
    
	get<0>(v_last_seq_pos[seq_id]) = 0;					// ignore this value in compression (only necessary in decompression)
	get<1>(v_last_seq_pos[seq_id]) = seq_pos;
	get<2>(v_last_seq_pos[seq_id]) = encoded_symbols;
    
	// Match len
	int ml_pref, ml_ctx_suf;
	unsigned char *ref_seq = file_data[seq_id].first;
	int ml = reduce_match_len(ref_seq+seq_pos, match_len);
    
	if(ml < len_2m_model_sizes[0])
	{
		ml_pref = 0;		ml_ctx_suf = 0;
	}
	else if(ml < len_2m_model_sizes[0]+len_2m_model_sizes[1])
	{
		ml -= len_2m_model_sizes[0];		ml_pref = 1;		ml_ctx_suf = 1;
	}
	else if(ml < len_2m_model_sizes[0]+len_2m_model_sizes[1]+len_2m_model_sizes[2])
	{
		ml -= len_2m_model_sizes[0]+len_2m_model_sizes[1];		ml_pref = 2;		ml_ctx_suf = 2;
	}
	else if(ml < len_2m_model_sizes[0]+len_2m_model_sizes[1]+len_2m_model_sizes[2]+len_2m_model_sizes[3])
	{
		ml -= len_2m_model_sizes[0]+len_2m_model_sizes[1]+len_2m_model_sizes[2];		ml_pref = 3;		ml_ctx_suf = 3;
	}
	else
	{
		ml -= len_2m_model_sizes[0]+len_2m_model_sizes[1]+len_2m_model_sizes[2]+len_2m_model_sizes[3];		ml_pref = 4;		ml_ctx_suf = 4;
	}
    
	// Len prefix
	qsgetfreq(&qsm_len_m2_pref, ml_pref, &syfreq, &ltfreq);
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_LEN);
	qsupdate(&qsm_len_m2_pref, ml_pref);
	
	// Len suffix
	if(ml_pref < 4)
	{
		qsgetfreq(&qsm_len_m2_suf[ml_ctx_suf], ml , &syfreq, &ltfreq);
		encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_LEN);
		qsupdate(&qsm_len_m2_suf[ml_ctx_suf], ml);
	}
	else
		for(int i = 0; i < 3; ++i)
		{
			int b = (ml >> (i * 8)) & 0xff;
			qsgetfreq(&qsm_len_m2_suf[ml_ctx_suf+i], b , &syfreq, &ltfreq);
			encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_LEN);
			qsupdate(&qsm_len_m2_suf[ml_ctx_suf+i], b);
		}
	
	determine_expected_pos(file_data[seq_id].first + seq_pos, match_len);
    
	// Advance encoded symbols
	encoded_symbols += symbols_in_match(file_data[seq_id].first + seq_pos, match_len);
}

// Move symbols_to_advance symbols in encoded seq_id sequence from  seq_pos
// Calculate closest seq_pos from end pos
void Compressor2Stage::advance_in_sequence(int seq_id, int tuple_pos, int seq_pos, int symbols_to_advance, pair<int, int> &expected_seq_pos)
{
	int r = seq_pos;
    
	unsigned char *p = file_data[seq_id].first + seq_pos;
	int seq_size = file_data[seq_id].second;
    
	int ref_pos;
	int len;
        
	while(symbols_to_advance > 0 && r < seq_size)
	{
        if(p[0] < 128)				// 1st level match
		{
			decode_match_1l_par(p, ref_pos, len);
			if(len <= 2*symbols_to_advance)					// find closest tuple
			{
				r += ref_match_len;		p += ref_match_len;
				++tuple_pos;
			}
			symbols_to_advance -= len;
		}
		else
		{
			r += ref_literal_len;		p += ref_literal_len;
			--symbols_to_advance;
			++tuple_pos;
		}
	}
    
	expected_seq_pos.first  = tuple_pos;
	expected_seq_pos.second = r;
}

//**********************************************************
// Calculates distance in steps between pos1 and pos2 in sequence seq_id
int Compressor2Stage::tuple_distance(int seq_id, int pos1, int pos2)
{
	int sign = 1;
    
	if(pos1 > pos2)
	{
		swap(pos1, pos2);
		sign = -1;
	}
    
	int dist = 0;
	unsigned char *p = file_data[seq_id].first;
    
	while(pos1 < pos2)
	{
		++dist;
		if(p[pos1] < 128)					// 1st level match
			pos1 += ref_match_len;
		else
			pos1 += ref_literal_len;
	}
    
	return dist * sign;
}

//**********************************************************
// Calculates number of symbols (real) in match
int Compressor2Stage::symbols_in_match(unsigned char *p, int match_len)
{
	int r = 0;
    
	for(int i = 0; i < match_len; ++i)
	{
		if(p[0] < 128)		// 1st level match
		{
			int ref_pos;
			int len;
			decode_match_1l_par(p, ref_pos, len);
			r += len;
			p += ref_match_len;
		}
		else				// 1st level literal
		{
			++r;
			p += ref_literal_len;
		}
	}
    
	return r;
}

//Checking how many tuples from s, but omitting the ones present in hashing_len_bytes
int Compressor2Stage::reduce_match_len(unsigned char *s, int len)
{
	int pos = 0;
    
	while(pos < hashing_len_bytes)
	{
		if(s[pos] < 128)			// match
			pos += ref_match_len;
		else
			pos += ref_literal_len;
		--len;
	}
    
	return len;
}

void Compressor2Stage::determine_expected_pos(unsigned char *p, int match_len)
{
	for(int i = 0; i < match_len; ++i)
	{
		if(p[0] < 128)	// match 1st level
		{
			int ref_pos;
			int len;
			decode_match_1l_par(p, ref_pos, len);
			expected_match_pos = ref_pos + len;
			p += ref_match_len;
		}
		else			// literal 1st level;
		{
			p += ref_literal_len;
			expected_match_pos++;
		}
	}
}

void Compressor2Stage::add_to_ht(int seq_id, unsigned char *data, int size)
{
	for(int i = 0; i < size;)
	{
		ht_insert(seq_id, i);
		if(data[i] < 128)			// 1st level match
			i += ref_match_len;
		else						// 1st level literal
			i += ref_literal_len;
	}
}

// Delete range coder models
void Compressor2Stage::delete_rc_models(void)
{
	// Flags
	for(int i = 0; i < RC_NO_FLAGS_MODEL; ++i)
		deleteqsmodel(&qsm_flags[i]);
    
	for(int i = 0; i < RC_NO_SYMBOLS_MODEL; ++i)
		deleteqsmodel(&qsm_symbols[i]);
    
	deleteqsmodel(&qsm_len_m1_pref);
	for(int i = 0; i < 6; ++i)
		deleteqsmodel(&qsm_len_m1_suf[i]);
    
	deleteqsmodel(&qsm_pos_m1_pref);
	for(int i = 0; i < 5; ++i)
		deleteqsmodel(&qsm_pos_m1_suf[i]);
    
	deleteqsmodel(&qsm_len_m2_pref);
	for(int i = 0; i < 7; ++i)
		deleteqsmodel(&qsm_len_m2_suf[i]);
    
	deleteqsmodel(&qsm_seq_pos_pref);
	for(int i = 0; i < 8; ++i)
		deleteqsmodel(&qsm_seq_pos_suf[i]);
    
	deleteqsmodel(&qsm_seq_id_pref);
	for(int i = 0; i < n_seq_ids_suf; ++i)
		deleteqsmodel(&qsm_seq_id_suf[i]);
	delete[] qsm_seq_id_suf;
}

void Compressor2Stage::putDataToQueue(uchar* d, int d_size, string filename, uchar eol_type, char new_line_at_end_of_file, uint32 line_width, uint32 *seq_sizes, vector<string> &_seq_names)
{
	QElem elem;
	elem.data = d;
	elem.dataSize = d_size;
	elem.filename = filename;
    elem.eol_type = eol_type;
    elem.line_width = line_width;
    elem.seq_sizes = seq_sizes;
    elem.seq_names.swap(_seq_names);
    elem.new_line_at_end_of_file = new_line_at_end_of_file;
	elementQueue->Push(elem);
}

void Compressor2Stage::start_decompression()
{
#ifdef COMPRESS_DESC
	gzFile_s *in_desc;
#else
	FILE *in_desc;
#endif

	// Read archive description
#ifdef WIN32
#ifdef COMPRESS_DESC
	in_desc = gzopen((archive_name + ".gdc2_desc").c_str(), "rb");
#else
	fopen_s(&in_desc, (archive_name + ".gdc2_desc").c_str(), "rb");
#endif
#else
#ifdef COMPRESS_DESC
	in_desc = gzopen((archive_name + ".gdc2_desc").c_str(), "rb");
#else
	in_desc = fopen((archive_name + ".gdc2_desc").c_str(), "rb");
#endif
#endif
    if(!in_desc)
    {
        cout << "Cannot open " << archive_name << ".gdc2_desc" << endl;
        exit(1);
    }
    string s, name;
    
	int i = 0;
	int c;
	int size = 0;
    uchar eol_type;
    char new_line_at_end_of_file;
    uint32 line_width;
    
    my_fread(&compression_level, 1, 1, in_desc);
    
    
	while ((c = my_getc(in_desc)) != EOF)
	{
		s.push_back(c);
        i++;
		if (c == 0 && i > 0)
		{
			my_fread(&size, 1, 4, in_desc);
			my_fread(&eol_type, 1, 1, in_desc);  //ad
            my_fread(&new_line_at_end_of_file, 1, 1, in_desc);  //ad
            my_fread(&line_width, 1, 4, in_desc);     //ad
            uint32 * _seq_sizes = new uint32[no_seq];
            vector<string> _seq_names;
            for(uint32 i = 0; i < no_seq; i++)
            {
                do
				{
					c = my_getc(in_desc);
                    name.push_back(c);
                }while(c != 0);
                _seq_names.push_back((name.c_str()));
                name.clear();
				my_fread(&_seq_sizes[i], 1, 4, in_desc);
            }
			my_getc(in_desc); //get end indicator (1)
            
            seq_names.push_back(_seq_names);
            seq_sizes.push_back(_seq_sizes);
			output_names.push_back(make_pair((s.c_str()), size));
            eol_types.push_back(eol_type);
            new_line_at_end_of_files.push_back(new_line_at_end_of_file);
            line_widths.push_back(line_width);
			s.clear();
            i = 0;
		}
	}
    
	my_fclose(in_desc);
	if (output_names.empty())
		return;
    
#ifdef WIN32
	fopen_s(&*rc_file, (archive_name + ".gdc2_rc").c_str(), "rb");
#else
	*rc_file = fopen((archive_name + ".gdc2_rc").c_str(), "rb");
#endif
	if (!*rc_file)
	{
		cout << "No archive file\n";
		return;
	}
	input_size = output_names.size();
    
    init_rc_models(0);
	start_decoding(&rc);
    
    //if all files should be decompressed
    if(inputFilenameVector->size()  == 0 || inputFilenameVector->size() == output_names.size() + 1)
    {
        // Decompress 1st sequence
#ifdef AD     
        cout << "*1*\n";
#endif
        
        int data_size = output_names[0].second;
        uchar *data = new uchar[data_size];
        unstore_1st_seq(data, data_size, true);

        putDataToQueue(data, data_size, output_names[0].first, eol_types[0], new_line_at_end_of_files[0], line_widths[0], seq_sizes.at(0), seq_names.at(0));
        file_data.push_back(make_pair(data, data_size));
        // Decompress next sequences
        cur_file_size = 0;
        cur_file_pos = 0;
        
        for (cur_file_id = 1; cur_file_id < output_names.size(); ++cur_file_id)
        {
#ifdef AD
            cout << "*" <<  cur_file_id+1 << "*\n";
#endif
            unprocess_data(data, data_size, true);
            putDataToQueue(data, data_size, output_names[cur_file_id].first, eol_types[cur_file_id], new_line_at_end_of_files[cur_file_id], line_widths[cur_file_id], seq_sizes.at(cur_file_id), seq_names.at(cur_file_id));
            file_data.push_back(make_pair(data, data_size));
        }
    }
    else    // only listed files should be decompressed
    {
        // determine files to decompress
        vector<uint32> seqToDecompress;
        string seqToDec;
        uint32 i, j;
        for(i = 0; i < inputFilenameVector->size(); i++)
        {
            seqToDec = inputFilenameVector->at(i);
            for(j = 0; j < output_names.size(); j++)
            {
                if(seqToDec.compare(output_names.at(j).first) == 0)
                {
                    seqToDecompress.push_back(j);
                    break;
                }
            }
            if(j == output_names.size() && (int32)i != ref_id_in_list)
                cout << "WARNING!!! Provided file name \"" << seqToDec << "\" is not present in the archive " << archive_name <<  "!!! [filename will be ommited]"<<endl;
        }
        seqToDecompress.push_back(output_names.size()+1);
        std::sort(seqToDecompress.begin(), seqToDecompress.end());
     
        uint32 id = 0;
        uint32 nextFile = seqToDecompress.at(id++);
        compression_level_threshold = (int)output_names.size()*compression_level*0.1;
        
        if(nextFile < output_names.size())
        {
        // Decompress 1st sequence
#ifdef AD
            cout << "*1*\n";
#endif
        
            int data_size = output_names[0].second;
            uchar *data = new uchar[data_size];
            
            if(0 < compression_level_threshold || nextFile == 0)
            {
                unstore_1st_seq(data, data_size, true);
                
                file_data.push_back(make_pair(data, data_size));
                if(nextFile == 0)
                {
                    putDataToQueue(data, data_size, output_names[0].first, eol_types[0], new_line_at_end_of_files[0], line_widths[0], seq_sizes.at(0), seq_names.at(0));
                
                    nextFile = seqToDecompress.at(id++);
                }
            }
            else
            {
                unstore_1st_seq(data, data_size, false);
            }
        
            // Decompress next sequences
            cur_file_size = 0;
            cur_file_pos = 0;
            for (cur_file_id = 1; cur_file_id < output_names.size() && nextFile < output_names.size(); ++cur_file_id)
            {
                if(cur_file_id < compression_level_threshold || nextFile == cur_file_id)
                {
    #ifdef AD
                    cout << "*" <<  cur_file_id+1 << "*\n";
    #endif
                    unprocess_data(data, data_size, true);
                    file_data.push_back(make_pair(data, data_size));
                
                    if(nextFile == cur_file_id)
                    {
                        putDataToQueue(data, data_size, output_names[cur_file_id].first, eol_types[cur_file_id], new_line_at_end_of_files[cur_file_id], line_widths[cur_file_id], seq_sizes.at(cur_file_id), seq_names.at(cur_file_id));
                        
                        nextFile = seqToDecompress.at(id++);
                    }
                }
                else
                {
                    unprocess_data(data, data_size, false);
                }
            }
        }
    }
    
	done_decoding(&rc);
	delete_rc_models();
    fclose(*rc_file);
	elementQueue->MarkCompleted();
}

void Compressor2Stage::unprocess_data(unsigned char*&data, int &data_size, bool ht)
{
	cur_file_size = output_names[cur_file_id].second;
    
	data_size = cur_file_size;
	if(ht)
        data = new unsigned char[data_size];
    
	v_last_seq_pos.clear();
	v_last_seq_pos.resize(cur_file_id);
    
	// Decompress sequence
	int encoded_symbols = 0;
	int encoded_tuples = 0;
	cur_file_pos = 0;
	expected_match_pos = 0;
	ctx_flag = 0;
    
	while (cur_file_pos < cur_file_size)
	{
		int syfreq, ltfreq;
		int flag;
        
		ltfreq = decode_culshift(&rc, RC_SHIFT_FLAGS);
		flag = qsgetsym(&qsm_flags[ctx_flag], ltfreq);
		qsgetfreq(&qsm_flags[ctx_flag], flag, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_FLAGS);
		qsupdate(&qsm_flags[ctx_flag], flag);
		ctx_flag = UPD_FLAG_CTX(ctx_flag, flag);
        
		if (flag == FLAG_MATCH)			// match
			unstore_match(data, encoded_symbols, encoded_tuples, ht);
		else if (flag == FLAG_LITERAL_MATCH)	// literal_1m
			unstore_literal_1m(data, encoded_symbols, encoded_tuples, ht);
		else						// literal_1l
			unstore_literal_1l(data, encoded_symbols, encoded_tuples, ht);
	}
    
	// Add mapping tuple_pos onto seq_pos
	if(ht)
        add_tuple_seq_pos_mapping(data, data_size, encoded_tuples);
}

void Compressor2Stage::unstore_match(unsigned char *data, int &encoded_symbols, int &encoded_tuples, bool ht)
{
	int syfreq, ltfreq;
    
	int seq_id;
	int seq_pos;
	int match_len;
    
	int si_pref, si_suf;
    
	ltfreq = decode_culshift(&rc, RC_SHIFT_SEQ_ID);
	si_pref = qsgetsym(&qsm_seq_id_pref, ltfreq);
	qsgetfreq(&qsm_seq_id_pref, si_pref, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SEQ_ID);
	qsupdate(&qsm_seq_id_pref, si_pref);
    
	ltfreq = decode_culshift(&rc, RC_SHIFT_SEQ_ID);
	si_suf = qsgetsym(&qsm_seq_id_suf[si_pref], ltfreq);
	qsgetfreq(&qsm_seq_id_suf[si_pref], si_suf, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SEQ_ID);
	qsupdate(&qsm_seq_id_suf[si_pref], si_suf);
    
	seq_id = (si_pref - 1) * 256 + si_suf;
    
	unsigned char *ref_data = file_data[seq_id].first;
	unsigned int tmp_seq_pos = 0;
    
	// Seq. pos.
	int sp_pref;
    
	ltfreq = decode_culshift(&rc, RC_SHIFT_SEQ_POS);
	sp_pref = qsgetsym(&qsm_seq_pos_pref, ltfreq);
	qsgetfreq(&qsm_seq_pos_pref, sp_pref, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SEQ_POS);
	qsupdate(&qsm_seq_pos_pref, sp_pref);
    
	// Len suf
	int sp_add = 0;
	if (sp_pref == 0)
		sp_add = 0;
	else if (sp_pref == 1)
		sp_add = 1;
	else if (sp_pref == 2)
		sp_add = 1;
	else if (sp_pref == 3)
		sp_add = 16;
	else if (sp_pref == 4)
		sp_add = 16;
    
	if (sp_pref == 0)
		tmp_seq_pos = 0;
	else if (sp_pref < 5)
	{
		int b;
		ltfreq = decode_culshift(&rc, RC_SHIFT_SEQ_POS);
		b = qsgetsym(&qsm_seq_pos_suf[sp_pref - 1], ltfreq);
		qsgetfreq(&qsm_seq_pos_suf[sp_pref - 1], b, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SEQ_POS);
		qsupdate(&qsm_seq_pos_suf[sp_pref - 1], b);
        
		tmp_seq_pos = b + sp_add;
	}
	else
        for (int i = 0; i < 4; ++i)
        {
            int b;
            
            ltfreq = decode_culshift(&rc, RC_SHIFT_SEQ_POS);
            b = qsgetsym(&qsm_seq_pos_suf[4 + i], ltfreq);
            qsgetfreq(&qsm_seq_pos_suf[4 + i], b, &syfreq, &ltfreq);
            decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SEQ_POS);
            qsupdate(&qsm_seq_pos_suf[4 + i], b);
            
            tmp_seq_pos += b << (i * 8);
        }
    
	int rel_seq_pos = (int)tmp_seq_pos;
	if (sp_pref > 0 && sp_pref % 2 == 0)
		rel_seq_pos *= -1;
    
	pair<int, int> expected_seq_pos;
	advance_in_sequence(seq_id, get<0>(v_last_seq_pos[seq_id]), get<1>(v_last_seq_pos[seq_id]), encoded_symbols - get<2>(v_last_seq_pos[seq_id]), expected_seq_pos);
    
	int tuple_pos = expected_seq_pos.first + rel_seq_pos;		// position in sequence in tuples;
	seq_pos = v_tuple_pos_map[seq_id][tuple_pos];
    
	get<0>(v_last_seq_pos[seq_id]) = tuple_pos;
	get<1>(v_last_seq_pos[seq_id]) = seq_pos;
	get<2>(v_last_seq_pos[seq_id]) = encoded_symbols;
    
	// Len pref
	int ml_pref;
    
	ltfreq = decode_culshift(&rc, RC_SHIFT_2M_LEN);
	ml_pref = qsgetsym(&qsm_len_m2_pref, ltfreq);
	qsgetfreq(&qsm_len_m2_pref, ml_pref, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_2M_LEN);
	qsupdate(&qsm_len_m2_pref, ml_pref);
    
	// Len suf
	int ml_inc;
	if (ml_pref == 0)
		ml_inc = 0;
	else if (ml_pref == 1)
		ml_inc = len_2m_model_sizes[0];
	else if (ml_pref == 2)
		ml_inc = len_2m_model_sizes[0] + len_2m_model_sizes[1];
	else if (ml_pref == 3)
		ml_inc = len_2m_model_sizes[0] + len_2m_model_sizes[1] + len_2m_model_sizes[2];
	else
		ml_inc = len_2m_model_sizes[0] + len_2m_model_sizes[1] + len_2m_model_sizes[2] + len_2m_model_sizes[3];
    
	match_len = 0;
	if (ml_pref < 4)
	{
		int b;
		ltfreq = decode_culshift(&rc, RC_SHIFT_2M_LEN);
		b = qsgetsym(&qsm_len_m2_suf[ml_pref], ltfreq);
		qsgetfreq(&qsm_len_m2_suf[ml_pref], b, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_2M_LEN);
		qsupdate(&qsm_len_m2_suf[ml_pref], b);
        
		match_len = b;
	}
	else
        for (int i = 0; i < 3; ++i)
        {
            int b;
            
            ltfreq = decode_culshift(&rc, RC_SHIFT_2M_LEN);
            b = qsgetsym(&qsm_len_m2_suf[ml_pref + i], ltfreq);
            qsgetfreq(&qsm_len_m2_suf[ml_pref + i], b, &syfreq, &ltfreq);
            decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_2M_LEN);
            qsupdate(&qsm_len_m2_suf[ml_pref + i], b);
            
            match_len += b << (i * 8);
        }
    
	match_len += ml_inc;
	unsigned char *ref_seq = file_data[seq_id].first;
	match_len = increase_match_len(ref_seq + seq_pos, match_len);
	determine_expected_pos(&ref_data[seq_pos], match_len);
	encoded_symbols += symbols_in_match(file_data[seq_id].first + seq_pos, match_len);
	encoded_tuples += match_len;
    
	int byte_pos = 0;
	for (int i = 0; i < match_len; ++i)
	{
		if (ref_data[seq_pos] < 128)				// 1st level match
		{
			if(ht)
                copy(&ref_data[seq_pos], &ref_data[seq_pos] + ref_match_len, &data[cur_file_pos + byte_pos]);
			seq_pos += ref_match_len;
			byte_pos += ref_match_len;
		}
		else									// 1st level literal
		{
			if(ht)
                copy(&ref_data[seq_pos], &ref_data[seq_pos] + ref_literal_len, &data[cur_file_pos + byte_pos]);
			seq_pos += ref_literal_len;
			byte_pos += ref_literal_len;
		}
	}
    
	cur_file_pos += byte_pos;
}

// Checking how many tuples from s, but omitting the ones present in hashing_len_bytes
int Compressor2Stage::increase_match_len(unsigned char *s, int len)
{
	int pos = 0;
    
	while (pos < hashing_len_bytes)
	{
		if (s[pos] < 128)			// match
			pos += ref_match_len;
		else
			pos += ref_literal_len;
		++len;
	}
    
	return len;
}

void Compressor2Stage::unstore_1st_seq(unsigned char* data, int size, bool ht)
{
	cur_file_size = output_names[0].second;
	cur_file_pos = 0;
	ctx_flag = 0;
	int tmp;
	int encoded_tuples = 0;
    
	while (cur_file_pos < cur_file_size)
	{
		int syfreq, ltfreq;
		int flag;
        
		ltfreq = decode_culshift(&rc, RC_SHIFT_FLAGS);
		flag = qsgetsym(&qsm_flags[ctx_flag], ltfreq);
		qsgetfreq(&qsm_flags[ctx_flag], flag, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_FLAGS);
		qsupdate(&qsm_flags[ctx_flag], flag);
		ctx_flag = UPD_FLAG_CTX(ctx_flag, flag);
        
		if (flag == FLAG_LITERAL_MATCH)	// literal_1m
			unstore_literal_1m(data, tmp, tmp, ht);
		else						// literal_1l
			unstore_literal_1l(data, tmp, tmp, ht);
		++encoded_tuples;
	}
    
	if(ht)
        add_tuple_seq_pos_mapping(data, size, encoded_tuples);
}

void Compressor2Stage::unstore_literal_1m(unsigned char *data, int &encoded_symbols, int &encoded_tuples, bool ht)
{
	int syfreq, ltfreq;
	int pref;
	int ref_pos;
    
	// Len pref
	ltfreq = decode_culshift(&rc, RC_SHIFT_LEN);
	pref = qsgetsym(&qsm_len_m1_pref, ltfreq);
	qsgetfreq(&qsm_len_m1_pref, pref, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_LEN);
	qsupdate(&qsm_len_m1_pref, pref);
    
	// Len suf
	int ctx_suf = 0;
	int len = 0;
    
	if (pref == 0)
		ctx_suf = 0;
	else if (pref == 1)
	{
		ctx_suf = 1;
		len = 0x100;
	}
	else if (pref == 2)
	{
		ctx_suf = 3;
		len = 0x10000 + 0x100;
	}
    
	for (int i = 0; i <= pref; ++i)
	{
		int b;
        
		ltfreq = decode_culshift(&rc, RC_SHIFT_LEN);
		b = qsgetsym(&qsm_len_m1_suf[ctx_suf + i], ltfreq);
		qsgetfreq(&qsm_len_m1_suf[ctx_suf + i], b, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_LEN);
		qsupdate(&qsm_len_m1_suf[ctx_suf + i], b);
		len += b << (8 * i);
	}
    
    if(ht)
    {
        data[cur_file_pos + 4] = (len >> 16);
        data[cur_file_pos + 5] = (len >> 8) & 0xff;
        data[cur_file_pos + 6] = (len)& 0xff;
    }
    
	int pos_to_code;
	ltfreq = decode_culshift(&rc, RC_SHIFT_1M_POS);
	pref = qsgetsym(&qsm_pos_m1_pref, ltfreq);
	qsgetfreq(&qsm_pos_m1_pref, pref, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_1M_POS);
	qsupdate(&qsm_pos_m1_pref, pref);
    
	if (pref == 0)
		pos_to_code = 0;
	else if (pref == 1)
	{
		int b;
		ltfreq = decode_culshift(&rc, RC_SHIFT_1M_POS);
		b = qsgetsym(&qsm_pos_m1_suf[0], ltfreq);
		qsgetfreq(&qsm_pos_m1_suf[0], b, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_1M_POS);
		qsupdate(&qsm_pos_m1_suf[0], b);
        
		if (b >= 2 * pos_short_limits - 2)
		{
			if (b == 2 * pos_short_limits - 2)
				pos_to_code = PSEUDO_MATCH_N;
			else
				pos_to_code = PSEUDO_MATCH_n;
		}
		else if (b < pos_short_limits - 1)
			pos_to_code = b + 1;
		else
			pos_to_code = b - (2 * pos_short_limits - 2);
	}
	else
	{
		unsigned int p = 0;
		for (int i = 0; i < 4; ++i)
		{
			int b;
			ltfreq = decode_culshift(&rc, RC_SHIFT_1M_POS);
			b = qsgetsym(&qsm_pos_m1_suf[i + 1], ltfreq);
			qsgetfreq(&qsm_pos_m1_suf[i + 1], b, &syfreq, &ltfreq);
			decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_1M_POS);
			qsupdate(&qsm_pos_m1_suf[i + 1], b);
            
			p += b << (8 * i);
		}
		pos_to_code = (int)p;
	}
    
	if (pos_to_code != PSEUDO_MATCH_N && pos_to_code != PSEUDO_MATCH_n)
	{
		ref_pos = expected_match_pos - pos_to_code;
		expected_match_pos = ref_pos + len;
	}
	else
	{
		ref_pos = pos_to_code;
		expected_match_pos += len;
	}
    
	if(ht)
    {
        data[cur_file_pos + 0] = (ref_pos >> 24) & 0xff;
        data[cur_file_pos + 1] = (ref_pos >> 16) & 0xff;
        data[cur_file_pos + 2] = (ref_pos >> 8)  & 0xff;
        data[cur_file_pos + 3] = (ref_pos)       & 0xff;
    }
    
	encoded_symbols += len;
	++encoded_tuples;
    
	cur_file_pos += ref_match_len;
}

void Compressor2Stage::unstore_literal_1l(unsigned char *data, int &encoded_symbols, int &encoded_tuples, bool ht)
{
	int syfreq, ltfreq;
	int symbol;
    
	ltfreq = decode_culshift(&rc, RC_SHIFT_SYMBOLS);
	symbol = qsgetsym(&qsm_symbols[ctx_symbol], ltfreq);
	qsgetfreq(&qsm_symbols[ctx_symbol], symbol, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SYMBOLS);
	qsupdate(&qsm_symbols[ctx_symbol], symbol);
    
	if(ht)
        data[cur_file_pos] = (unsigned char)(symbol + 0x80 + 0x40);
    
	ctx_symbol = ctx_from_symbol(symbol + 0x40);
    
	++expected_match_pos;
	++encoded_symbols;
	++encoded_tuples;
    
	cur_file_pos += ref_literal_len;
}

// Adding mapping tuple_pos -> seq_pos
void Compressor2Stage::add_tuple_seq_pos_mapping(unsigned char *p, int size, int encoded_tuples)
{
	v_tuple_pos_map.push_back(vector<int>());
	vector<int> &mapping = v_tuple_pos_map.back();
	mapping.reserve(encoded_tuples);
    
	for (int i = 0; i < size;)
	{
		mapping.push_back(i);
		if (p[i] < 128)			// 1st level match
			i += ref_match_len;
		else
			i += ref_literal_len;
	}
}