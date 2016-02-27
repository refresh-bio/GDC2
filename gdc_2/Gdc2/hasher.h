/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef _HASHER_H
#define _HASHER_H
#include "defs.h"

#include <vector>

using namespace std;

class CHasher
{
public:
	vector<uchar*> seqs_ref;
	vector<uint64> seqs_size;
	vector<uint64> seqs_alloc;
	int32 refs_no;

	vector<int*> hts_ref;
	vector<int64> hts_size;
	vector<int64> hts_elems;
	vector<uint64> hash_vals_ref;

	int *cur_hts_ref;
	int64 cur_hts_size;
	int64 cur_hts_elems;
	uint64 *cur_hash_vals_ref;

	uint64 hash_value;
	int32 hashing_step;
	int32 hashing_colisions;
	int32 cur_hashing_step;
	int32 cur_hashing_seq;
	int32 min_match_ext;
	int32 min_match_len;

	static double ht_load_factor;
	int32 hash_key_len;

	inline uint64 Hash1(uchar *p, uint64 ht_size, uint32 &num_N);
	inline uint64 Hash1Inc(uchar *p, uint64 ht_size, uint32 &num_N);
	inline uint64 Hash1Inc1(uchar *p, uint64 ht_size, uint32 &num_N);

	bool MakeRefHashing(int col_id);

	/* Concurrent */
	inline uint64 ConcurrentHash1(uchar *p, uint64 ht_size, uint32 &num_N, uint64 &hash_value);
	inline uint64 ConcurrentHash1Inc1(uchar *p, uint64 ht_size, uint32 &num_N, uint64 &hash_value);
	/* End concurrent */

public:
	CHasher(int32 _hashing_step = 1, int32 _hash_key_len = 15);
	~CHasher();
	void SetParams(int32 _hashing_step = 1, int32 _hashing_colisions = 0, int32 _hash_key_len = 15, int32 _min_match_ext = 0);

	bool SetRefSeq(uchar *_seq_ref, uint64 _size);
	inline bool FindFirst(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, bool full_comp);
	inline bool FindNext(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id);

	/* Concurrent */
	inline bool ConcurrentFindFirst(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, bool full_comp,
										 int32 &cur_hashing_step, int32 &cur_hashing_seq, uint64* &cur_hash_vals_ref, int* &cur_hts_ref,
										 int64 &cur_hts_size, int64 &cur_hts_elems, uint64 &hash_value, vector<uint64> &hash_vals_ref);
	inline bool ConcurrentFindNext(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, int* &cur_hts_ref,
										uint64* &cur_hash_vals_ref, int32 &cur_hashing_step, int32 &cur_hashing_seq, int64 &cur_hts_size, 
										int64 &cur_hts_elems, uint64 &hash_value, vector<uint64> &hash_vals_ref);
	/* End concurrent */

};

// ********************************************************************
inline uint64 CHasher::Hash1(uchar *p, uint64 ht_size, uint32 &num_N)
{
	hash_value = 0;
	num_N = 0;

	for(int32 i = 0; i < hash_key_len; ++i, ++p)
	{
		hash_value <<= 2;
		hash_value ^= (*p << 9) - *p;
		num_N += (*p == 'N') || (*p == 'n');
	}

	return ((hash_value << 7) + hash_value) % ht_size;
}

inline uint64 CHasher::ConcurrentHash1(uchar *p, uint64 ht_size, uint32 &num_N, uint64 &hash_value)
{
	hash_value = 0;
	num_N = 0;

	for(int32 i = 0; i < hash_key_len; ++i, ++p)
	{
		hash_value <<= 2;
		hash_value ^= (*p << 9) - *p;
		num_N += (*p == 'N') || (*p == 'n');
	}

	return ((hash_value << 7) + hash_value) % ht_size;
}

// ********************************************************************
inline uint64 CHasher::Hash1Inc(uchar *p, uint64 ht_size, uint32 &num_N)
{
	for(int32 i = 0; i < hashing_step; ++i)
	{
		hash_value ^= (p[-(hashing_step-i)] * (uint64) 511) << (uint64) (hash_key_len*2-2);
		num_N -= p[-(hashing_step-i)] == 'N' || p[-(hashing_step-i)] == 'n';
		hash_value <<= 2;
		hash_value ^= (p[hash_key_len - hashing_step+i] << 9) - p[hash_key_len - hashing_step+i];
		num_N += p[hash_key_len - hashing_step+i] == 'N' || p[hash_key_len - hashing_step + i] == 'n';
	}
	
	return ((hash_value << 7) + hash_value) % ht_size;
}

// ********************************************************************
inline uint64 CHasher::Hash1Inc1(uchar *p, uint64 ht_size, uint32 &num_N)
{
	hash_value ^= (p[-1] * (uint64) 511) << (uint64) (hash_key_len*2-2);
	num_N -= p[-1] == 'N' || p[-1] == 'n';
	hash_value <<= 2;
	hash_value ^= p[hash_key_len-1] * (uint64) 511;
	num_N += p[hash_key_len-1] == 'N' || p[hash_key_len-1] == 'n';
	
	return ((hash_value << 7) + hash_value) % ht_size;
}

inline uint64 CHasher::ConcurrentHash1Inc1(uchar *p, uint64 ht_size, uint32 &num_N, uint64 &hash_value)
{
	hash_value ^= (p[-1] * (uint64) 511) << (uint64) (hash_key_len*2-2);
	num_N -= p[-1] == 'N' || p[-1] == 'n';
	hash_value <<= 2;
	hash_value ^= p[hash_key_len-1] * (uint64) 511;
	num_N += p[hash_key_len-1] == 'N' || p[hash_key_len-1] == 'n';
	
	return ((hash_value << 7) + hash_value) % ht_size;
}

// ********************************************************************
inline bool CHasher::FindFirst(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, bool full_comp)
{
	uint32 tmp;

	cur_hashing_step = 0;
	cur_hashing_seq  = 0;		// aug. seq.

	cur_hash_vals_ref = &hash_vals_ref[0];
	*cur_hash_vals_ref = Hash1(p, hts_size[0], tmp);

	cur_hts_ref   = hts_ref[0];
	cur_hts_size  = hts_size[0];
	cur_hts_elems = hts_elems[0];

	return FindNext(p, max_errors, pos, lens, col_id);
}

inline bool CHasher::ConcurrentFindFirst(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, bool full_comp,
										 int32 &cur_hashing_step, int32 &cur_hashing_seq, uint64* &cur_hash_vals_ref, int* &cur_hts_ref,
										 int64 &cur_hts_size, int64 &cur_hts_elems, uint64 &hash_value, vector<uint64> &hash_vals_ref)
{
	uint32 tmp;

	cur_hashing_step = 0;
	cur_hashing_seq  = 0;		// aug. seq.

	cur_hash_vals_ref = &hash_vals_ref[0];
	*cur_hash_vals_ref = ConcurrentHash1(p, hts_size[0], tmp, hash_value);

	cur_hts_ref   = hts_ref[0];
	cur_hts_size  = hts_size[0];
	cur_hts_elems = hts_elems[0];

	return ConcurrentFindNext(p, max_errors, pos, lens, col_id, cur_hts_ref, cur_hash_vals_ref, cur_hashing_step, cur_hashing_seq, cur_hts_size,
		cur_hts_elems, hash_value, hash_vals_ref);
}

// ********************************************************************
inline bool CHasher::FindNext(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id)
{
	pos = cur_hts_ref[*cur_hash_vals_ref];	lens.clear();

	while(!pos)
	{
		uint32 tmp;
		++cur_hashing_step;
		if(cur_hashing_step >= hashing_step)
		{
			if(cur_hashing_seq+1 < refs_no)			
			{
				cur_hashing_seq++;
				cur_hash_vals_ref = &hash_vals_ref[cur_hashing_seq];
				cur_hts_ref   = hts_ref[cur_hashing_seq];
				cur_hts_size  = hts_size[cur_hashing_seq];
				cur_hts_elems = hts_elems[cur_hashing_seq];

				cur_hashing_step = 0;
				*cur_hash_vals_ref = Hash1(p, cur_hts_size, tmp);
				pos = cur_hts_ref[*cur_hash_vals_ref];
				continue;
			}
			else
				return false;
		}
		*cur_hash_vals_ref = Hash1Inc1(p+cur_hashing_step, cur_hts_size, tmp);
		pos = cur_hts_ref[*cur_hash_vals_ref];
	}

	if(pos < (uint32) cur_hashing_step)
	{
		col_id = cur_hashing_seq;
		(*cur_hash_vals_ref)++;
		return true;
	}
	pos -= cur_hashing_step;

	int32 n_err = 0;
	int32 len;
	uchar *s;
	s = seqs_ref[cur_hashing_seq];
	int32 last_error_pos = -1;
	int32 max_mismatches = max_errors;

    
	for(len = 0; n_err <= max_mismatches; ++len)
	{
		if(p[len] == s[pos+len])
			continue;

		if(len-last_error_pos-1 < (int32) min_match_ext)
			break;
		if(n_err == 0 && len-last_error_pos-1 < (int32) min_match_len)
			break;
		lens.push_back(len-last_error_pos-1);
		last_error_pos = len;

		n_err++;
	}

	(*cur_hash_vals_ref)++;
	if((int64) *cur_hash_vals_ref >= cur_hts_size)
		*cur_hash_vals_ref = 0;

	col_id = cur_hashing_seq;
    
	return true;
}

inline bool CHasher::ConcurrentFindNext(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, int* &cur_hts_ref,
										uint64* &cur_hash_vals_ref, int32 &cur_hashing_step, int32 &cur_hashing_seq, int64 &cur_hts_size, 
										int64 &cur_hts_elems, uint64 &hash_value, vector<uint64> &hash_vals_ref)
{
	pos = cur_hts_ref[*cur_hash_vals_ref];	lens.clear();

	while(!pos)
	{
		uint32 tmp;
		++cur_hashing_step;
		if(cur_hashing_step >= hashing_step)
		{
			if(cur_hashing_seq+1 < refs_no)			
			{
				cur_hashing_seq++;
				cur_hash_vals_ref = &hash_vals_ref[cur_hashing_seq];
				cur_hts_ref   = hts_ref[cur_hashing_seq];
				cur_hts_size  = hts_size[cur_hashing_seq];
				cur_hts_elems = hts_elems[cur_hashing_seq];

				cur_hashing_step = 0;
				*cur_hash_vals_ref = ConcurrentHash1(p, cur_hts_size, tmp, hash_value);
				pos = cur_hts_ref[*cur_hash_vals_ref];
				continue;
			}
			else
				return false;
		}
		*cur_hash_vals_ref = ConcurrentHash1Inc1(p+cur_hashing_step, cur_hts_size, tmp, hash_value);
		pos = cur_hts_ref[*cur_hash_vals_ref];
	}

	if(pos < (uint32) cur_hashing_step)
	{
		col_id = cur_hashing_seq;
		(*cur_hash_vals_ref)++;
		return true;
	}
	pos -= cur_hashing_step;

	int32 n_err = 0;
	int32 len;
	uchar *s;
	s = seqs_ref[cur_hashing_seq];
	int32 last_error_pos = -1;
	int32 max_mismatches = max_errors;

    
	for(len = 0; n_err <= max_mismatches; ++len)
	{
		if(p[len] == s[pos+len])
			continue;

		if(len-last_error_pos-1 < (int32) min_match_ext)
			break;
		if(n_err == 0 && len-last_error_pos-1 < (int32) min_match_len)
			break;
		lens.push_back(len-last_error_pos-1);
		last_error_pos = len;

		n_err++;
	}

	(*cur_hash_vals_ref)++;
	if((int64) *cur_hash_vals_ref >= cur_hts_size)
		*cur_hash_vals_ref = 0;

	col_id = cur_hashing_seq;
    
	return true;
}

#endif
