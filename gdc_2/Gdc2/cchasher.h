/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef CONCURRENT_CHASHER_H
#define CONCURRENT_CHASHER_H

#include <vector>
#include "defs.h"
#include "hasher.h"

using namespace std;

class ConcurrentCHasher
{
	CHasher* _hasher;

	int *cur_hts_ref;
	int64 cur_hts_size;
	int64 cur_hts_elems;
	uint64 *cur_hash_vals_ref;
	int32 cur_hashing_step;
	int32 cur_hashing_seq;
	uint64 hash_value;

	vector<uint64> hash_vals_ref;
public:
	ConcurrentCHasher(CHasher* hasher)
	{
		_hasher = hasher;
		
		hash_vals_ref.reserve(hasher->hash_vals_ref.size());
		copy(hasher->hash_vals_ref.begin(),hasher->hash_vals_ref.end(),back_inserter(hash_vals_ref));
	}

	inline bool FindFirst(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, bool full_comp);
	inline bool FindNext(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id);
};

inline bool ConcurrentCHasher::FindFirst(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id, bool full_comp)
{
	return _hasher->ConcurrentFindFirst(p, max_errors, pos, lens, col_id, full_comp, cur_hashing_step, cur_hashing_seq, 
		cur_hash_vals_ref, cur_hts_ref, cur_hts_size, cur_hts_elems, hash_value, hash_vals_ref);
}

inline bool ConcurrentCHasher::FindNext(uchar *p, int32 max_errors, uint32 &pos, vector<int32> &lens, int32 &col_id)
{
	return _hasher->ConcurrentFindNext(p, max_errors, pos, lens, col_id, cur_hts_ref, cur_hash_vals_ref, cur_hashing_step, 
		cur_hashing_seq, cur_hts_size, cur_hts_elems, hash_value, hash_vals_ref);
}

#endif