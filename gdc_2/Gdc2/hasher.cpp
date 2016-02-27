/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#include "stdafx.h"


#include <algorithm>
#include <iostream>
#include "hasher.h"

using namespace std;

double CHasher::ht_load_factor = 0.65;

// ********************************************************************
CHasher::CHasher(int32 _hashing_step, int32 _hash_key_len) :  hashing_step(_hashing_step), hash_key_len(_hash_key_len)
{
	hashing_colisions = 1 << 30;
	refs_no = 0;
}

// ********************************************************************
CHasher::~CHasher()
{
	if(!seqs_ref.empty() && seqs_ref[0])
		delete[] seqs_ref[0];
	// Do not delete seq_ref, since this sequence is maintained elsewhere

	for(size_t i = 0; i < hts_ref.size(); ++i)
		if(hts_ref[i])
			delete[] hts_ref[i];
}

// ********************************************************************
void CHasher::SetParams(int32 _hashing_step, int32 _hashing_colisions, int32 _hash_key_len, int32 _min_match_ext)
{
	hashing_step      = _hashing_step;
	hashing_colisions = _hashing_colisions;
	hash_key_len	  = _hash_key_len;
	min_match_ext	  = _min_match_ext;
	min_match_len	  = _hash_key_len;

	if(hashing_colisions == 0)
		hashing_colisions = 1 << 30;
}

// ********************************************************************
// Important: 0th position of seq_ref is for special purposes, so must be unused
bool CHasher::SetRefSeq(uchar *_seq_ref, uint64 _size)
{
	if(!refs_no)
	{
		// Add empty aug. seq.
		seqs_alloc.push_back(1 << 20);
		seqs_size.push_back(1);
		seqs_ref.push_back(new uchar[seqs_alloc[0]]);

		hts_size.push_back((uint64) ((seqs_alloc[0]+hashing_step) / hashing_step / ht_load_factor));
		hts_elems.push_back(0);
		hts_ref.push_back(new int32[hts_size[0]]);
		fill_n(hts_ref[0], hts_size[0], 0);
		hash_vals_ref.push_back(0);

		refs_no = 1;
	}

	seqs_alloc.push_back(_size);
	seqs_size.push_back(_size);
	seqs_ref.push_back(_seq_ref);

	hts_size.push_back((uint64) ((seqs_alloc.back()+hashing_step) / hashing_step / ht_load_factor));
	hts_elems.push_back(0);
	hts_ref.push_back(new int32[hts_size.back()]);
	hash_vals_ref.push_back(0);
	fill_n(hts_ref.back(), hts_size.back(), 0);

	_seq_ref[0] = 255;

	return MakeRefHashing(refs_no++);
}

// ********************************************************************
bool CHasher::MakeRefHashing(int32 col_id)
{
	uint32 num_N = 0;
	uint64 h = Hash1(seqs_ref[col_id]+1, hts_size[col_id], num_N);
	hts_elems[col_id] = 0;
	uint32 next_rap_pos = 0;
	uint32 rap_step = (1 << 20) * hashing_step;

	for(uint64 i = 1; i + hash_key_len <= seqs_size[col_id];)
	{
		if(i >= next_rap_pos)
		{
			next_rap_pos += rap_step;
		}

		if(num_N > 3)
		{
			i += hashing_step;
			h = Hash1Inc(seqs_ref[col_id]+i, hts_size[col_id], num_N);
			continue;
		}

		uint64 h0 = h;
		for(; hts_ref[col_id][h]; h = ((int64) h < hts_size[col_id]-1) ? h+1 : 0)
			;

		if(h-h0 < (uint64)hashing_colisions)
		{
			hts_ref[col_id][h] = (int32) i;
			hts_elems[col_id]++;
		}

		i += hashing_step;
		h = Hash1Inc(seqs_ref[col_id]+i, hts_size[col_id], num_N);
	}

	return true;
}
