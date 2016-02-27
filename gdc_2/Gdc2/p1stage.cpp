/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#include "p1stage.h"

#define MAX(a, b) (a >= b ? a : b)
#define POS_MAX ((unsigned int)(1<<31)-3)
#define LEN_MAX ((1<<24)-1)

void Producer1Stage::start()
{
    string filename;
    uint32 max_len;
	while (!file_input_queue->IsCompleted())
	{
		if (file_input_queue->Pop(filename))
		{
            if (readSeq(filename))
			{
				uchar* outData = NULL;
				int outDataSize;
			
                addCollectionDif(filename, data, outData, outDataSize);
                if(data.size() != no_seq)
                {
                    cout << "Number of sequences in " << filename << " does not match number of sequences in reference! (" << data.size() << " != " << no_seq << ")" << endl;
                    exit(1);
                }
                
                uint32 * seq_sizes = new uint32[no_seq];
                vector<string> seq_names;
                seq_names.clear();
                
                for(uint32 i = 0; i < no_seq; i++){
                    seq_sizes[i] = data[i].size;
                    seq_names.push_back(string(data[i].seq_name));
                }
                
                max_len = 0;
                
                for(size_t i = 0; i < data.size(); i++)
                {
                    max_len = max_len > data[i].line_len ? max_len : data[i].line_len;
                }

                
                addOutputDataDoQueue(outData, outDataSize, filename, data[0].eol_type, max_len, seq_sizes, seq_names, fasta->new_line_at_end_of_file);
                
              
                
                fasta->Close();
			}
            else
            {
                cout << "Cannot open file from the list: " << filename << endl;
                exit(1);
            }
		}
	}
	_elementQueue->MarkCompleted();
}

void Producer1Stage::addOutputDataDoQueue(uchar* outData, int outDataSize, string filename, uchar eol_type, uint32 line_width, uint32 *_seq_sizes, vector<string> &_seq_names, char _new_line_at_end_of_file)
{
	QElem e;
	e.filename = filename;
	e.data = outData;
	e.dataSize = outDataSize;
    e.eol_type = eol_type;
    e.seq_sizes = _seq_sizes;
    e.line_width = line_width;
    e.seq_names.swap(_seq_names);
    e.new_line_at_end_of_file = _new_line_at_end_of_file;
	_elementQueue->Push(e);
}

bool Producer1Stage::readSeq(string filename)
{
	if(!fasta->Open(filename))
	{
		return false;
	}
	fasta->Read(data);
	return true;
}

bool Producer1Stage::addCollectionDif(string &col_name, vector<CFastaSequence> &data, uchar* &outData, int &outDataSize)
{
	vector<uchar> rv; // result vector
	uint32 tot_len = 0;
	for(size_t i = 0; i < data.size(); ++i)
		tot_len += data[i].size;
    
    uchar *seq_alt;
    seq_alt = new uchar[tot_len+14];
	A_memcpy(seq_alt+1, data[0].raw_data, tot_len);
    
	seq_alt[tot_len+2] = '\0';
	fill_n(seq_alt+tot_len+3, 10, 1);

    for(uint32 i = 1  ; i <= tot_len;)
	{
        uint32 tmp;
        uint32 pos;
        int32 col_id;
        vector<int32> lens;
        vector<int32> best_lens;
        uint32 literals_start = 0;
		uint32 best_pos = 0, best_len = 0;
		uint32 best_len1 = 0;
		best_lens.clear();
        
        
        bool go;
        uint32 next_pos, ref_pos;
        uint32 after_snp_len, after_ins_len, after_del_len, after_ins2_len, after_del2_len, max_len;
        
        
		int32 N_run_len = NRunLen(seq_alt, i);
        uint32 tmplen;
		if(N_run_len)
		{
			best_pos = 0;
			best_len = best_len1 = N_run_len;
		}
        
        if(N_run_len)
		{
			if(N_run_len > 0)
			{
				tmp = POS_MAX + (uint32)2;
                writeBE(tmp, rv);
                tmplen = N_run_len;
                while(writeBElen(tmplen, rv))
                {
                    tmplen = tmplen - LEN_MAX;
                    writeBE(tmp, rv);
                }
                i += (uint32) N_run_len;
			}
			else
			{
                tmp = POS_MAX + (uint32)1;
                writeBE(tmp, rv);
                
                N_run_len = - N_run_len;
                tmplen = N_run_len;
                while(writeBElen(tmplen, rv))
                {
                    tmplen = tmplen - LEN_MAX;
                    writeBE(tmp, rv);
                }
                i += (uint32) N_run_len;
			}
		}
        else if(hasher->FindFirst(seq_alt+i, 0, pos, lens, col_id, literals_start == 0))
		{
			do
			{
				uint32 sum_len = accumulate(lens.begin(), lens.end(), 0);
                
                if(sum_len >= min_match_1st_len && sum_len > best_len) /// total_len changed to min_match_1st_len!
                {
					best_len  = sum_len;
					best_pos  = pos;
					best_lens.swap(lens);
				}
			} while(hasher->FindNext(seq_alt+i, 0, pos, lens, col_id));
            
			if(best_lens.empty())
			{
                tmp = seq_alt[i]|128;
				rv.push_back((uchar)(tmp & 0x000000ff)); // do zmiany
                ++i;
			}
			else
			{
                tmp = 0;
                if(best_pos > POS_MAX)
                {
                    cout << "Position out of range: " << best_pos << endl;
                    exit(2);
                }
                writeBE(best_pos, rv);
                
                tmp = best_pos;
                tmplen = best_len;
                while(writeBElen(tmplen, rv))
                {
                    tmplen = tmplen - LEN_MAX;
                    tmp = tmp + LEN_MAX;
                    writeBE(tmp, rv);
                }
                
				i += (uint32) (accumulate(best_lens.begin(), best_lens.end(), 0));
                
                next_pos = i;
                ref_pos = best_pos+best_len;
                
                do
                {
                    go = false;
                    after_snp_len = getMatchLen(seq_ref+ref_pos+1, seq_alt+next_pos+1);
                    after_ins_len = getMatchLen(seq_ref+ref_pos, seq_alt+next_pos+1);
                    after_del_len = getMatchLen(seq_ref+ref_pos+1, seq_alt+next_pos);

                    if (indel2)
                    {
                        after_ins2_len = getMatchLen(seq_ref+ref_pos, seq_alt+next_pos+2);
                        after_del2_len = getMatchLen(seq_ref+ref_pos+2, seq_alt+next_pos);
                    }
                    else
                    {
                        after_ins2_len = after_del2_len = 0;
                    }
                    
                    max_len = MAX(MAX(MAX(MAX(after_snp_len, after_ins_len), after_del_len), after_ins2_len), after_del2_len);
                    
                    if (max_len >= min_after_match)
                    {
                        if (max_len == after_snp_len)
                        {
                            tmp = seq_alt[next_pos]|128;
                            
							rv.push_back((uchar)(tmp & 0x000000ff));
                            tmp = ref_pos+1;
                            if(tmp > POS_MAX)
                            {
                                cout << "Position out of range: " << best_pos << endl;
                                exit(2);
                            }
                            writeBE(tmp, rv);
                            
                            tmplen = after_snp_len;
                            while(writeBElen(tmplen, rv))
                            {
                                tmplen = tmplen - LEN_MAX;
                                tmp = tmp + LEN_MAX;
                                writeBE(tmp, rv);
                            }
                            
                            ref_pos = ref_pos + 1 + after_snp_len ;
                            next_pos = next_pos + 1 + after_snp_len ;
                            go = true;
                            
                        }
                        else if (max_len == after_ins_len)
                        {
    
                            tmp = seq_alt[next_pos]|128;
                            
							rv.push_back((uchar)(tmp & 0x000000ff));
							
                            tmp = 0;
                            if(ref_pos > POS_MAX)
                            {
                                cout << "Position out of range: " << best_pos << endl;
                                exit(2);
                            }
                            writeBE(ref_pos, rv);
                            
                            tmp = ref_pos;
                            tmplen = after_ins_len;
                            while(writeBElen(tmplen, rv))
                            {
                                tmplen = tmplen - LEN_MAX;
                                tmp = tmp + LEN_MAX;
                                writeBE(tmp, rv);
                            }
                            
                            
                            ref_pos = ref_pos + after_ins_len ;
                            next_pos = next_pos + 1 + after_ins_len ;
                            go = true;
                        }
                        else if (max_len == after_del_len)
                        {
                            tmp = ref_pos+1;
                            if(tmp > POS_MAX)
                            {
                                cout << "Position out of range: " << best_pos << endl;
                                exit(2);
                            }
                            writeBE(tmp, rv);
                            
                            tmplen = after_del_len;
                            while(writeBElen(tmplen, rv))
                            {
                                tmplen = tmplen - LEN_MAX;
                                tmp = tmp + LEN_MAX;
                                writeBE(tmp, rv);
                            }
                            
                            ref_pos = ref_pos + 1 + after_del_len ;
                            next_pos = next_pos + after_del_len ;
                            go = true;
                            
                        }
                        else if (max_len == after_ins2_len)
                        {
                            tmp = seq_alt[next_pos]|128;
							rv.push_back((uchar)(tmp & 0x000000ff));
                            
                            tmp = seq_alt[next_pos+1]|128;
							rv.push_back((uchar)(tmp & 0x000000ff));
                            
                            tmp = 0;
                            if(ref_pos > POS_MAX)
                            {
                                cout << "Position out of range: " << best_pos << endl;
                                exit(2);
                            }
                            writeBE(ref_pos, rv);
                            tmp = ref_pos;
                            tmplen = after_ins2_len;
                            while(writeBElen(tmplen, rv))
                            {
                                tmplen = tmplen - LEN_MAX;
                                tmp = tmp + LEN_MAX;
                                writeBE(tmp, rv);
                            }
                            
                            
                            ref_pos = ref_pos + after_ins2_len ;
                            next_pos = next_pos + 2 + after_ins2_len ;
                            go = true;
                        }
                        else // (max_len == after_del2_len)
                        {
                            tmp = ref_pos+2;
                            if(tmp > POS_MAX)
                            {
                                cout << "Position out of range: " << best_pos << endl;
                                exit(2);
                            }
                            writeBE(tmp, rv);
                           
                            tmplen = after_del2_len;
                            while(writeBElen(tmplen, rv))
                            {
                                tmplen = tmplen - LEN_MAX;
                                tmp = tmp + LEN_MAX;
                                writeBE(tmp, rv);
                            }
                            ref_pos = ref_pos + 2 + after_del2_len ;
                            next_pos = next_pos + after_del2_len ;
                            go = true;
                        }
                    }
                } while (go);
                
                i = next_pos;
            }
		}
        else
		{
            tmp = seq_alt[i]|128;
			rv.push_back((uchar)(tmp & 0x000000ff));
            ++i;
		}
    }
    delete [] seq_alt;
    
	outData = new uchar[rv.size()];
	outDataSize = rv.size();
	A_memcpy(outData, &rv[0], outDataSize);
    
	return true;
}

int32 Producer1Stage::NRunLen(uchar* seq, uint32 pos)
{
	uint32 r = 0;
	uchar cN = seq[pos];
    
	if(cN != 'N' && cN != 'n')
		return 0;
    
	for(; seq[pos+r] == cN; ++r)
		;
    
	if(r < (uint32) min_match_1st_len)
		return 0;
	else if(cN == 'N')
		return (int32) r;
	else
		return - ((int32) r);
}

void Producer1Stage::writeBE(uint32 tmp, vector<uchar> &rv)
{
    uchar b1, b2, b3, b4;
    
    b1 = (tmp >> 24) & 0x000000FF;
	rv.push_back(b1);
    
    b2 = (tmp >> 16) & 0x000000FF;
    rv.push_back(b2);
    
    b3 = (tmp >> 8) & 0x000000FF;
	rv.push_back(b3);
    
    b4 = tmp & 0x000000FF;
	rv.push_back(b4);
}

int Producer1Stage::writeBElen(uint32 tmp, vector<uchar> &rv)
{
    uchar byte;
    
    if(tmp <= LEN_MAX)
    {
        byte = (tmp >> 16) & 0x000000FF;
		rv.push_back(byte);
        byte = (tmp >> 8) & 0x000000FF;
        rv.push_back(byte);
        byte = tmp & 0x000000FF;
        rv.push_back(byte);
        return 0;
    }
    else
    {
        byte = 0xFF;
        rv.push_back(byte);
        rv.push_back(byte);
        rv.push_back(byte);
        return 1;
    }
}

uint32 Producer1Stage::getMatchLen(uchar * seq1, uchar * seq2)
{
    uint32 i = 0;
    while (seq1[i] == seq2[i])
    {
        i++;
    }
    if(i>0 && (seq1[i-1] == '\0' || seq2[i-1] == '\0'))
        i--;
    return i;
}