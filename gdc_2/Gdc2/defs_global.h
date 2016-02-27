/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */
#ifndef _DEFS_H
#define _DEFS_H

#ifndef WIN32
#define _tmain main
#define _TCHAR char
#endif

#define _DEVELOPMENT_MODE

#include <string>
using namespace std;

#define MAPPER_VERSION		"0.1"

#define BOUNDARY_N_LEN		10240
#define SEPARATE_N_LEN		64

#ifdef WIN32
typedef unsigned long long uint64_t;
typedef long long int64_t;
typedef unsigned int uint32_t;
typedef int int32_t;
#endif
typedef unsigned char uchar_t;

typedef uint32_t ref_pos_t;

const char sym_code_A      = 0;
const char sym_code_C      = 1;
const char sym_code_G      = 2;
const char sym_code_T	   = 3;
const char sym_code_N_read = 4;
const char sym_code_N_ref  = 5;

const uint32_t lut_short_prefix_len = 8;
const uint32_t lut_long_prefix_len  = 12;

const uint32_t bin_total_bits     = 40;
const uint32_t bin_sub_block_bits = 3;
const uint32_t bin_in_block_bits  = 21;

const uint64_t fastq_block_size = (1 << 25);

enum class genome_t {direct, rev_comp};
enum class mapping_mode_t {first, second, all};

const uint32_t error_unknown = 255;

typedef uint64_t read_id_t;
const read_id_t empty_read_id = ~((read_id_t) 0);

const uint32_t max_fastq_rec_length = 10240;

// Index extenstions
const string EXT_DEF_IDX		 = ".trm";
const string EXT_REF_SEQ_DESC    = EXT_DEF_IDX + ".ref_seq_desc";
const string EXT_REF_SEQ_DIR     = EXT_DEF_IDX + ".ref_seq_dir";
const string EXT_REF_SEQ_RC      = EXT_DEF_IDX + ".ref_seq_rc";
const string EXT_REF_SEQ_DIR_PCK = EXT_DEF_IDX + ".ref_seq_dir_pck";
const string EXT_REF_SEQ_RC_PCK  = EXT_DEF_IDX + ".ref_seq_rc_pck";
const string EXT_SA_DIR          = EXT_DEF_IDX + ".sa_dir";
const string EXT_SA_RC           = EXT_DEF_IDX + ".sa_rc";
const string EXT_LUT_SHORT_DIR   = EXT_DEF_IDX + ".lut_short_dir";
const string EXT_LUT_SHORT_RC    = EXT_DEF_IDX + ".lut_short_rc";
const string EXT_LUT_LONG_DIR    = EXT_DEF_IDX + ".lut_long_dir";
const string EXT_LUT_LONG_RC     = EXT_DEF_IDX + ".lut_long_rc";

// Project extensions
const string EXT_DEF_PRO		 = ".trm_pro";
const string EXT_READS_BIN		 = EXT_DEF_PRO + ".reads_bin";
const string EXT_ID_STORE        = EXT_DEF_PRO + ".id_store";
const string EXT_BIN_STATS       = EXT_DEF_PRO + ".bin_stats";
const string EXT_MAPPING_GROUP   = EXT_DEF_PRO + ".mapping_group";


// Index markers
const string MARKER_DEF_IDX         = "MAP1";
const string MARKER_REF_SEQ_DESC    = MARKER_DEF_IDX + "RSDE"; 
const string MARKER_REF_SEQ_DIR     = MARKER_DEF_IDX + "RSDI"; 
const string MARKER_REF_SEQ_RC      = MARKER_DEF_IDX + "RSRC"; 
const string MARKER_REF_SEQ_DIR_PCK = MARKER_DEF_IDX + "RSDI_PCK"; 
const string MARKER_REF_SEQ_RC_PCK  = MARKER_DEF_IDX + "RSRC_PCK"; 
const string MARKER_SA_DIR          = MARKER_DEF_IDX + "SA_DIR"; 
const string MARKER_SA_RC           = MARKER_DEF_IDX + "SA_RC"; 
const string MARKER_LUT_SHORT_DIR   = MARKER_DEF_IDX + "LUT_SHORT_DIR"; 
const string MARKER_LUT_SHORT_RC    = MARKER_DEF_IDX + "LUT_SHORT_RC"; 
const string MARKER_LUT_LONG_DIR    = MARKER_DEF_IDX + "LUT_LONG_DIR"; 
const string MARKER_LUT_LONG_RC     = MARKER_DEF_IDX + "LUT_LONG_RC"; 

// Project markers
const string MARKER_DEF_PRO         = "MAP1_PRO";
const string MARKER_READS_BIN       = MARKER_DEF_PRO + "READS_BIN";
const string MARKER_ID_STORE        = MARKER_DEF_PRO + "ID_STORE";
const string MARKER_BIN_STATS       = MARKER_DEF_PRO + "BIN_STATS";
const string MARKER_MAPPING_GROUP   = MARKER_DEF_PRO + "MAPPING_GROUP";


// Running stats
const uint32_t STAT_READ_SPLITTING	  = 0x1000;

const uint32_t STAT_PRELOADING		  = 0x1800;

const uint32_t STAT_SORTING_BASE	  = 0x2000;
const uint32_t STAT_SORTING_TOTAL	  = STAT_SORTING_BASE + 0xfff;

const uint32_t STAT_BINS_WRITER_BASE  = 0x3000;
const uint32_t STAT_BINS_WRITER_TOTAL = STAT_BINS_WRITER_BASE + 0xfff;

const uint32_t STAT_MAPPING_GROUP_WRITER_BASE  = 0x4000;
const uint32_t STAT_MAPPING_GROUP_WRITER_TOTAL = STAT_MAPPING_GROUP_WRITER_BASE + 0xfff;

const uint32_t STAT_TIME_BASE	  = 0x5000;
const uint32_t STAT_TIME_TOTAL	  = STAT_TIME_BASE + 0xfff;
const uint32_t STAT_TIME_SPLIT	  = STAT_TIME_BASE + 0xffe;

const uint32_t STAT_DIR_EXACT_MATCH	          = 0x6000;
const uint32_t STAT_RC_EXACT_MATCH			  = 0x6010;
const uint32_t STAT_DIR_OR_RC_EXACT_MATCH	  = 0x6020;
const uint32_t STAT_DIR_AND_RC_EXACT_MATCH	  = 0x6030;

const uint32_t STAT_MISMATCHES_TO_DIR_BASE	  = 0x7000;
const uint32_t STAT_MISMATCHES_TO_RC_BASE	  = 0x8000;

const uint32_t STAT_SELECTED_READS_BASE		  = 0X9000;

const uint32_t STAT_TIME_THR_RES_WRITER        = 0x10000;
const uint32_t STAT_TIME_THR_FASTQ_READER      = 0x10001;
const uint32_t STAT_TIME_THR_SPLITTER          = 0x10002;
const uint32_t STAT_TIME_THR_BIN_WRITER_BASE   = 0x11000;
const uint32_t STAT_TIME_THR_BIN_READER_BASE   = 0x12000;
const uint32_t STAT_TIME_THR_MAPPING_CORE_BASE = 0x13000;
const uint32_t STAT_TIME_THR_FASTQ_READER_PP   = 0x14000;
const uint32_t STAT_TIME_THR_RES_READER        = 0x14001;
const uint32_t STAT_TIME_THR_SAM_GENERATOR     = 0x14002;
const uint32_t STAT_TIME_THR_SAM_SORTING       = 0x14003;
const uint32_t STAT_TIME_THR_SAM_PROCESSING    = 0x14004;

const uint32_t STAT_PUSH_RESULTS               = 0x15000;
const uint32_t STAT_PUSH_RESULTS_UNQ           = 0x15001;
const uint32_t STAT_PUSH_RESULTS_SEND_BYTES    = 0x15002;
const uint32_t STAT_PUSH_RESULTS_ADDED_BYTES   = 0x15003;

# define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))	
# define MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))	


#ifdef _DEBUG
#define A_memcpy	memcpy
#define A_memset	memset
#endif

#endif