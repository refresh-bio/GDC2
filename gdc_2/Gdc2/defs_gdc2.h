/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef DEFS_GDC2_H
#define DEFS_GDC2_H

enum t_working_mode {wm_compress, wm_decompress, wm_list};
#define POS_MAX ((unsigned int)(1<<31)-3)
#define LEN_MAX ((1<<24)-1)

#define BUF_SIZE 50000000

#endif