#ifndef MOOCDATA_H
#define MOOCDATA_H

#include<iostream>
#include<stdlib.h>
#include"worker.h"
#include"pool.h"
#include<string>
#include<math.h>
#include"test.h"
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<cassert>

void get_mooc_data(int& arity, int& n, int& wn, bool** &attempted, int** &ans, int* &true_ans);

void get_mooc_data_binary(int& n, int& wn, bool** &attempted, bool** &ans, bool* &true_ans);

void get_mooc_data_regular(int& n, int& wn, bool** &ans, bool* &true_ans);

#endif
