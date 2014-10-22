#ifndef WORDSIMDATA_H
#define WORDSIMDATA_H

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

void get_wordsim_data(int &arity, int& n, int& wn, bool** &attempted, int** &ans, int* &true_ans);

#endif
