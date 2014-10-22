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
#include<Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::JacobiSVD;

#define f(i,n) for(int i=0;i<n;i++)

using namespace std;

void get_imana_data(int& n, int& wn, bool** &attempted, bool** &ans, bool* &true_ans){
	double readIn;
	wn = 402;
	n = 60;
	true_ans = new bool[n];
	attempted = new bool*[wn];
	ans = new bool*[wn];
	f(i,wn){
		attempted[i] = new bool[n];
		f(j,n){
			attempted[i][j] = true;
		}
		ans[i] = new bool[n];
	}
	
	ifstream data;
	data.open("imana_data.txt");
	f(i,wn){
		f(j,n){
			if(j==36){
				data >> readIn >> readIn >> readIn >> readIn; //skip data for the 4 counting tasks
			}
			data >> readIn;
			ans[i][j] = readIn;
		}
	}
	
	//true ans being set as majority. so not good idea to use this in majority comparison
	f(j,n){
		int count = 0;
		f(i,wn){
			if(ans[i][j]){
				count++;
			}
		}
		if(count > wn/2){
			true_ans[j] = true;
		}
		else{
			true_ans[j] = false;
		}
	} 

}
