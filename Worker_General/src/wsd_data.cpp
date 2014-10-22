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

void get_wsd_data(int &arity, int& n, int& wn, bool** &attempted, int** &ans, int* &true_ans){
	arity = 3;
	n = 0;
	wn = 0;
	ifstream wsddata;
	wsddata.open("Annotation Data/wsd.standardized.tsv");
	string ss;
	getline(wsddata,ss,'\n');
	map<string,int> WorkerMap;
	map<string,int> TaskMap;
 	map<string,int>::iterator it1,it2;
	string ss1,ss2,ss3,ss4,ss5;
	while(wsddata >> ss1){
		wsddata >> ss2 >> ss3 >> ss4 >> ss5;
		it1 = WorkerMap.find(ss2);
		if(it1==WorkerMap.end()){
			WorkerMap.insert(pair<string,int>(ss2,wn));
			wn++;
		}
		it2 = TaskMap.find(ss3);
		if(it2==TaskMap.end()){
			TaskMap.insert(pair<string,int>(ss3,n));
			n++;
		}
	}
	wsddata.close();
	
	true_ans = new int[n];
	attempted = new bool*[wn];
	ans = new int*[wn];
	f(i,wn){
		attempted[i] = new bool[n];
		f(j,n){
			attempted[i][j] = false;
		}
		ans[i] = new int[n];
	}
	
	wsddata.open("Annotation Data/wsd.standardized.tsv");
	getline(wsddata,ss,'\n');
	while(wsddata >> ss1){
		wsddata >> ss2 >> ss3 >> ss4 >> ss5;
		it1 = WorkerMap.find(ss2);
		it2 = TaskMap.find(ss3);
		int WorkerNo = it1->second;
		int TaskNo = it2->second;
		true_ans[TaskNo] = atoi(ss5.c_str())-1;
		attempted[WorkerNo][TaskNo] = true;
		ans[WorkerNo][TaskNo] = atoi(ss4.c_str())-1;
	}
	wsddata.close();
}

