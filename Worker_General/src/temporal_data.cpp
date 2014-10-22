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

void get_temporal_data(int& n, int& wn, bool** &attempted, bool** &ans, bool* &true_ans){
	n = 0;
	wn = 0;
	ifstream tempdata;
	tempdata.open("Annotation Data/temp.standardized.tsv");
	string ss;
	getline(tempdata,ss,'\n');
	map<string,int> WorkerMap;
	map<string,int> TaskMap;
 	map<string,int>::iterator it1,it2;
	string ss1,ss2,ss3,ss4,ss5;
	while(tempdata >> ss1){
		tempdata >> ss2 >> ss3 >> ss4 >> ss5;
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
	tempdata.close();
	
	true_ans = new bool[n];
	attempted = new bool*[wn];
	ans = new bool*[wn];
	f(i,wn){
		attempted[i] = new bool[n];
		f(j,n){
			attempted[i][j] = false;
		}
		ans[i] = new bool[n];
	}
	
	tempdata.open("Annotation Data/temp.standardized.tsv");
	getline(tempdata,ss,'\n');
	while(tempdata >> ss1){
		tempdata >> ss2 >> ss3 >> ss4 >> ss5;
		it1 = WorkerMap.find(ss2);
		it2 = TaskMap.find(ss3);
		int WorkerNo = it1->second;
		int TaskNo = it2->second;
		true_ans[TaskNo] = (ss5=="1");
		attempted[WorkerNo][TaskNo] = true;
		ans[WorkerNo][TaskNo] = (ss4=="1");
	}
	tempdata.close();
}

void get_temporal_data_regular(int& n, int& wn, bool** &ans, bool* &true_ans){
	n = 0;
	wn = 0;
	ifstream tempdata;
	tempdata.open("Annotation Data/temp.standardized.tsv");
	string ss;
	getline(tempdata,ss,'\n');
	map<string,int> WorkerMap;
	map<string,int> TaskMap;
 	map<string,int>::iterator it1,it2;
	string ss1,ss2,ss3,ss4,ss5;
	while(tempdata >> ss1){
		tempdata >> ss2 >> ss3 >> ss4 >> ss5;
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
	tempdata.close();
	
	bool * true_tempans = new bool[n];
	bool ** tempattempted = new bool*[wn];
	bool ** tempans = new bool*[wn];
	int counts[wn];
	f(i,wn){
		tempattempted[i] = new bool[n];
		f(j,n){
			tempattempted[i][j] = false;
		}
		tempans[i] = new bool[n];
		counts[i] = 0;
	}
	
	tempdata.open("Annotation Data/temp.standardized.tsv");
	getline(tempdata,ss,'\n');
	while(tempdata >> ss1){
		tempdata >> ss2 >> ss3 >> ss4 >> ss5;
		it1 = WorkerMap.find(ss2);
		it2 = TaskMap.find(ss3);
		int WorkerNo = it1->second;
		int TaskNo = it2->second;
		true_tempans[TaskNo] = (ss5=="1");
		tempattempted[WorkerNo][TaskNo] = true;
		counts[WorkerNo]++;
		tempans[WorkerNo][TaskNo] = (ss4=="1");
	}
	tempdata.close();

	int taskThresh = 250;
	int workers[wn];
	int tasks[n];
	int newn = 0;
	int newwn = 0;
	
	f(i,wn){//make list of workers who have done >= 100 tasks
		if(counts[i] >= taskThresh){
			workers[newwn] = i;
			newwn++;
		}
	}
	f(j,n){//select tasks which have been attempted by every worker in above list
		bool select = true;
		f(i,newwn){
			if(!tempattempted[workers[i]][j]){
				select = false;
			}
		}
		
		if(select){
			tasks[newn] = j;
			newn++;
		}
	}
	
	n = newn;
	wn = newwn;
	
	true_ans = new bool[n];
	ans = new bool*[wn];
	f(i,wn){
		ans[i] = new bool[n];
	}
	f(i,wn){
		f(j,n){
			ans[i][j] = tempans[workers[i]][tasks[j]];
		}
	}
	f(j,n){
		true_ans[j] = true_tempans[tasks[j]];
	}
	return;
}
