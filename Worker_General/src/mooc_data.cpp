#include<iostream>
#include<fstream>
#include<map>
#include<string>
#include<stdlib.h>
#include<math.h>
#include<string>
#include<sstream>
#include<cassert>
#include"worker.h"
#include"pool.h"
#include"test.h"
#include<Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::JacobiSVD;


#define f(i,n) for(int i=0;i<n;i++)

using namespace std;

void get_mooc_data(int &arity, int& n, int& wn, bool** &attempted, int** &ans, int* &true_ans){
	arity = 6;
	ifstream data;
	string SubId,AsgnId,GradId,GType,Grade,PeerGrade,StaffGrade;
	map<string,int> WorkerMap;
	map<string,int> TaskMap;
	data.open("MOOC_Data/evaluation_data.csv");
	n = 0;
	while(getline(data,SubId,',')){
		getline(data,AsgnId,',');
		getline(data,GradId,',');
		getline(data,GType,',');
		getline(data,Grade,',');
		getline(data,PeerGrade,',');
		getline(data,StaffGrade,'\n');
		if(GType=="\"staff\""){
			map<string,int>::iterator it = TaskMap.find(SubId);
			if(it==TaskMap.end()){
				TaskMap.insert(pair<string,int>(SubId,n));
				n++;
			}
		}
	}
	data.close();
	true_ans = new int[n];
	
	data.open("MOOC_Data/evaluation_data.csv");
	wn = 0;
	while(getline(data,SubId,',')){
		getline(data,AsgnId,',');
		getline(data,GradId,',');
		getline(data,GType,',');
		getline(data,Grade,',');
		getline(data,PeerGrade,',');
		getline(data,StaffGrade,'\n');
		map<string,int>::iterator it = TaskMap.find(SubId);
		if(it!=TaskMap.end()){
			if(GType=="\"student\""){
				map<string,int>::iterator it = WorkerMap.find(GradId);
				if(it==WorkerMap.end()){
					WorkerMap.insert(pair<string,int>(GradId,wn));
					wn++;
				}
			}
			else{
				int TaskNo = it->second;
				true_ans[TaskNo] = Grade[1] - '0';
			}
		}
		
	}
	data.close();

	attempted = new bool*[wn];
	ans = new int*[wn];
	f(i,wn){
		attempted[i] = new bool[n];
		f(j,n){
			attempted[i][j] = false;
		}
		ans[i] = new int[n];
	}

	data.open("MOOC_Data/evaluation_data.csv");
	while(getline(data,SubId,',')){
		getline(data,AsgnId,',');
		getline(data,GradId,',');
		getline(data,GType,',');
		getline(data,Grade,',');
		getline(data,PeerGrade,',');
		getline(data,StaffGrade,'\n');
		map<string,int>::iterator it2 = TaskMap.find(SubId);
		if(it2 != TaskMap.end() && GType=="\"student\""){
			map<string,int>::iterator it1 = WorkerMap.find(GradId);
			int WorkerNo = it1->second;
			int TaskNo = it2->second;
			attempted[WorkerNo][TaskNo] = true;
			ans[WorkerNo][TaskNo] = Grade[1] - '0';
		}
	}
	data.close();
	return;
}

void get_mooc_data_binary(int& n, int& wn, bool** &attempted, bool** &ans, bool* &true_ans){
	ifstream data;
	string SubId,AsgnId,GradId,GType,Grade,PeerGrade,StaffGrade;
	map<string,int> WorkerMap;
	map<string,int> TaskMap;
	data.open("MOOC_Data/evaluation_data.csv");
	n = 0;
	while(getline(data,SubId,',')){
		getline(data,AsgnId,',');
		getline(data,GradId,',');
		getline(data,GType,',');
		getline(data,Grade,',');
		getline(data,PeerGrade,',');
		getline(data,StaffGrade,'\n');
		if(GType=="\"staff\""){
			map<string,int>::iterator it = TaskMap.find(SubId);
			if(it==TaskMap.end()){
				TaskMap.insert(pair<string,int>(SubId,n));
				n++;
			}
		}
	}
	data.close();
	true_ans = new bool[n];
	
	data.open("MOOC_Data/evaluation_data.csv");
	wn = 0;
	while(getline(data,SubId,',')){
		getline(data,AsgnId,',');
		getline(data,GradId,',');
		getline(data,GType,',');
		getline(data,Grade,',');
		getline(data,PeerGrade,',');
		getline(data,StaffGrade,'\n');
		map<string,int>::iterator it = TaskMap.find(SubId);
		if(it!=TaskMap.end()){
			if(GType=="\"student\""){
				map<string,int>::iterator it = WorkerMap.find(GradId);
				if(it==WorkerMap.end()){
					WorkerMap.insert(pair<string,int>(GradId,wn));
					wn++;
				}
			}
			else{
				int TaskNo = it->second;
				true_ans[TaskNo] = (Grade[1] - '0') > 2;
			}
		}
		
	}
	data.close();
	
	attempted = new bool*[wn];
	ans = new bool*[wn];
	f(i,wn){
		attempted[i] = new bool[n];
		f(j,n){
			attempted[i][j] = false;
		}
		ans[i] = new bool[n];
	}

	data.open("MOOC_Data/evaluation_data.csv");
	while(getline(data,SubId,',')){
		getline(data,AsgnId,',');
		getline(data,GradId,',');
		getline(data,GType,',');
		getline(data,Grade,',');
		getline(data,PeerGrade,',');
		getline(data,StaffGrade,'\n');
		map<string,int>::iterator it2 = TaskMap.find(SubId);
		if(it2 != TaskMap.end() && GType=="\"student\""){
			map<string,int>::iterator it1 = WorkerMap.find(GradId);
			int WorkerNo = it1->second;
			int TaskNo = it2->second;
			attempted[WorkerNo][TaskNo] = true;
			ans[WorkerNo][TaskNo] = (Grade[1] - '0') > 2;
		}
	}
	data.close();
	return;
}

