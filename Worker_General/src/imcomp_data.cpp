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

void get_imcomp_data(int& n, int& wn, bool** &attempted, bool** &ans, bool* &true_ans){
	ifstream finGold("imcomp_gym_Gold.txt");
	ifstream finRes("imcomp_gym_data.txt");
	ifstream finQue("imcomp_gym_questions.txt");
	const int numImages = 42;
	n = 48;
	wn = 19;
	attempted = new bool*[wn];
	ans = new bool*[wn];
	f(i,wn){
		attempted[i] = new bool[n];
		f(j,n){
			attempted[i][j] = true;
		}
		ans[i] = new bool[n];
	}
	
	
	map<int,int> Person;
	int personId = 0;
	string input;
	srand ( time(NULL) );
	while(getline(finGold,input,'\n')){
		stringstream ss(input);
		int photoId;
		while(ss >> photoId){
			Person.insert(pair<int,int>(photoId,personId));
		}
		personId++;
	}
	finGold.close();
	int q1[n],q2[n];
	true_ans = new bool[n];
	f(i,n){
		int im1,im2;
		finQue >> im1 >> im2;
		q1[i] = im1;
		q2[i] = im2;
		true_ans[i] = ((Person.find(im1))->second == (Person.find(im2))->second);
	}
	finQue.close();
	f(i,wn){
		f(j,n){
			int inp;
			finRes >> inp;
			ans[i][j] = (2-inp);
		}
	}
	finRes.close();
	
	
}
