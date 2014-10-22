#include<iostream>
#include<stdlib.h>
#include"worker.h"
#include"pool.h"
#include<string>
#include<math.h>
#include"test.h"
#include"imcomp_data.h"
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<cassert>


#define f(i,n) for(int i=0;i<n;i++)

using namespace std;

int main(){
	int n,wn;
	bool* true_ans;
	bool** ans;
	bool** attempted;
	get_imcomp_data(n, wn, attempted, ans, true_ans);
		
	f(i,n){
		int corr = 0;
		f(j,wn){
			if(ans[j][i] == true_ans[i]){
				corr++;
			}
		}	
		//cout << i << '\t' << q1[i] << ',' << q2[i] << '\t' << ((double)corr)/wn << endl;
		double acc = ((double)corr)/wn;
		if(acc < 0.5){
			//cout << i << '\t' << q1[i] << ',' << q2[i] << '\t' << ((double)corr)/wn << endl;	
		}
	}

	double Wacc[wn];

	f(j,wn){
		int corr = 0;
		f(i,n){
			if(ans[j][i] == true_ans[i]){
				corr++;
			}
		}
		Wacc[j] = ((double)corr)/n;
		//cout << j << '\t' << acc << endl;	
	}
	
	
	for(double conf = 0.1; conf < 1.0; conf+=0.05){
		int matches = 0;
		int total_comp = 0;
		f(w1,wn){for(int w2=w1+1;w2<wn;w2++){for(int w3=w2+1;w3<wn;w3++){
			bool found = false;
			int testW = 3;
			assert(testW <= wn);
			int w[testW];
			double est[testW],epsilon[testW];
			/*
			f(i,testW){
				int found = false;
				while(!found){
					found = true;
					w[i] = rand()%wn;
					f(j,i){
						if(w[j]==w[i]){
							found = false;
						}
					}
				}
			}
			*/
			w[0] = w1;w[1] = w2;w[2] = w3;
			
			bool ** answers;
			answers = new bool*[testW];
			f(i,testW){
				answers[i] = new bool[n];
			}
			f(i,testW){
				f(j,n){
					answers[i][j] = ans[w[i]][j];
				}
			}
			Three_Worker(answers, n, est, epsilon,conf);
			f(i,testW){
				total_comp++;
				//cout << "(" << est[i]-epsilon[i]<<", "<<est[i]+epsilon[i]<<")\t"<<(1-Wacc[w[i]])<<"\t";
				if(fabs(est[i]-(1-Wacc[w[i]])) < epsilon[i]){
					matches++;
				}
			}
		}}}
		cout << conf << '\t' << ((double)matches)/total_comp<<endl;
	}
	
	
	/*
		int j1=rand()%wn,j2=rand()%wn;
		int corr1 = 0,corr2 = 0, corr12 = 0;
		f(i,n){
			if(ans[j1][i] == true_ans[i]){
				corr1++;
			}
			if(ans[j2][i] == true_ans[i]){
				corr2++;
			}
			if(ans[j1][i] == ans[j2][i]){
				corr12++;
			}
		}
		double acc1 = ((double)corr1)/n, acc2 = ((double)corr2)/n, acc12 = ((double)corr12)/n;
		cout << acc1 << '\t' << acc2 << '\t' << (acc1*acc2+(1-acc1)*(1-acc2)) << '\t' << acc12 << endl;
	*/
	
	/*
	double est1[wn],epsilon1[wn];
	double est2[wn],epsilon2[wn];
	const int numTrain = n/2;
	const int numTest = n-numTrain;
	bool ans1[numTrain];
	bool correct1[numTrain];
	bool ans2[numTest];
	bool correct2[numTest];
	double conf = 0.9;
	int numIters = wn;
	ofstream fout("MTurk_diff");
	for(conf = 0.1;conf<1;conf+=0.05){
		int wrongcount = 0;
		f(iter,numIters){
			int i = iter;
			//int hardcount = 0;
			int cTrain=0,cTest=0;
			int corrTrain = 0, corrTest = 0;
			f(j,n){
				if(false){
					continue;
				}
				bool toTrain; 
				if(cTrain==numTrain){
					toTrain = false;  
				}
				else if(cTest==numTest){
					toTrain = true;
				}
				else{
					toTrain = true;//((rand()%numTasks)<numTrain); 
				}
				if(toTrain){
					ans1[cTrain] = ans[i][j];
					correct1[cTrain] = true_ans[j];
					if(ans1[cTrain] == correct1[cTrain])
						corrTrain++;
					cTrain++;
				}
				else{
					ans2[cTest] = ans[i][j];
					correct2[cTest] = true_ans[j];
					if(ans2[cTest] == correct2[cTest])
						corrTest++;
					cTest++;
				}
			}
			
			gold_standard_test(ans1,correct1,numTrain, est1[i], epsilon1[i], conf);
			gold_standard_test(ans2,correct2,numTest , est2[i], epsilon2[i], conf);
			
			
			if(fabs(est1[i]-est2[i]) > (epsilon1[i]+epsilon2[i])){
				wrongcount++;
			}
			
			if(conf < 0.11 && conf > 0.09){
				fout << est1[i]-est2[i] << endl;
			}
		}
		cout << conf << '\t' << ((double)wrongcount)/numIters<<endl;
	}	
	*/
	
	
}
