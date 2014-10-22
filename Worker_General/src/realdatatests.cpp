#include<iostream>
#include<stdlib.h>
#include"worker.h"
#include"pool.h"
#include<string>
#include<math.h>
#include"test.h"
#include"datasets.h"
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

int main(){
	srand(time(0));
	int n,wn;
	bool* true_ans;
	bool** ans;
	bool** attempted;
	get_imcomp_data(n, wn, attempted, ans, true_ans);
	int majcount = 0;
	bool maj_ans[n];
	f(j,n){
		int no = 0;
		int yes = 0;
		f(i,wn){
			if(attempted[i][j]){
				if(ans[i][j]){
					yes++;
				}	
				else{
					no++;
				}
			}
		}
		maj_ans[j] = yes >= no;
	}
	
	double werrmaj[wn];
	double werr[wn];
	double terr[n];
	int wtot[wn];
	int ttot[n];
	
	int count[wn][wn];
	double q[wn][wn];
	double estq[wn][wn];
	
	f(j1,wn){
		f(j2,wn){
			count[j1][j2] = 0;
			q[j1][j2] = 0.0;
		}
	}
	
	f(i,n){
		f(j1,wn){
			if(!attempted[j1][i]){
				continue;
			}
			f(j2,wn){
				if(attempted[j2][i]){
					count[j1][j2]++;
					if(ans[j1][i] == ans[j2][i]){
						q[j1][j2] += 1;
					}
				}
			}
		}
	}
	
	f(j1,wn){
		f(j2,wn){
			estq[j1][j2] = q[j1][j2]/count[j1][j2];
		}
	}
	
	
	f(i,wn){
		werrmaj[i] = 0;
		werr[i] = 0;
		wtot[i] = 0;
	}
	f(j,n){
		terr[j] = 0;
		ttot[j] = 0;
	}
	f(i,wn){
		f(j,n){
			if(attempted[i][j]){
				wtot[i] ++;
				ttot[j] ++;
				if(ans[i][j] != true_ans[j]){
					werr[i] += 1;
					terr[j] += 1;
				}
				if(ans[i][j] != maj_ans[j]){
					werrmaj[i] += 1;
				}
			}
		}
	}
	f(i,wn){
		werr[i] /= wtot[i];
		werrmaj[i] /= wtot[i];
	}
	f(j,n){
		terr[j] /= ttot[j];
	}
	
	f(i,wn){
		if(werr[i] > 0.25){
			//cout << werr[i] << '\t' << werrmaj[i] << endl;
		}		
	}
	
	double werrEasy[wn];
	double werrHard[wn];
	int wtotEasy[wn];
	int wtotHard[wn];
	f(i,wn){
		werrEasy[i] = 0;
		werrHard[i] = 0;
		wtotEasy[i] = 0;
		wtotHard[i] = 0;
	}
	bool easy[n];
	f(j,n){
		if(terr[j] < 0.9){///Easy/Hard threshold set here 
			easy[j] = true;
		}	
		else{//cout << terr[j] << endl;
			easy[j] = false;
		}
	}
	f(i,wn){
		f(j,n){
			if(attempted[i][j]){
				if(easy[j]){
					wtotEasy[i]++;
				}
				else{
					wtotHard[i]++;
				}
				if(ans[i][j] != true_ans[j]){
					if(easy[j]){
						werrEasy[i] += 1;
					}
					else{
						werrHard[i] += 1;
					}
				}
			}
		}
	}
	f(i,wn){
		werrEasy[i] /= wtotEasy[i];
		werrHard[i] /= wtotHard[i];
	}
	
	
	double conf = 0.8;
	double est[wn],epsilon[wn];
	
	int Iters = 50;
	bool** tempans = new bool* [wn];
	bool** tempattempted = new bool* [wn];
	f(i,wn){
		tempans[i] = new bool[n];
		tempattempted[i] = new bool[n];
	}
	bool taken[wn];
	int tempwn=11;
	int tempWorkers[wn];
	int tempn;
	int tempTasks[n];
	double uncertainty;
	int total;
	
	double density = 0.8;
	
	for(conf = 0.1; conf < 1; conf += 0.05){
		total = 0;
		double corrcount = 0;
		uncertainty = 0;
		f(numIter,Iters){
			f(i,wn){
				taken[i] = false;
			}
			tempwn = 0;
			f(i,wn){
				if(werrmaj[i] < 1.4){
					tempWorkers[tempwn] = i;
					tempwn++;
				}
			}
			
			/*
			f(i,tempwn){
				bool done = false;
				while(!done){
					int workerNo = i;//rand()%wn;
					if(taken[workerNo] ==  false){
						done = true;
						taken[workerNo] = true;
						tempWorkers[i] = workerNo;
					}
				}
			}
			*/
			
			tempn = 0;
			f(j,n){
				if(true || easy[j]){
					tempTasks[tempn] = j;
					tempn++;
				}
			}
			
			f(i,tempwn){
				f(j,tempn){
					tempans[i][j] = ans[tempWorkers[i]][tempTasks[j]];
					tempattempted[i][j] = attempted[tempWorkers[i]][tempTasks[j]];
					//below is for imcomp only, since it's regular, we are de-regularizing it
					if(rand()%100 >= density*100){
						tempattempted[i][j] = false;
					}
				}
			}
			//cout << tempwn << '\t' << tempn << " out of " << n << endl;
			
			N_Worker_nonreg(tempans, tempattempted, tempwn, tempn, est, epsilon, conf, false); 
			//N_Worker_nonreg(nans, nattempted, wn, nn, est, epsilon, conf); 
			f(i,tempwn){//cout << est[i] << '\t' << epsilon[i] << '\t' << werr[i] << endl;
				//cout << tempWorkers[i] << '\t' << werr[tempWorkers[i]] << "\t" << est[i] << '\t' << epsilon[i] << endl;
				
				if(est[i]+epsilon[i]>=0 || est[i]+epsilon[i] < 0){
					total++;//that means est and epsilon are not nan.
					uncertainty += epsilon[i];	
					if(fabs(est[i] - werrEasy[tempWorkers[i]]) <= epsilon[i]){
						corrcount += 1;
					}
					else{
						if(conf > 0.9){
							//cout << "("<<est[i]-epsilon[i]<<", "<<est[i]+epsilon[i]<<")\t"<<werrEasy[tempWorkers[i]]<<endl;
						}
					}
				}
			}
		}
		uncertainty /= total;
		corrcount /= total;
		cout << conf << '\t' << corrcount << '\t' << uncertainty << endl;//cerr << conf << endl;
	}
	
	return 0;
}
