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
	int arity;
	bool* btrue_ans;
	int* true_ans;
	bool** bans;
	int ** ans;
	bool** attempted;
	get_wordsim_data(arity, n, wn, attempted, ans, true_ans);
	arity = 1 + (arity-1)/6;
	f(j,n){
		true_ans[j] /= 6;
		f(i,wn){
			ans[i][j] /= 6;
		}
	}
	/*
	ans = new int*[wn];
	f(i,wn){
		ans[i] = new int[n];
		f(j,n){
			ans[i][j] = bans[i][j];
		}
	}
	true_ans = new int[n];
	f(j,n){
		true_ans[j] = btrue_ans[j];
	}
	*/
	MatrixXd werr[wn];
	MatrixXd terr[n];
	MatrixXd wsubtot[wn];
	MatrixXd s(arity,1);
	f(a1,arity){
		s(a1,0) = 0;
	}
	f(i,wn){
		wsubtot[i].resize(arity,1);
		werr[i].resize(arity,arity);
	}
	f(j,n){
		terr[j].resize(arity,1);
	}
	int wtot[wn];
	int ttot[n];
	f(i,wn){
		f(a1,arity){
			wsubtot[i](a1,0) = 0;
			f(a2,arity){
				werr[i](a1,a2) = 0;
			}
		}
		wtot[i] = 0;
	}
	f(j,n){
		f(a1,arity){
			terr[j](a1,0) = 0;
		}
		ttot[j] = 0;
		s(true_ans[j],0) += 1;
	}
	
	f(a1,arity){
		s(a1,0) /= n;
	}
	
	f(i,wn){
		f(j,n){
			if(attempted[i][j]){
				wtot[i] ++;
				wsubtot[i](true_ans[j],0) += 1;
				ttot[j] ++;
				terr[j](ans[i][j],0) += 1;
				werr[i](true_ans[j],ans[i][j]) += 1;
			}
		}
	}
	
	f(i,wn){
		f(a1,arity){
			f(a2,arity){
				werr[i](a1,a2) /= wsubtot[i](a1,0);
			}
		}
	}
	f(j,n){
		terr[j] /= ttot[j];
	}
	double conf = 0.8;
	MatrixXd est[wn],epsilon[wn];
	
	int Iters = 200;
	int** tempans = new int* [wn];
	bool** tempattempted = new bool* [wn];	
	bool taken[wn];
	int tempwn=3;
	int tempWorkers[tempwn];
	int taskThresh = 30;
	double uncertainty;
	double total;
	for(conf = 0.1; conf < 0.99; conf += 0.05){
		double corrcount = 0;
		total = 0;
		uncertainty = 0;
		f(numIter,Iters){
			int select[3];
			bool workersFound = false;
			while(!workersFound){
				f(i,3){
					select[i] = 0;
				}
				f(i,wn){
					taken[i] = false;
				}
			
				f(i,tempwn){
					bool done = false;
					while(!done){
						int workerNo = rand()%wn;
						if(taken[workerNo] ==  false){
							done = true;
							taken[workerNo] = true;
							tempWorkers[i] = workerNo;
							tempans[i] = ans[workerNo];
							tempattempted[i] = attempted[workerNo];
						}
					}
				}
				int taskcount = 0;
				f(j,n){
					bool skip = false;
					f(i,tempwn){
						if(!attempted[tempWorkers[i]][j]){
							skip = true;
						}
					}
					if(!skip){
						taskcount++;select[true_ans[j]]++;
					}
				}
				if(taskcount >= taskThresh && select[1] >= 1){
					workersFound = true;
				}
			}//cout << select[0] << '\t' << select[1] << '\t' << select[2] << endl;
			Three_Worker_nonreg_Kary(arity, tempattempted, tempans, n, est, epsilon, conf);
			
			MatrixXd sel = MatrixXd::Zero(arity,1);
			f(i,tempwn){
				f(a1,arity){
					f(a2,arity){
						sel(a1,0) += est[i](a1,a2);	
					}
				}
			}
			f(a1,arity){
				sel(a1,0) /= tempwn;
				sel(a1,0) = sel(a1,0)*sel(a1,0);
			}
			
			f(i,tempwn){//cout << est[i] << '\t' << epsilon[i] << '\t' << werr[i] << endl;
				MatrixXd temp = sel.cwiseSqrt().asDiagonal() * werr[tempWorkers[i]];
				
				f(a1,arity){
					f(a2,arity){if(epsilon[i](a1,a2) > sel(a1,0)){epsilon[i](a1,a2) = sel(a1,0);}
						if(epsilon[i](a1,a2)<=0 || epsilon[i](a1,a2)>=0){
							
							if(fabs(est[i](a1,a2) - temp(a1,a2)) <= epsilon[i](a1,a2)){/////HAVE to Multiple by selectivity matrix OR normalise rows
								corrcount += 1;
							}
							else{
								if(conf > 0.89){
									//cout << "new\n"<< temp << endl << est[i]-epsilon[i] << endl << est[i]+epsilon[i] << endl;
								}
							}
							uncertainty += epsilon[i](a1,a2)/sel(a1,0);	
						}
						total += 1;
					}
				}
			}
		}
		//corrcount /= (tempwn*Iters*arity*arity);
		//uncertainty /= (tempwn*Iters*arity*arity);
		corrcount /= total;
		uncertainty /= total;
		cout << conf << '\t' <<  corrcount << '\t' << uncertainty << endl;//cerr << conf << endl;
	}
	
	return 0;
}
