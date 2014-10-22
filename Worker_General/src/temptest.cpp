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
	int n;
	int wn;
	bool* true_ans;
	bool** attempted;
	bool** ans;
	get_temporal_data(n, wn, attempted, ans, true_ans);

	int wtots[wn];
	double wacc[wn];
	int ttots[n];
	double tacc[n];
	
	f(i,wn){
		wtots[i] = 0;
		wacc[i] = 0;
	}
	f(j,n){
		ttots[j] = 0;
		tacc[j] = 0;
	}
	
	f(i,wn){
		f(j,n){
			if(attempted[i][j]){
				wtots[i]++;
				ttots[j]++;
				if(ans[i][j] == true_ans[j]){
					wacc[i] += 1;
					tacc[j] += 1.0;
				}
			}
		}
	}
	
	int numCats = 11;
	double t_hist[numCats];
	double w_hist[numCats];
	
	double histacc[numCats][numCats];
	int histtot[numCats][numCats];
	
	f(Windex,numCats){
		w_hist[Windex] = 0;
	}
	f(Tindex,numCats){
		t_hist[Tindex] = 0;	
	}
	
	f(Windex,numCats){
		f(Tindex,numCats){
			histtot[Windex][Tindex] = 0;
			histacc[Windex][Tindex] = 0;
		}
	}
	
	f(i,wn){
		wacc[i] /= wtots[i];
		w_hist[(int)(wacc[i]*(numCats-1))] += 1;
	}
	f(j,n){
		tacc[j] /= ttots[j];
		t_hist[(int)(tacc[j]*(numCats-1))] += 1;
	}
	
	f(i,wn){
		//cout << wacc[i] << '\t' << wtots[i] << endl;
	}
	f(j,n){
		//cout << tacc[j] << '\t' << ttots[j] << endl;
	}
	
	f(i,numCats){
		cout << i << '\t' << w_hist[i] << '\t' << t_hist[i] << endl;
	}
	cout << endl;
	
	
	f(i,wn){
		f(j,n){
			if(attempted[i][j]){
				int Windex = (int)(wacc[i]*(numCats-1));
				int Tindex = (int)(tacc[j]*(numCats-1));
				histtot[Windex][Tindex]++;
				if(ans[i][j] == true_ans[j]){
					histacc[Windex][Tindex] += 1;
				}
			}
		}
	}
	
	double func1[numCats][numCats],func2[numCats][numCats], func3[numCats][numCats], func4[numCats][numCats];
	
	f(Windex,numCats){
		f(Tindex,numCats){
			histacc[Windex][Tindex] /= histtot[Windex][Tindex];
			func1[Windex][Tindex] = log(1/histacc[Windex][Tindex] - 1);
			func2[Windex][Tindex] = log(2*histacc[Windex][Tindex] - 1);
			func4[Windex][Tindex] = log(-log(1/histacc[Windex][Tindex] - 1));
			func3[Windex][Tindex] = log(-log(histacc[Windex][Tindex]));
		}
	}
	
	
	for(int Windex=numCats/2; Windex < numCats; Windex++){
		for(int Tindex=numCats/2;Tindex < numCats;Tindex++){
			double f = (func3[Windex][Tindex]-func3[Windex-1][Tindex]-func3[Windex][Tindex-1]+func3[Windex-1][Tindex-1])/(func3[Windex][Tindex]+func3[Windex-1][Tindex]+func3[Windex][Tindex-1]+func3[Windex-1][Tindex-1]);
			printf("%.2f,%d   \t",f,histtot[Windex][Tindex]);
			//cout << histacc[Windex][Tindex] << "," << histtot[Windex][Tindex] << "  ";
		}
		cout << endl;
	}
	cout << endl;
	
	double conf;
	return 0;
}
