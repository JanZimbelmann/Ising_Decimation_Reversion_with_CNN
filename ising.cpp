////////////////////////////////
//Code written by Jan Zimbelmann
////////////////////////////////
//Monte Carlo Simulation for  //
//calculating obvservables and//
//creating spin configurations//
////////////////////////////////

//libraries
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <ctime>
#include <sys/stat.h>

using namespace std;

///////////////////////////////////////
//initiate system boundary conditions//
///////////////////////////////////////
//costumization: 
//
//step > 0 does renormalize the beta to
//fit the according step of system
//enlargement
///////////////////////////////////////

//boundary conditions
const int L = 32; //1D system length
int N = L; //spin number
double J = 1; //coupling constant
double Beta[10] = {1,2,3,4,5,6,7,8,9,10}; //inverse temperature array
double factor = 0.3; //beta is later on multiplied by factor
int spin[L]; //spin array
double r = 0; // uniform random number
int x = 0; // random position
int steps = 0; //enlarging steps, part of a future code

//observables and flipping rules
double condition = 0;
int beforeE = 0;
int afterE = 0;
int difE = 0;
int E = 0;
int M = 0;
double savedE = 0;
double savedM = 0;
double savedG = 0;
int distG = 3;
int data_points = 10;

//saving conditions (requires manual tweeking if counter is 'not reached')
int sweep = 0;
int start = N * 1000;
int outSize = 20000;
int counter = 0;
int iteration = (outSize * 10 * N )+start;

/////////////
//functions//
/////////////

//function for applying periodic boundary condition
int add1(int target){
	if(target<(L-1)){return target + 1;}
	else{return 0;}
}

int sub1(int target){
	if(target>0){return target - 1;}
	else{return L-1;}
}

//function for plotting the ising configuration (for testing purposes)
void printDisplay(){
	for ( const auto x : spin){
		if(x==1){
			cout << "\033[1;32m"<< " " << x << "\033[0m";
		}
		else{
			cout << "\33[;31m" << x << "\033[0m";
		}
	};
	cout << endl;
}

//functions for calculating the observables
//energy of the entire configuration
double totE(){
	double energ = 0;
	for(int i = 0; i<L; ++i){
		energ += spin[i]*spin[add1(i)];
	}
	energ *= -J;
	return energ;
}
//energy of the local spin flip
double localE(int xs){
	double dif = 0;
	dif += spin[xs]*spin[add1(xs)];
	dif += spin[xs]*spin[sub1(xs)];
	dif = -dif * J;
	return dif;
}
//absolute magnetization of the entire configuration
int totM(){
	int m = 0;
	for(int i = 0; i<L; ++i){
		m += spin[i];
	}
	return abs(m);
}

//spin-spin correlation function of the entire configuration
double totG(int distance){
	double g = 0;
	for(int i = 0; i<L; ++i){
		int iNext = i;
		for(int j = 0; j<distance; ++j){
			iNext = add1(iNext);
		}
		g+= spin[i] * spin[iNext];
	}
	return g/L;
}


//////////////////
//main algorithm//
//////////////////

int main(){
	clock_t begin = clock();
	cout << L << endl;
	srand(time(NULL));
	//create a folder under linux
	string folder = "configurations/";
	mkdir("configurations", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);	
	//initiate observables files 
	ofstream mcConfFile;
	ofstream mrConfFile;
	ofstream mcFile;
	string name = string("z") + to_string(steps) + string("Mc") + string("L") + to_string(L) + string(".csv");
	mcFile.open(folder+name);
	//randomize spins, plotting and calculating observables
	for(int i=0; i<L; ++i){spin[i]=((round((float)rand()/((float)RAND_MAX)))-0.5)*2;}
	printDisplay();
	E = totE();
	M = totM();
	//loop over all inverse temperatures of interest
	for(auto& B : Beta){
		B = B*factor;
		//initiate configuration files
		string name = string("z") + to_string(steps) + string("Mc") + to_string(int(round(B*10))) + string("L") + to_string(L) + string(".csv");
		mcConfFile.open(folder + name);
		name = string("z") + to_string(steps) + string("Mr") + to_string(int(round(B*10))) + string("L") + to_string(L) + string(".csv");
		mrConfFile.open(folder + name);
		//counter for surveil the configuration size
		counter = 0;
		//for loop over the iteration amount
		for(int i = 0; i < iteration; ++i){		
			//random spin position
			x=rand()%L;
			//calculate Energy difference and flip spin
			beforeE = localE(x);
			spin[x] *= -1;
			afterE = localE(x);
			difE = afterE - beforeE;
			//check spin acceptance condition
			r = ((double) rand() / (RAND_MAX));
			condition = exp(-difE *B);
			if(difE < 0||r <= condition){
				E += difE;
			}
			else{
				spin[x] *= -1;
			}
			//store observables
			if(i>=start){
				savedE += E;
				savedM += totM();
				savedG += totG(distG);
			}
			//save configurations
			//if(i>= start && i%(10*N) ==0 && counter < outSize){
			if(i>= start && i%(10*N) ==0){
				counter += 1;
				for(int k = 0; k<L; ++k){
					if(!(k == L-1)){
						mcConfFile << (spin[k]+1)/2 << ",";
						if(k%2==0){
							mrConfFile << (spin[k]+1)/2 << ",";
						}
						else{
							mrConfFile << (spin[k-1]+1)/2 << ",";
						}
					}
					else{
						mcConfFile << (spin[k]+1)/2 << endl;
						mrConfFile << (spin[k-1]+1)/2 << endl;
					}	
				}
			}
		}
		//average stored observables
		savedE /= iteration-start;
		savedM /= iteration-start;
		savedG /= iteration-start;
		//saving observables
		mcFile << B << "," << savedM << "," << savedE << "," << savedG << endl;	
		cout << "B : " << B << endl << "E : " << savedE << endl << "M : " << savedM << endl;
		printDisplay();
		//closing configuration files
		mcConfFile.close();
		mrConfFile.close();
	}
	cout << "Counter has reached: " << counter << ". It was expected to reach: " << outSize << "." << endl;
	//closing observable file
	mcFile.close();
	//end timer
	clock_t end = clock();
	double elapsed_secs = (double)(end - begin) / CLOCKS_PER_SEC;
	cout << "The exeuction of this code took: " << elapsed_secs << " seconds." << endl;
	return 0;
}
