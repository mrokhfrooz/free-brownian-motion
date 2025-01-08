#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <ctime>
#pragma warning(disable : 4996) // https://stackoverflow.com/questions/13550864/error-c4996-ctime-this-function-or-variable-may-be-unsafe
using namespace std;


double dotProduct(vector<double> a, vector<double> b) {  // input must be 2 3-elements vectors  
	double result_dot = 0;
	result_dot = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
	return result_dot;
}

vector<double> crossProduct(vector<double> a, vector<double> b) { 
	vector<double> result(3);

	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];

	return result;
}




int main() {

	cout << "############################# Background information #############################" << endl;
	cout << endl;
	cout << "According to Stokes-Einstein Equation, diffusivity of a spherical nanoparticle" << endl;
	cout << "in a solvent ,like water, can be estimated using this formula:" << endl;
	cout << endl;
	cout << "			D= (kT)/(3*pi*mu*d)" << endl;
	cout << "where" << endl;	
	cout << "k= The Boltzmann constant= 1.3806e-23[J/K]" << endl;
	cout << "T= Tempreature (Kelvin)" << endl;
	cout << "pi= pi number 3.1415" << endl;
	cout << "mu= Viscosity (Pa.S)" << endl;
	cout << "d= NP diameter (nm)" << endl;
	
	cout << endl;
	cout << "More information about this formula can be found here:" << endl;
	cout << "https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory) " << endl;
	cout << endl;
	cout << "############################# Brownian Dynamics simulation #############################" << endl;
	cout << endl;
	cout << "Here, I perform Brownian dynamics simulation to compute nanoparticle diffusivity from their trajectory." << endl;
	cout << "More information about Brownian dynamics can be found here:" << endl;
	cout << "https://en.wikipedia.org/wiki/Brownian_dynamics" << endl;
	cout << endl;
	cout << "############################# Start of simulation #############################" << endl;
	const int N = 500;//1000; // 674; // 246*5;// Number of particles
	const int step = 5000; //500000;// 100000; //100000; // 00; // 1000; // 100000; // number of steps
	const int Saveeverystep = 10; 

	int SaveN = 0;
	for (int kk = 0; kk < step + 1; kk++) {
		if (fmod(kk, Saveeverystep) == 0) {
			SaveN = SaveN + 1;
		}
	}
	//cout << "SaveN: # of rows of saved data : " << SaveN << endl;

	 srand(time(nullptr));  // Seed the random number generator with the current time random num generation for initial position // Attention
	auto start = chrono::high_resolution_clock::now(); // chat gpb this part

	ofstream myfile125("zamanha-console.txt"); // print in file all times


	


	const double pi = 3.14159265358979323846;

	double k_B = 1.3806e-23; //1.380658e-23; // unit J / K
	double T ;
	cout << "Please enter Tempreature in Kelvin (for example, 298): " ;
	cin >> T;   //310.15;  //Temp [kelvin]

	double b = 200e-9; //140 nm mesh size[meter]
	//cout << k_B * T << endl;


	double mu ;
	cout << "Please enter Water viscosity in Pa.Sec (for example, water viscosity is 0.001): " ;
	cin >> mu;   //0.001; // Pa·s% viscosity

	double p ;
	cout << "Please enter Nanoparticle diameter [nm] (for example, 200): " ;
	cin >> p;   //particle diameter

		// get the current time
	time_t currentTime = time(0);
	// convert the current time to a string
	char* timeString = ctime(&currentTime);
	// print the current time
	cout  << endl;
	cout  << endl;
	cout << "start of simulation time is: " << timeString << endl;
	myfile125 << "start of simulation time is: " << timeString << endl;
		
	 p = p*1e-9; // particle diameter
	double radius= p/2.0; // particle radius

	double mu_not=1e-6; // The rescaled timestep, hansing paper, unitless time step
double a = 0 * 1e-9;// mucin diameter[meter]
double s=p+a; // 


double dt=(mu_not*b*b*3*pi*mu*p)/(k_B*T); // time step [sec]

	double D0 = (k_B*T)/(3*pi*p*mu) ;// stokes einsitien diffusiity
	// cout << "Diffusivity from the Stokes-Einstein relation is: " << D0 << "[m^2/s]" << endl;
double eta=3*pi*mu*p; // stokes drag coefficient 
double jazr= sqrt(2*k_B*T*eta/dt); // browni force 





	double zamanscale1 = b*b/D0; // equation 5 in hansing 
	double endzaman = step * dt;  // not exported

//	cout << " *** " << endl;
//	cout << "dt= " << dt << endl;
//	cout << "endzaman " << endzaman << endl;
//	cout << " *** " << endl;


	vector <double>  zaman(step + 1);
	vector <double>  zamanD(step + 1); // D means dimensionless

	vector <double>  zaman_internal(SaveN);  // Number of datapoints for internal sampling can be different
	vector <double>  zamanD_internal(SaveN); // dimensionless times for internal sampling

	zaman[0] = 0;

	for (int i = 1; i < step + 1; i++) {
		zaman[i] = zaman[i - 1] + dt;
		//cout << zaman[i] << endl;
	}

	for (int i = 1; i < SaveN; i++) {
		zaman_internal[i] = zaman_internal[i - 1] + Saveeverystep * dt;
		zamanD_internal[i] = zaman_internal[i] / zamanscale1;
	}

	//for (int i = 0; i < step+1; i++) {
	//	cout << zaman[i] << endl;
	//}


	for (int i = 0; i < step + 1; i++) {
		zamanD[i] = zaman[i] / zamanscale1;
		//cout << zamanD[i] << endl;
	}

	vector <int> counter(step + 1);
	vector <double>  msd(step + 1);
	vector <double>  msd_internal(step + 1);
	vector <double>  msd_theory(step + 1);

	vector <double>  slope(step + 1);
	vector <double>  slope_internal(step + 1);

	vector <double>   Dtransratio(step + 1);
	vector <double>   Dtransratio_internal(step + 1);

	vector <double>   error(step + 1); // error means (D-D_SE)/ D_SE : deviations from stokes-einstin diffusivity 
	for (int i = 0; i <= step; i++) {
		counter[i] = 0;
		msd_theory[i] = 0;
		msd[i] = 0;
		msd_internal[i] = 0;
		error[i] = 0;
		slope[i] = 0;
		slope_internal[i] = 0;
		Dtransratio[i] = 0;
		Dtransratio_internal[i] = 0;
	}

	double initial_cord = b/2; //10 * 1e-9; // initial location of all particles 
	// cout << initial_cord*1e9 <<endl;


	//cout << endzaman << endl;


	static  double qx[2][N];
	static	double qx_total[2][N];
	static	double qy[2][N];
	static	double qy_total[2][N];
	static	double qz[2][N];
	static	double qz_total[2][N];
	static	double qx0[1][N];
	static	double qy0[1][N];
	static	double qz0[1][N];

	vector <vector<double >> save_qx(SaveN, vector<double>(N));
	vector <vector<double >> save_qx_total(SaveN, vector<double>(N));
	vector <vector<double >> save_qy(SaveN, vector<double>(N));
	vector <vector<double >> save_qy_total(SaveN, vector<double>(N));
	vector <vector<double >> save_qz(SaveN, vector<double>(N));
	vector <vector<double >> save_qz_total(SaveN, vector<double>(N));



	for (int i = 0; i < 2; i++) {

		for (int j = 0; j < N; j++) {
			qx[i][j] = initial_cord;
			qy[i][j] = initial_cord;
			qz[i][j] = initial_cord;
		}
	}

	//????????????????????????????????????????????????????????????????????? While  when


	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "Initialization succusfully done time is: " << timeString << endl;

	myfile125 << "Initialization succusfully done time is: " << timeString << endl;



	/////////////////////////////////////////////////////// fixing initial positions 



	for (int i = 0; i < 2; i++) {

		for (int j = 0; j < N; j++) {
			qx_total[i][j] = qx[i][j];
			qy_total[i][j] = qy[i][j];
			qz_total[i][j] = qz[i][j];
		}
	}
	for (int i = 0; i < 1; i++) {

		for (int j = 0; j < N; j++) {
			qx0[i][j] = qx[i][j];
			qy0[i][j] = qy[i][j];
			qz0[i][j] = qz[i][j];
			save_qx[i][j] = qx[i][j];
			save_qy[i][j] = qy[i][j];
			save_qz[i][j] = qz[i][j];
			save_qx_total[i][j] = qx[i][j];
			save_qy_total[i][j] = qy[i][j];
			save_qz_total[i][j] = qz[i][j];
		}
	}





	double progress[20]; // To monitor the progress of the code easily 
	progress[0] = round(0.05 * step);
	for (int i = 1; i < 20; i++) {
		progress[i] = round(progress[i - 1] + (0.05 * step));
	}


   //.................................................... Print  all constants

	ofstream ourfile("parameters.txt");

	ourfile << "Number of particles" << ";" << N << endl;
	ourfile << "step" << ";" << step << endl;
	ourfile << "Saveeverystep" << ";" << Saveeverystep << endl;
	ourfile << "#ofSavedDataPoints" << ";" << SaveN << endl;
	ourfile << "b-mesh-size [m]" << ";" << b << endl;
	ourfile << "k_B [J/K]" << ";" << k_B << endl;
	ourfile << "T [K]" << ";" << T << endl;
	ourfile << "mu [Pa.s]" << ";" << mu << endl;
    ourfile << "mu_not [unitless]" << ";" << mu_not << endl;
	ourfile << "p-partcile-diameter [m]" << ";" << p << endl;
	ourfile << "a-mucin-diameter [m]" << ";" << a << endl;
	ourfile << "dt_selected [s]" << ";" << dt << endl;
	ourfile << "endtime [s]" << ";" << endzaman << endl;
	ourfile << "D0_stokes einstein  diffusiity [m^2/s]" << ";" << D0 << endl;
   ourfile << endl;

    
    
    ourfile << "*** Pi numbers, dimensionless numbers***" << ";" << endl;
    ourfile << "mu_not:dimensionless time scale " << ";" << mu_not << endl;
	ourfile.close();

	////////////////////////////////////////////////////////////////////// print omega Beofre  main loop 

	// print qx
	ofstream yourfile4("before_qx.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			yourfile4 << qx[k][j] << ";";


		}
		yourfile4 << endl;
	}
	yourfile4.close();

	// print qy
	ofstream yourfile5("before_qy.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			yourfile5 << qy[k][j] << ";";


		}
		yourfile5 << endl;
	}
	yourfile5.close();

	// print qz
	ofstream yourfile6("before_qz.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			yourfile6 << qz[k][j] << ";";


		}
		yourfile6 << endl;
	}
	yourfile6.close();




	ofstream yourfile13("before_qx-total.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			yourfile13 << qx_total[k][j] << ";";


		}
		yourfile13 << endl;
	}
	yourfile13.close();

	ofstream yourfile14("before_qy-total.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			yourfile14 << qy_total[k][j] << ";";


		}
		yourfile14 << endl;
	}
	yourfile14.close();

	ofstream yourfile15("before_qz-total.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			yourfile15 << qz_total[k][j] << ";";


		}
		yourfile15 << endl;
	}
	yourfile15.close();

	ofstream yourfile16("before_qx0.txt");
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {

			yourfile16 << qx0[k][j] << ";";


		}
		yourfile16 << endl;
	}
	yourfile16.close();

	ofstream yourfile17("before_qy0.txt");
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {

			yourfile17 << qy0[k][j] << ";";


		}
		yourfile17 << endl;
	}
	yourfile17.close();

	ofstream yourfile18("before_qz0.txt");
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {

			yourfile18 << qz0[k][j] << ";";


		}
		yourfile18 << endl;
	}
	yourfile18.close();




	ofstream yourfile20("before_all_coordinates.txt");
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile20 << qx0[k][j] << ";";
		}
		yourfile20 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile20 << qy0[k][j] << ";";
		}
		yourfile20 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile20 << qz0[k][j] << ";";
		}
		yourfile20 << endl;
	}
	yourfile20.close();






	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





	default_random_engine generator;
	normal_distribution<double> qq1(0, 1); // source: https://cplusplus.com/reference/random/normal_distribution/
	normal_distribution<double> qq2(0, 1); // https://cplusplus.com/reference/random/normal_distribution/
	normal_distribution<double> qq3(0, 1); // https://cplusplus.com/reference/random/normal_distribution/



	int i = 0;   // do not change i anymore, i=0 because c++ basis is 0 not 1



	ofstream myfile121("progress.txt"); // to print progress in txt file
	ofstream myfile900("0save_qx.txt");
	ofstream myfile901("0save_qy.txt");
	ofstream myfile902("0save_qz.txt");

	ofstream myfile910("0save_qx_total.txt");
	ofstream myfile911("0save_qy_total.txt");
	ofstream myfile912("0save_qz_total.txt");


	ofstream myfile("Results.txt");

	myfile << "zaman[k]" << ";" << "zamanD" << ";" << "msd_theory[k]" << ";" << "msd[k]"  << ";" << "slope[k]" << ";" << "error% " << "; " << "Dtransratio " << "; " << endl;

int k_save = 0;
int shomarande=0;

// Steric force parameters

double delta=1e-10;  // small number for Differentiation
double FS1x,FS2x,FS3x,FS4x,FS5x,FS6x,FS7x,FS8x,FS9x,FS10x,FS11x,FS12x;   // Steric force - x component 
double FS1y,FS2y,FS3y,FS4y,FS5y,FS6y,FS7y,FS8y,FS9y,FS10y,FS11y,FS12y;  // Steric force - y component 
double FS1z,FS2z,FS3z,FS4z,FS5z,FS6z,FS7z,FS8z,FS9z,FS10z,FS11z,FS12z;  // Steric force - z component 
double FSx, FSy, FSz;  // total steric forces from all Rods




vector <double> v1(3); // starting point of a rod
vector <double> v2(3); // ending point of a rod 
vector <double> R_i (3);  // store positon of each particle [x,y,z]
vector <double> R_ix (3);  // used for derivation 
vector <double> R_iy (3); // used for derivation 
vector <double> R_iz (3); // used for derivation 
double rand_gauss1; // random number for brownian force
double rand_gauss2; // random number for brownian force
double rand_gauss3; // random number for brownian force
double dis, disx,disy,disz; // distance of particle from a rod 





	
	///////////////\//////////////////////////////////////////////////////////////// heart of the code
	
	
	
	
	
	
	
	
	
	
	
	
	
for (int k = 0; k < step; k++) {

		for (int j = 0; j < N; j++) {
			R_i = { qx[i][j] , qy[i][j], qz[i][j] };
			R_ix = { qx[i][j]+delta , qy[i][j], qz[i][j] };
			R_iy = { qx[i][j] , qy[i][j]+delta, qz[i][j] };
			R_iz = { qx[i][j] , qy[i][j], qz[i][j]+delta };

     ///////////////////////////////////////////////////////////////////////////////////  



////////////////////////////////////////////////////////////////////////////// Elec and Steric force calculation done
			 rand_gauss1 = qq1(generator);
			 rand_gauss2 = qq1(generator);
			 rand_gauss3 = qq2(generator);
			

			qx[i + 1][j] = qx[i][j] + (rand_gauss1*jazr*dt/eta) ; // integration
			qy[i + 1][j] = qy[i][j] + (rand_gauss2*jazr*dt/eta) ;
			qz[i + 1][j] = qz[i][j] + (rand_gauss3*jazr*dt/eta) ;




			// Peirodic boundary condition
			// x-axis
			if (qx[i + 1][j] > b)
			{
				qx[i + 1][j] = qx[i + 1][j] - b;
				qx_total[i + 1][j] = qx_total[i][j] + (qx[i + 1][j] - qx[i][j] + b);
				counter[k] = counter[k] + 1;
				
			}
			else if (qx[i + 1][j] < 0) {
				qx[i + 1][j] = qx[i + 1][j] + b;
				qx_total[i + 1][j] = qx_total[i][j] + (qx[i + 1][j] - qx[i][j] - b);
				counter[k] = counter[k] + 1;
				
			}
			else {
				qx_total[i + 1][j] = qx_total[i][j] + (qx[i + 1][j] - qx[i][j]);
			}


			// y-axis
			if (qy[i + 1][j] > b)
			{
				qy[i + 1][j] = qy[i + 1][j] - b;
				qy_total[i + 1][j] = qy_total[i][j] + (qy[i + 1][j] - qy[i][j] + b);
				counter[k] = counter[k] + 1;
				
			}
			else if (qy[i + 1][j] < 0) {
				qy[i + 1][j] = qy[i + 1][j] + b;
				qy_total[i + 1][j] = qy_total[i][j] + (qy[i + 1][j] - qy[i][j] - b);
				counter[k] = counter[k] + 1;
				
			}
			else {
				qy_total[i + 1][j] = qy_total[i][j] + (qy[i + 1][j] - qy[i][j]);
			}

			// z-axis
			if (qz[i + 1][j] > b)
			{
				qz[i + 1][j] = qz[i + 1][j] - b;
				qz_total[i + 1][j] = qz_total[i][j] + (qz[i + 1][j] - qz[i][j] + b);
				counter[k] = counter[k] + 1;
				
			}
			else if (qz[i + 1][j] < 0) {
				qz[i + 1][j] = qz[i + 1][j] + b;
				qz_total[i + 1][j] = qz_total[i][j] + (qz[i + 1][j] - qz[i][j] - b);
				counter[k] = counter[k] + 1;
				
			}
			else {
				qz_total[i + 1][j] = qz_total[i][j] + (qz[i + 1][j] - qz[i][j]);
			}



			// End of periodic BC

			msd[k + 1] = msd[k + 1] + pow(qx_total[i + 1][j] - qx0[0][j], 2) + pow(qy_total[i + 1][j] - qy0[0][j], 2) + pow(qz_total[i + 1][j] - qz0[0][j], 2);

			// cout << msd[k + 1] << endl;
			qx[i][j] = qx[i + 1][j];
			qx_total[i][j] = qx_total[i + 1][j];
		

			qy[i][j] = qy[i + 1][j];
			qy_total[i][j] = qy_total[i + 1][j];


			qz[i][j] = qz[i + 1][j];
			qz_total[i][j] = qz_total[i + 1][j];


		
			if (fmod(k + 1, Saveeverystep) == 0) {
			
				save_qx[k_save + 1][j] = qx[i][j];
				save_qy[k_save + 1][j] = qy[i][j];
				save_qz[k_save + 1][j] = qz[i][j];

				save_qx_total[k_save + 1][j] = qx_total[i][j];
				save_qy_total[k_save + 1][j] = qy_total[i][j];
				save_qz_total[k_save + 1][j] = qz_total[i][j];


				// Aim : save these parameters inside the loop
				//for (int j = 0; j < N; j++) {
				// Time = 0 does not saved now, not important much at all.
					myfile900 << save_qx[k_save + 1][j] << ";";
					myfile901 << save_qy[k_save + 1][j] << ";";
					myfile902 << save_qz[k_save + 1][j] << ";";

					myfile910 << save_qx_total[k_save + 1][j] << ";";
					myfile911 << save_qy_total[k_save + 1][j] << ";";
					myfile912 << save_qz_total[k_save + 1][j] << ";";


					//}

				if (j == N - 1) {  
					k_save++;
				}
			}



			if (k == progress[shomarande])
			{
				cout << "progress is: %" << progress[shomarande] * 100 / step << endl;

				myfile121 << "progress is: % "  << progress[shomarande] * 100 / step << ";" <<"step is " << k << " ; "<< endl;

				shomarande = shomarande + 1;
			}





		}            // End to loop N particles


		if (fmod(k + 1, Saveeverystep) == 0) {
			myfile900 << endl;
			myfile901 << endl;
			myfile902 << endl;

			myfile910 << endl;
			myfile911 << endl;
			myfile912 << endl;


		}

		// results are calculated inside the loop
		// only draw back- last data point will be losed, not important much , 1 data point
		msd[k] = msd[k] / N;
		msd_theory[k] = 6 * D0 * zaman[k];
		//cout << msd_theory [k+1]<< endl;
		slope[k] = msd[k] / (6 * zaman[k]);

		Dtransratio[k] = slope[k] / D0;
		error[k] = ((slope[k] - D0) / D0) * 100; // update kardam

		// print above results in txt file:


		myfile << zaman[k] << ";" << zamanD[k] << ";" << msd_theory[k] << ";" << msd[k]  << ";" << slope[k] << ";" << error[k] << ";" << Dtransratio[k] << ";";
		myfile << endl;




}
///////////////\//////////////////////////////////////////////////////////////// heart of the code is finished 
cout << "*******************" << endl;

	
	myfile121.close();

	myfile900.close();
	myfile901.close();
	myfile902.close();

	myfile910.close();
	myfile911.close();
	myfile912.close();




	myfile.close();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// End ot main asli

	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "End of main loop : " << timeString << endl;
	myfile125 << "End of main loop : " << timeString << endl;




	// print qx_total
	ofstream myfile5("qx_total.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			myfile5 << qx_total[k][j] << ";";


		}
		myfile5 << endl;
	}
	myfile5.close();

	string bob;

	////////////////////////////////////////////
	// print qy_total
	ofstream myfile6("qy_total.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			myfile6 << qy_total[k][j] << ";";


		}
		myfile6 << endl;
	}
	myfile6.close();

	ofstream myfile7("qz_total.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			myfile7 << qz_total[k][j] << ";";


		}
		myfile7 << endl;
	}
	myfile7.close();

	// print qx
	ofstream myfile8("qx.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			myfile8 << qx[k][j] << ";";


		}
		myfile8 << endl;
	}
	myfile8.close();

	// print qy
	ofstream myfile9("qy.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			myfile9 << qy[k][j] << ";";


		}
		myfile9 << endl;
	}
	myfile9.close();

	// print qz
	ofstream myfile10("qz.txt");
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < N; j++) {

			myfile10 << qz[k][j] << ";";


		}
		myfile10 << endl;
	}
	myfile10.close();


	ofstream myfile11("counter.txt");
	for (int k = 0; k < step + 1; k++) {


		myfile11 << counter[k] << ";";


		myfile11 << endl;
	}
	myfile11.close();


	ofstream yourfile21("after_all_coordinates.txt");
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile21 << qx[k][j] << ";";
		}
		yourfile21 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile21 << qy[k][j] << ";";
		}
		yourfile21 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile21 << qz[k][j] << ";";
		}
		yourfile21 << endl;
	}
	yourfile21.close();


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Internal sampling
		// internal sampling
	// Print current time Start
// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "start of internal sampling  time is: " << timeString << endl;
	// Print current time  End of block
	myfile125 << "start of internal sampling  time is: " << timeString << endl;


	int numberOfdeltaT = save_qx_total.size();
//	cout << "numberOfdeltaT " << numberOfdeltaT << endl;
	int left = numberOfdeltaT - floor(numberOfdeltaT / 2);

//	cout << "left " << left << endl;
	// cout << "numberOfdeltaT:  " << numberOfdeltaT << endl;
	for (int dt = 1; dt < numberOfdeltaT; dt++) {  // dt=0 cherte mahz
		double msd_dt = 0.0;
		double mdr_dt = 0.0;
		int numElements = 0;

		if (dt < (numberOfdeltaT / 2)) {
			//numElements = left - 1;
			//cout << " numElements " << numElements << endl;

			for (int i = dt; i < left + dt; i++) {  // shak
				numElements = numElements + 1;
				for (int j = 0; j < N; j++) {
					//    cout << "save_qx_total[i - dt][j] " << save_qx_total[i - dt][j] * 1e9 << endl;
					//	cout << "save_qx_total[i][j] " << save_qx_total[i][j] *1e9 << endl;

					msd_dt = msd_dt + (pow(save_qx_total[i][j] - save_qx_total[i - dt][j], 2) + pow(save_qy_total[i][j] - save_qy_total[i - dt][j], 2) + pow(save_qz_total[i][j] - save_qz_total[i - dt][j], 2));
					//cout << "msd_dt " << msd_dt << endl;
				}
				//cout << " *** " << endl;
			}
		}
		/*if (dt < (numberOfdeltaT / 2)) {  // idea behin this block was wrong, resulted in merge Hansing1000 and OS1000 at the middle point
			numElements = numberOfdeltaT - dt-left+1;  // shak   // +1 at the end of line: found and corrected by printing the qx and numElements when N=1 and step=10;
		   // cout << " numElements " << numElements << endl;
		   for (int i = dt; i < left; i++) {  // shak
			   for (int j = 0; j < N; j++) {

				   mdr_dt = mdr_dt + ((save_omega[3 * i][j] * save_omega[3 * (i - dt)][j]) + (save_omega[3 * i + 1][j] * save_omega[3 * (i - dt) + 1][j]) + (save_omega[3 * i + 2][j] * save_omega[3 * (i - dt) + 2][j]));
				   msd_dt = msd_dt + (pow(save_qx_total[i][j] - save_qx_total[i - dt][j], 2) + pow(save_qy_total[i][j] - save_qy_total[i - dt][j], 2) + pow(save_qz_total[i][j] - save_qz_total[i - dt][j], 2));
				   //cout << "msd_dt " << msd_dt << endl;
			   }
			   //	cout << " *** " << endl;
		   }
	   } */
		else {
			numElements = numberOfdeltaT - dt;
			for (int i = dt; i < numberOfdeltaT; i++) {
				for (int j = 0; j < N; j++) {
				

					msd_dt = msd_dt + (pow(save_qx_total[i][j] - save_qx_total[i - dt][j], 2) + pow(save_qy_total[i][j] - save_qy_total[i - dt][j], 2) + pow(save_qz_total[i][j] - save_qz_total[i - dt][j], 2));
					//cout << "msd_dt " << msd_dt << endl;
				}
				//	cout << " *** " << endl;
			}
		}
		//	cout << ">>>>>>" << endl;
		msd_internal[dt] = msd_dt / numElements;

	}



	// Calculate results for internal sampling
	for (int k = 1; k < SaveN; k++) {

		msd_internal[k] = msd_internal[k] / N;

		slope_internal[k] = msd_internal[k] / (6 * zaman_internal[k]);


		Dtransratio_internal[k] = slope_internal[k] / D0;
	}

	ofstream myfile99("Internal_TranslationRotation.txt");
	myfile99 << "zamanD_internal" << ";" << " Dtransratio_internal" << "; " << endl;
	for (int k = 0; k < SaveN; k++) {


		myfile99 << zamanD_internal[k] << ";"  << Dtransratio_internal[k] << ";";


		myfile99 << endl;
	}
	myfile99.close();

	// Print current time
// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time

	// Print current time

	cout << "End of internal sampling time is: " << timeString << endl;
	myfile125 << "End of internal sampling time is: " << timeString << endl;
	///////////////////////////////////////////////////////////////////////////////??????????????????????????????????????????? internal sampling finished SAVE ha

		// print qX-save
	ofstream myfile90("save_qx.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile90 << save_qx[k][j] << ";";
		}
		myfile90 << endl;
	}
	myfile90.close();

	// print qY-save
	ofstream myfile91("save_qy.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile91 << save_qy[k][j] << ";";
		}
		myfile91 << endl;
	}
	myfile91.close();

	// print qZ-save
	ofstream myfile92("save_qz.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile92 << save_qz[k][j] << ";";
		}
		myfile92 << endl;
	}
	myfile92.close();

	// print save_qx_total
	ofstream myfile93("save_qx_total.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile93 << save_qx_total[k][j] << ";";
		}
		myfile93 << endl;
	}
	myfile93.close();

	// print save_qy_total
	ofstream myfile94("save_qy_total.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile94 << save_qy_total[k][j] << ";";
		}
		myfile94 << endl;
	}
	myfile94.close();

	// print save_qz_total
	ofstream myfile95("save_qz_total.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile95 << save_qz_total[k][j] << ";";
		}
		myfile95 << endl;
	}
	myfile95.close();


	///////////////////////////////////////////////////////////////////////////////???????????????????????????????????????????  SAVE ha

	// Print current time
	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "The End of simulation  time is: " << timeString << endl;
	myfile125 << "The End of simulation  time is: " << timeString << endl;
	// Print current time


 	double slope_avg_result;
 	for (int uu =1 ; uu< step ; uu++)
 	{
 		slope_avg_result= slope_avg_result+slope[uu]; // Do summation
	 }
	 slope_avg_result=slope_avg_result/step ;  // Do averaging
 	
	auto end = chrono::high_resolution_clock::now();  // // chat gpb this part

	auto duration = chrono::duration_cast<chrono::seconds>(end - start);

	cout << "Execution time of simulation: " << duration.count() << " seconds\n";

	myfile125.close();
	
	cout << "Diffusivity from the Stokes-Einstein relation is: " << D0 << "[m^2/s]" << endl;
	cout << "Note: This result can be checked with this website: https://www.vcalc.com/wiki/stokes-einstein-diffusion-coefficient" << endl;
	cout << "Diffusivity of our Brownian dynamics simulation is: " << slope_avg_result << "[m^2/s]" << endl;
	cout << "Error between two results is: " << (slope_avg_result-D0 )*100/(D0) << endl;

	return 0;

}
// Changes:

// zaman= time, I did not use the "time" keyword as some C++ libraries  use this keyword too
