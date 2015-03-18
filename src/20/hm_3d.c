/*
Project: 3D Heisenberg Model MCMA Simulation
         Note: mu=b=B=k_b=1 (in magnitude only, of course)
Author: Michael Conroy
Class: PHY 471
Professor: Dr. Enjalran
Date: April 2014
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/******************************************************
These are the paramters that are changed to run an MCMA
simulation of the 3D Heisenberg Model. No changes are
necessary in the source code itself.
******************************************************/
#define LENGTH 20 //Length of lattice (as doible for angle types)
#define N 8000 //LENGTH^3; number of spins/sites
#define J 1.0 //Exhange constant (as double for arithmetic)
#define MIN_TEMP 0.1 //Minimum temperature in Kelvin
#define MAX_TEMP 5.0 //Max temperature in Kelvin
#define MAX_TEMP_STEPS 10 //Max number divisions for temps (0,MAX_TEMP] (as a double for arithmetic)
#define Y_SIZE 2 ////For 3d array declaration; index corresponds to E, E^2 (3D), or theta and phi (4D)
#define M_SIZE 4 //For 3d array declaration for magnetization
#define PI M_PI //The number pi set to constant in math.h
#define RUNS 250 //Num equilibration and measurement phase iterations
#define SKIPS 30 //Num skipps between measurements
#define SWEEPS 250 //Num sweeps between each skip
#define SEED 12345678 //Seed parameters
#define PRINT_FILE 1 //Whether or not to print to file
#define PRINT_TO_SCREEN 0 //Whether or not print data to screen
#define DATA_DEBUG 0 //Whether or not to output data for debugging
#define WHERE_IS 1 //Where the program is while runnin

double rng(void);
double rng_theta(void);
double rng_phi(void);
void temp_init(double temps[]);
void initialize(double ****lattice);
double array_sum_1(double ***a,int t);
double array_sum_2(double ***a,int t);
double ecalc(double ****lattice);
double mcalc(double ****lattice);
void equilibration(double ****lattice,double temperature);
void measurement(double ****lattice, double data[][13], double temperature, double ***energy, double ***magnetization, int t);
//void data_calc(double ***energy, double ***magnetization,double data[][13],double temperature,int t)
void print_to_file(double data[][13]);

int main(void)
{
    printf("\n\n");
    printf("Project: Heisenberg Model\n");
    printf("Author: Michael Conroy\n");
    printf("Class: PHY 471\n");
    printf("Date: February 2014\n\n");
    printf("Simulation Parameters:\n");
    printf("Lattice Dimensions: %d\n",LENGTH);
    printf("Number of Sites: %d\n",N);
    printf("Min Temp (K): %1.1f\n",MIN_TEMP);
    printf("Max Temp (K): %1.1f\n",MAX_TEMP);
    printf("Runs: %d\n",RUNS);
    printf("Skips (N): %d\n",SKIPS);
    printf("Sweeps: %d\n",SWEEPS);
    printf("Print to file: ");
    printf(PRINT_FILE ? "true" : "false");
    printf("\n\n");

    int t; //Temp counter
    int i,j,k; //3D and 4D array counters
    //double lattice[LENGTH][LENGTH][LENGTH][2] = {}; //The Lattice[x loc][y loc][z loc][1=theta/2=phi]
    //double theta, phi;
    //double temperature=999.0;
    //double dT = MAX_TEMP/MAX_TEMP_STEPS; //Temp precision
    double temps[(int)MAX_TEMP_STEPS]; //Holds temperatures; temps[99] = MAX_TEMP
    double temperature;
    //double angles[1001][2]; //Holds angles; indexed from 'angles[0]' to 'angles[1000]'
    //double energy[(int)MAX_TEMP_STEPS][2][RUNS]; //Energies taken during measurement phase
    //double magnetization[(int)MAX_TEMP_STEPS][2][RUNS]; //Magnetization taken during measurement phase
    double data[(int)MAX_TEMP_STEPS][13]; //Holds calculated data

    /*--------------------------------------------*/
    /*--------------Declaring Arrays--------------*/
    /*--------------------------------------------*/
    //Declare 3D Arrays
    printf("Declaring 3D arrays...\n\n");
    //energy
    double ***energy;
    energy = (double***)malloc(MAX_TEMP_STEPS * sizeof(double **));
        if (energy == NULL)
        {
            printf("1: Out of memory!\n");
            exit(1);
        }
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        energy[i] = (double**)malloc(Y_SIZE * sizeof(double *));
            if (energy == NULL)
            {
                printf("2: Out of memory!\n");
                exit(2);
            }
        for(j = 0; j < Y_SIZE; j++)
            {
                energy[i][j] = (double*)malloc(RUNS * sizeof(double));
                    if (energy == NULL)
                    {
                        printf("3: Out of memory!\n");
                        exit(3);
                    }
            }
    }
    //magnetization
    double ***magnetization;
    magnetization = (double***)malloc(MAX_TEMP_STEPS * sizeof(double **));
        if (magnetization == NULL)
        {
            printf("4: Out of memory!\n");
            exit(4);
        }
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        magnetization[i] = (double**)malloc(M_SIZE * sizeof(double *));
            if (magnetization == NULL)
            {
                printf("5: Out of memory!\n");
                exit(5);
            }
        for(j = 0; j < Y_SIZE; j++)
            {
                magnetization[i][j] = (double*)malloc(RUNS * sizeof(double));
                    if (magnetization == NULL)
                    {
                        printf("6: Out of memory!\n");
                        exit(6);
                    }
            }
    } //Declare 3D arrays
    //Declare 4D arrays
    printf("Declaring 4D arrays...\n\n");
    //lattice
    //double lattice[LENGTH][LENGTH][LENGTH][2]; //The Lattice[x loc][y loc][z loc][0=theta/1=phi]
    double ****lattice;
    lattice = (double****)malloc(LENGTH * sizeof(double ***));
        if (lattice == NULL)
        {
            printf("7: Out of memory!\n");
            exit(7);
        }
    for(i = 0; i < LENGTH; i++)
    {
    /*--------------------------------------------*/
        lattice[i] = (double***)malloc(LENGTH * sizeof(double**));
            if (lattice == NULL)
            {
                printf("8: Out of memory!\n");
                exit(8);
            }
        for(j = 0; j < LENGTH; j++)
        {
            lattice[i][j] = (double**)malloc(LENGTH * sizeof(double*));
                if (lattice == NULL)
                {
                    printf("9: Out of memory!\n");
                    exit(9);
                }
            for (k = 0; k < LENGTH; k++)
            {
                lattice[i][j][k] = (double*)malloc(Y_SIZE * sizeof(double));
                    if (lattice == NULL)
                    {
                        printf("10: Out of memory!\n");
                        exit(10);
                    }
            }
        }
    } //Declare 4D arrays
    /*--------------------------------------------*/
    /*--------------Declaring Arrays--------------*/
    /*--------------------------------------------*/

    //Seed random number generator
    printf("Seeding RNG...\n\n");
    srand(SEED);

    //Temperature array initialization
    printf("Initializing temperature array...\n\n");
    temp_init(temps);

    //Initialize lattice spins
    printf("Initializing lattice...\n\n");
    initialize(lattice);

    //MCMA Algorithm
    printf("Entering MCMA...\n\n");
    for(t = 0; t < MAX_TEMP_STEPS; t++)
    {
        //Indicates what step MCMA is at
        if (WHERE_IS){printf("\tStep: %d   Temp: %f\n",t,temps[t]);}
        //Calc temoper
        temperature = temps[t];
        //Equilibration Phase
            if (WHERE_IS){printf("\t\tEquilibration...\n");}
        equilibration(lattice,temperature);
        //Measurement Phase
            if (WHERE_IS){printf("\t\tMeasurement...\n");}
        measurement(lattice,data,temperature,energy,magnetization,t);
        //Print to screen
        if(PRINT_TO_SCREEN)
        {
            printf("data array (main): %f\t %f\t %f\t %f\t %f\t %f\t%f\t%f\t%f\t%f\t%f\t%f\n\n",
            data[t][0],data[t][1],data[t][2],data[t][4],data[t][5],data[t][6], data[t][7],data[t][8],data[t][9],data[t][10],data[t][11],data[t][12]);
        }
    }
    if (WHERE_IS){printf("\tEnd Algorithm...\n\n");}

    //Print to file
    if (PRINT_FILE) {print_to_file(data);}

    /*--------------------------------------------*/
    /*--------------Deallocate Arrays-------------*/
    /*--------------------------------------------*/
    //Deallocate 3D and 4D Arrays
    //energy
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        for(j = 0; j < Y_SIZE; j++)
        {
            free(energy[i][j]);
        }
        free(energy[i]);
    }
    free(energy);
    //magnetization
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        for(j = 0; j < Y_SIZE; j++)
        {
            free(magnetization[i][j]);
        }
        free(magnetization[i]);
    }
    free(magnetization);
    //lattice
    for(i = 0; i < LENGTH; i++)
    {
        for(j = 0; j < LENGTH; j++)
        {
            for (k = 0; k < LENGTH; k++)
            {
                free(lattice[i][j][k]);
            }
            free(lattice[i][j]);
        }
        free(lattice[i]);
    }
    free(lattice);
    /*--------------------------------------------*/
    /*--------------Deallocate Arrays-------------*/
    /*--------------------------------------------*/
    //End program
    printf("End of program...\n\n");
    //Return
    return 0;
}

//Random number [0,1]; regular random number generation
double rng(void)
{
    double r;
    r = (double)(rand()/ (double)RAND_MAX);
    return r;
}

//Random theta [0,pi]
double rng_theta(void)
{
    double r = 1.0/cos(1.0-2.0*rng()); //Exn 16.24, p. 399 Newman and Barkema
    return r;
}

//Random phi [0,2*pi]
double rng_phi(void)
{
    double r = 2*PI*rng();
    return r;
}

//Temperature array initialization
void temp_init(double temps[])
{
    int t; //Counters
    //Temperature array initialization
    for (t = 0; t < MAX_TEMP_STEPS; t++)
    {
        double dT = MAX_TEMP/MAX_TEMP_STEPS;
        temps[t] = MAX_TEMP-t*dT;
    }
}

//Initialize the lattice/spin array
void initialize(double ****lattice)
{
    int i,j,k; //Counters
    //Populate initialization lattice
    for (i = 0 ; i < LENGTH; i++)
    {
        for(j = 0; j < LENGTH; j++)
        {
            for (k = 0; k < LENGTH; k++)
            {
                lattice[i][j][k][0] = rng_theta();
                lattice[i][j][k][1] = rng_phi();
                //printf("Init latice: %f\t%f\n",lattice[i][j][k][1],lattice[i][j][k][2]);
            }

        }
    }
}

//Sum an array with number of elements eq
double array_sum_1(double ***a,int t)
{
    int i;
    double sum = 0;
    for(i = 0; i < RUNS; i++)
    {
        sum = sum + a[t][0][i];
    }
    return(sum);
}

//Sum an array with number of elements eq
double array_sum_2(double ***a,int t)
{
    int i;
    double sum = 0;
    for(i = 0; i < RUNS; i++)
    {
        sum = sum + a[t][1][i];
    }
    return(sum);
}

/*
//Calculate energy
double ecalc(double ****lattice)
{
    int i,j,k; //Counters
    int im,ip,jm,jp,km,kp; //Lattice indices
    int sw1, sw2, sw3; //Sweep counters
    double theta, phi; //Spin angles for current spin
    double e_storage = 0; //Holds sum of each spin's energy calculation
    double spinX,spinY,spinZ; //Components for current spin
    double spin1theta, spin1phi, spin2theta, spin2phi, spin3theta, spin3phi, //Spin angles for each n.n.
           spin4theta, spin4phi, spin5theta, spin5phi, spin6theta, spin6phi;
    double spinSum[3]; //Sum of n.n. components
    double dotProduct; //Valye of dot product
    //Energy sweep - calculate energy of lattice in single sweep
    for (sw1 = 0; sw1 < LENGTH; sw1++)
    {
    for (sw2 = 0; sw2 < LENGTH; sw2++)
    {
    for (sw3 = 0; sw3 < LENGTH; sw3++)
    {
        //Select random spin's location
        i = sw1; //rng_int_L(); //selects 0,1,2,3,4 for LENGTH=5
        j = sw2; //rng_int_L();
        k = sw3; //rng_int_L();
        //Six nearest neighbors
        ip = (i+LENGTH-1) % LENGTH;
        im = (i+1) % LENGTH;
        jp = (j+LENGTH-1) % LENGTH;
        jm = (j+1) % LENGTH;
        kp = (k+LENGTH-1) % LENGTH;
        km = (k+1) % LENGTH;
        //Current angles
        theta = lattice[i][j][k][0];
        phi = lattice[i][j][k][1];
        //Dot product calculation
        //E = -J (Si_current)*(S1 + S2+ S3+ S4+ S5+ S6)
        //NOTE: (Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
        //Left side of dot product
        spinX = sin(theta)*cos(phi);
        spinY = sin(theta)*sin(phi);
        spinZ = cos(theta);
        //Right side dot product
        spin1theta = lattice[im][j][k][0];
        spin1phi   = lattice[im][j][k][1];
        spin2theta = lattice[ip][j][k][0];
        spin2phi   = lattice[ip][j][k][1];
        spin3theta = lattice[i][jm][k][0];
        spin3phi   = lattice[i][jm][k][1];
        spin4theta = lattice[i][jp][k][0];
        spin4phi   = lattice[i][jp][k][1];
        spin5theta = lattice[i][j][km][0];
        spin5phi   = lattice[i][j][km][1];
        spin6theta = lattice[i][j][kp][0];
        spin6phi   = lattice[i][j][kp][1];
        //Putting the right side of the dot product together....
        //                     Spin 1                              Spin 2                            Spin 3                             Spin 4                            Spin 5                            Spin 6
        spinSum[0] = sin(spin1theta)*cos(spin1phi) +    sin(spin2theta)*cos(spin2phi) +     sin(spin3theta)*cos(spin3phi) +   sin(spin4theta)*cos(spin4phi) +   sin(spin5theta)*cos(spin5phi) +   sin(spin6theta)*cos(spin6phi);
        spinSum[1] = sin(spin1theta)*sin(spin1phi) +    sin(spin2theta)*sin(spin2phi) +     sin(spin3theta)*sin(spin3phi) +   sin(spin4theta)*sin(spin4phi) +   sin(spin5theta)*sin(spin5phi) +   sin(spin6theta)*sin(spin6phi);
        spinSum[2] = cos(spin1theta) +                  cos(spin2theta) +                   cos(spin3theta) +                 cos(spin4theta) +                 cos(spin5theta) +                 cos(spin6theta);
        //Compute dot product
        dotProduct = spinX*spinSum[0] + spinY*spinSum[1] + spinZ*spinSum[2];
        //Energy calculation
        e_storage = e_storage + -J*dotProduct; //

    } //sw1
    } //sw2
    } //sw3
    //Return energy value
    return e_storage/2.0;
}
*/

/*
//Calculate magnetization
double Xcalc(double ****lattice)
{
    int i,j,k; //Counters
    int sw1, sw2, sw3; //Sweep counters
    double theta, phi; //Spin angles for current spin
    double Mxsum = 0.0, Mysum = 0.0, Mzsum = 0.0;
    double X,Y,Z; //Components for current spin
    double Mx,My,Mz,M,m,msqr;
    //Energy sweep - calculate energy of lattice in single sweep
    for (sw1 = 0; sw1 < LENGTH; sw1++)
    {
    for (sw2 = 0; sw2 < LENGTH; sw2++)
    {
    for (sw3 = 0; sw3 < LENGTH; sw3++)
    {
        //Select random spin's location
        i = sw1; //rng_int_L(); //selects 0,1,2,3,4 for LENGTH=5
        j = sw2; //rng_int_L();
        k = sw3; //rng_int_L();
        //Current angles
        theta = lattice[i][j][k][0];
        phi = lattice[i][j][k][1];
        //Dot product calculation
        //m_rms = M = (1/N)sqrt(sum(S_i*S_i))
        //NOTE: (Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
        //Taken from Landau and Binder p.153, (5.16)
        X = sin(theta)*cos(phi);
        Y = sin(theta)*sin(phi);
        Z = cos(theta);
        Mxsum = Mxsum + X;
        Mysum = Mysum + Y;
        Mzsum = Mzsum + Z;
    } //sw1
    } //sw2
    } //sw3
    //Calc M_x, M_y, M_z
    //Ma = (1/N)*sum_i(Sia), where a = x, y, or z
    SxAvg = (1/(double)N)*Mxsum;
    SyAvg = (1/(double)N)*Mysum;
    SzAvg = (1/(double)N)*Mzsum;
    //Chi = beta*[<sx^2> + <sy^2> + <sz^2>  -  <sx> - <sy> - <sz> ]
    //where <sx> = 1/N sum_i Sx  and <sx^2> = (1/N sum_i Sx )^2
    //Calculate non-squared magnetization
    //m = Mx + My + Mz;
    //Calculate squared magnetization
    //msqr = Mx*Mx + My*My + Mz*Mz;
    //Calculate magnetization
    X = SxAvg*SxAvg + SyAvg*SyAvg + SzAvg*SzAvg - SxAvg - SyAvg - SzAvg
    //Return magnetization value
    return X;
}
*/
//Equilibration phase
void equilibration(double ****lattice,double temperature)
{
    int eq,sw,sw1,sw2,sw3; //Equil and sweep counter
    int i,j,k,ip,im,jp,jm,kp,km; //Lattice indices for spin location
    double theta, thetaPrime, phi, phiPrime;
    double spin1theta, spin1phi, spin2theta, spin2phi, spin3theta, spin3phi,
           spin4theta, spin4phi, spin5theta, spin5phi, spin6theta, spin6phi;
    double newSpin[2], currentSpin[2]; //New and old spin values
    double spinDiff[3], spinSum[3];
    double dotProduct;
    double deltaE; //Change in energy
    double prob; //Prob of selecting spin if deltaE > 0
    double randNum; //Random number for probability
    double beta = 1.0/temperature;

    //Equilibration loop - sweeps through lattice
    for (eq = 0; eq < RUNS; eq++)
    {
        //Sweeps
        for (sw = 0; sw < SWEEPS; sw++)
        {
            for (sw1 = 0; sw1 < LENGTH; sw1++)
            {
            for (sw2 = 0; sw2 < LENGTH; sw2++)
            {
            for (sw3 = 0; sw3 < LENGTH; sw3++)
            {
                //Select site location
                //Random
                //rng_int_L(); //selects 0,1,2,3,4 for LENGTH=5
                //rng_int_L();
                //rng_int_L();
                //One-by-one; sweep through each in succession
                i = sw1;
                j = sw2;
                k = sw3;
                //Six nearest neighbors
                ip = (i+LENGTH-1) % LENGTH;
                im = (i+1) % LENGTH;
                jp = (j+LENGTH-1) % LENGTH;
                jm = (j+1) % LENGTH;
                kp = (k+LENGTH-1) % LENGTH;
                km = (k+1) % LENGTH;
                //Current spin value
                currentSpin[0] = lattice[i][j][k][0];
                currentSpin[1] = lattice[i][j][k][1];
                //Propose new spin value
                newSpin[0] = rng_theta();
                newSpin[1] = rng_phi();
                //Current angles
                theta = currentSpin[0];
                phi = currentSpin[1];

                //MC for Theta
                //(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
                //deltaE = -J (Si_new – Si_old)*(S1 + S2+ S3+ S4+ S5+ S6)
                //Test only theta, so keep phiPrime the same
                thetaPrime = newSpin[0];
                phiPrime = phi; //MC for theta only; phi doesn't change
                //Left side of dot product
                spinDiff[0] = sin(thetaPrime)*cos(phiPrime) - sin(theta)*cos(phi);
                spinDiff[1] = sin(thetaPrime)*sin(phiPrime) - sin(theta)*sin(phi);
                spinDiff[2] = cos(thetaPrime) - cos(theta);
                //Right side dot product
                spin1theta = lattice[im][j][k][0];
                spin1phi   = lattice[im][j][k][1];
                spin2theta = lattice[ip][j][k][0];
                spin2phi   = lattice[ip][j][k][1];
                spin3theta = lattice[i][jm][k][0];
                spin3phi   = lattice[i][jm][k][1];
                spin4theta = lattice[i][jp][k][0];
                spin4phi   = lattice[i][jp][k][1];
                spin5theta = lattice[i][j][km][0];
                spin5phi   = lattice[i][j][km][1];
                spin6theta = lattice[i][j][kp][0];
                spin6phi   = lattice[i][j][kp][1];
                //Putting the right side of the dot product together....
                //                     Spin 1                              Spin 2                            Spin 3                             Spin 4                            Spin 5                            Spin 6
                spinSum[0] = sin(spin1theta)*cos(spin1phi) +    sin(spin2theta)*cos(spin2phi) +     sin(spin3theta)*cos(spin3phi) +   sin(spin4theta)*cos(spin4phi) +   sin(spin5theta)*cos(spin5phi) +   sin(spin6theta)*cos(spin6phi);
                spinSum[1] = sin(spin1theta)*sin(spin1phi) +    sin(spin2theta)*sin(spin2phi) +     sin(spin3theta)*sin(spin3phi) +   sin(spin4theta)*sin(spin4phi) +   sin(spin5theta)*sin(spin5phi) +   sin(spin6theta)*sin(spin6phi);
                spinSum[2] = cos(spin1theta) +                  cos(spin2theta) +                   cos(spin3theta) +                 cos(spin4theta) +                 cos(spin5theta) +                 cos(spin6theta);
                //Compute dot product
                dotProduct = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
                //Change in energy
                deltaE = -J*dotProduct;
                //deltaE = -J*(cos(newSpin)-cos(curretSpin));
                //Acceptance check
                if (deltaE <= 0)
                {
                    lattice[i][j][k][0] = thetaPrime;
                    lattice[i][j][k][1] = phiPrime; //Phi coordinate doesn't change
                }

                else //deltaE > 0
                {
                    prob = (sin(thetaPrime)/sin(theta))*exp(-beta*deltaE);
                    randNum = rng();
                    if (randNum <= prob)
                    {
                        lattice[i][j][k][0] = thetaPrime;
                        lattice[i][j][k][1] = phiPrime; //Phi coordinate doesn't change
                    }
                    //else //randNum > probability; no change
                    //    lattice[i][j][k][0] = currentSpin[0];
                    //    lattice[i][j][k][1] = currentSpin[2];
                }

                //MC for Phi
                //Test only phi, so keep thetaPrime the same
                thetaPrime = theta;
                phiPrime = newSpin[1];
                //Left side of dot product
                spinDiff[0] = sin(thetaPrime)*cos(phiPrime) - sin(theta)*cos(phi);
                spinDiff[1] = sin(thetaPrime)*sin(phiPrime) - sin(theta)*sin(phi);
                spinDiff[2] = cos(thetaPrime) - cos(theta);
                //Right side dot product
                spin1theta = lattice[im][j][k][0];
                spin1phi   = lattice[im][j][k][1];
                spin2theta = lattice[ip][j][k][0];
                spin2phi   = lattice[ip][j][k][1];
                spin3theta = lattice[i][jm][k][0];
                spin3phi   = lattice[i][jm][k][1];
                spin4theta = lattice[i][jp][k][0];
                spin4phi   = lattice[i][jp][k][1];
                spin5theta = lattice[i][j][km][0];
                spin5phi   = lattice[i][j][km][1];
                spin6theta = lattice[i][j][kp][0];
                spin6phi   = lattice[i][j][kp][1];
                //Putting the right side of the dot product together....
                //                     Spin 1                              Spin 2                            Spin 3                             Spin 4                            Spin 5                            Spin 6
                spinSum[0] = sin(spin1theta)*cos(spin1phi) +    sin(spin2theta)*cos(spin2phi) +     sin(spin3theta)*cos(spin4phi) +   sin(spin4theta)*cos(spin4phi) +   sin(spin5theta)*cos(spin5phi) +   sin(spin6theta)*cos(spin6phi);
                spinSum[1] = sin(spin1theta)*sin(spin1phi) +    sin(spin2theta)*sin(spin2phi) +     sin(spin3theta)*sin(spin4phi) +   sin(spin4theta)*sin(spin4phi) +   sin(spin5theta)*sin(spin5phi) +   sin(spin6theta)*sin(spin6phi);
                spinSum[2] = cos(spin1theta) +                  cos(spin2theta) +                   cos(spin3theta) +                 cos(spin4theta) +                 cos(spin5theta) +                 cos(spin6theta);
                //Compute dot product
                dotProduct = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
                //Change in energy
                deltaE = -J*dotProduct;
                //deltaE = -J*(cos(newSpin)-cos(curretSpin));
                //Acceptance check
                if (deltaE <= 0)
                {
                    lattice[i][j][k][0] = thetaPrime; //Theta coordinate doesn't change
                    lattice[i][j][k][1] = phiPrime;
                }

                else //deltaE > 0
                {
                    prob = (sin(thetaPrime)/sin(theta))*exp(-beta*deltaE);
                    randNum = rng();
                    if (randNum <= prob)
                    {
                        lattice[i][j][k][0] = thetaPrime;
                        lattice[i][j][k][1] = phiPrime;
                    }
                    //else //randNum > probability; no change
                    //    lattice[i][j][k][0] = currentSpin[0];
                    //    lattice[i][j][k][1] = currentSpin[2];
                }
            } //sw3
            } //sw2
            } //sw1
        } //sw
        //Data debug
        if (DATA_DEBUG){printf("\t\tEquilibration runs (0-%d): %d\n",RUNS,eq);}
    } //eq
}

//Measurement phase
void measurement(double ****lattice, double data[][13], double temperature, double ***energy, double ***magnetization, int t)
{
    int RUNS_ = (double)RUNS;
    int N_ = (double)N;
    int eq,sk,sw,sw1,sw2,sw3; //Equil and sweep counter
    int ii; //Array summation counter
    int i,j,k,ip,im,jp,jm,kp,km; //Lattice indices for spin location
    double theta, thetaPrime, phi, phiPrime;
    double spin1theta, spin1phi, spin2theta, spin2phi, spin3theta, spin3phi,
           spin4theta, spin4phi, spin5theta, spin5phi, spin6theta, spin6phi;
    double newSpin[2], currentSpin[2]; //New and old spin values
    double spinDiff[3], spinSum[3];
    double dotProduct;
    double deltaE; //Change in energy
    double prob; //Prob of selecting spin if deltaE > 0
    double randNum; //Random number for probability
    double beta = 1.0/temperature;
    double accept, test;
    double Mxsum = 0.0, Mysum = 0.0, Mzsum = 0.0;
    double spinX,spinY,spinZ; //Components for current spin
    double X,Y,Z; //Components for current spin
    double Mx[RUNS_],My[RUNS_],Mz[RUNS_];
    double m[RUNS_], m2[RUNS_], e[RUNS_], e2[RUNS_];
    double e_storage = 0; //Holds sum of each spin's energy calculation
    double eSum = 0,e2Sum = 0, mSum = 0, m2Sum = 0;
    double SxSum = 0, SySum = 0, SzSum = 0;
    double SxAvg, SyAvg, SzAvg;
    double SxAvgSqr, SyAvgSqr, SzAvgSqr;

    //Equilibration/Measurement loop - sweeps through lattice and records data
    for (eq = 0; eq < RUNS; eq++)
    {
        //Skips
        for (sk= 1; sk <= SKIPS; sk++)
        {
            //Sweeps
            for (sw = 0; sw < SWEEPS; sw++)
            {
                for (sw1 = 0; sw1 < LENGTH; sw1++)
                {
                for (sw2 = 0; sw2 < LENGTH; sw2++)
                {
                for (sw3 = 0; sw3 < LENGTH; sw3++)
                {
                    //Select site location
                    //Random
                    //rng_int_L(); //selects 0,1,2,3,4 for LENGTH=5
                    //rng_int_L();
                    //rng_int_L();
                    //One-by-one; sweep through each in succession
                    i = sw1;
                    j = sw2;
                    k = sw3;
                    //Six nearest neighbors
                    ip = (i+LENGTH-1) % LENGTH;
                    im = (i+1) % LENGTH;
                    jp = (j+LENGTH-1) % LENGTH;
                    jm = (j+1) % LENGTH;
                    kp = (k+LENGTH-1) % LENGTH;
                    km = (k+1) % LENGTH;
                    //Current spin value
                    currentSpin[0] = lattice[i][j][k][0];
                    currentSpin[1] = lattice[i][j][k][1];
                    //Propose new spin value
                    newSpin[0] = rng_theta();
                    newSpin[1] = rng_phi();
                    //Current angles
                    theta = currentSpin[0];
                    phi = currentSpin[1];

                    //MC for Theta
                    //(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
                    //deltaE = -J (Si_new – Si_old)*(S1 + S2+ S3+ S4+ S5+ S6)
                    //Test only theta, so keep phiPrime the same
                    thetaPrime = theta;
                    phiPrime = currentSpin[1];
                    //Left side of dot product
                    spinDiff[0] = sin(thetaPrime)*cos(phiPrime) - sin(theta)*cos(phi);
                    spinDiff[1] = sin(thetaPrime)*sin(phiPrime) - sin(theta)*sin(phi);
                    spinDiff[2] = cos(thetaPrime) - cos(theta);
                    //Right side dot product
                    spin1theta = lattice[im][j][k][0];
                    spin1phi   = lattice[im][j][k][1];
                    spin2theta = lattice[ip][j][k][0];
                    spin2phi   = lattice[ip][j][k][1];
                    spin3theta = lattice[i][jm][k][0];
                    spin3phi   = lattice[i][jm][k][1];
                    spin4theta = lattice[i][jp][k][0];
                    spin4phi   = lattice[i][jp][k][1];
                    spin5theta = lattice[i][j][km][0];
                    spin5phi   = lattice[i][j][km][1];
                    spin6theta = lattice[i][j][kp][0];
                    spin6phi   = lattice[i][j][kp][1];
                    //Putting the right side of the dot product together....
                    //                     Spin 1                              Spin 2                            Spin 3                             Spin 4                            Spin 5                            Spin 6
                    spinSum[0] = sin(spin1theta)*cos(spin1phi) +    sin(spin2theta)*cos(spin2phi) +     sin(spin3theta)*cos(spin3phi) +   sin(spin4theta)*cos(spin4phi) +   sin(spin5theta)*cos(spin5phi) +   sin(spin6theta)*cos(spin6phi);
                    spinSum[1] = sin(spin1theta)*sin(spin1phi) +    sin(spin2theta)*sin(spin2phi) +     sin(spin3theta)*sin(spin3phi) +   sin(spin4theta)*sin(spin4phi) +   sin(spin5theta)*sin(spin5phi) +   sin(spin6theta)*sin(spin6phi);
                    spinSum[2] = cos(spin1theta) +                  cos(spin2theta) +                   cos(spin3theta) +                 cos(spin4theta) +                 cos(spin5theta) +                 cos(spin6theta);
                    //Compute dot product
                    dotProduct = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
                    //Change in energy
                    deltaE = -J*dotProduct;

                    //Acceptance check
                    test++;
                    if (deltaE <= 0)
                    {
                        lattice[i][j][k][0] = thetaPrime;
                        lattice[i][j][k][1] = phiPrime;
                        accept++;
                    }

                    else //deltaE > 0
                    {
                        prob = (sin(thetaPrime)/sin(theta))*exp(-beta*deltaE);
                        randNum = rng();
                        if (randNum <= prob)
                        {
                            lattice[i][j][k][0] = thetaPrime;
                            lattice[i][j][k][1] = phiPrime;
                            accept++;
                        }
                        //else //randNum > probability
                        //    lattice[i][j][k] = oldSpin;
                    }

                    //MC for Phi
                    //Test only phi
                    thetaPrime = theta; //MC for phi only; theta doesn't change
                    phiPrime = newSpin[1];
                    //Left side of dot product
                    spinDiff[0] = sin(thetaPrime)*cos(phiPrime) - sin(theta)*cos(phi);
                    spinDiff[1] = sin(thetaPrime)*sin(phiPrime) - sin(theta)*sin(phi);
                    spinDiff[2] = cos(thetaPrime) - cos(theta);
                    //Right side dot product
                    spin1theta = lattice[im][j][k][0];
                    spin1phi   = lattice[im][j][k][1];
                    spin2theta = lattice[ip][j][k][0];
                    spin2phi   = lattice[ip][j][k][1];
                    spin3theta = lattice[i][jm][k][0];
                    spin3phi   = lattice[i][jm][k][1];
                    spin4theta = lattice[i][jp][k][0];
                    spin4phi   = lattice[i][jp][k][1];
                    spin5theta = lattice[i][j][km][0];
                    spin5phi   = lattice[i][j][km][1];
                    spin6theta = lattice[i][j][kp][0];
                    spin6phi   = lattice[i][j][kp][1];
                    //Putting the right side of the dot product together....
                    //                     Spin 1                              Spin 2                            Spin 3                             Spin 4                            Spin 5                            Spin 6
                    spinSum[0] = sin(spin1theta)*cos(spin1phi) +    sin(spin2theta)*cos(spin2phi) +     sin(spin3theta)*cos(spin4phi) +   sin(spin4theta)*cos(spin4phi) +   sin(spin5theta)*cos(spin5phi) +   sin(spin6theta)*cos(spin6phi);
                    spinSum[1] = sin(spin1theta)*sin(spin1phi) +    sin(spin2theta)*sin(spin2phi) +     sin(spin3theta)*sin(spin4phi) +   sin(spin4theta)*sin(spin4phi) +   sin(spin5theta)*sin(spin5phi) +   sin(spin6theta)*sin(spin6phi);
                    spinSum[2] = cos(spin1theta) +                  cos(spin2theta) +                   cos(spin3theta) +                 cos(spin4theta) +                 cos(spin5theta) +                 cos(spin6theta);
                    //Compute dot product
                    dotProduct = spinDiff[0]*spinSum[0] + spinDiff[1]*spinSum[1] + spinDiff[2]*spinSum[2];
                    //Change in energy
                    deltaE = -J*dotProduct;
                    //deltaE = -J*(cos(newSpin)-cos(curretSpin));
                    //Acceptance check
                    test++;
                    if (deltaE <= 0)
                    {
                        lattice[i][j][k][0] = thetaPrime; //Theta coordinate doesn't change
                        lattice[i][j][k][1] = phiPrime;
                        accept++;
                    }

                    else //deltaE > 0
                    {
                        prob = (sin(thetaPrime)/sin(theta))*exp(-beta*deltaE);
                        randNum = rng();
                        if (randNum <= prob)
                        {
                            lattice[i][j][k][0] = thetaPrime; //Theta coordinate doesn't change
                            lattice[i][j][k][1] = phiPrime;
                            accept++;
                        }
                        //else //randNum > probability; no change
                        //    lattice[i][j][k][0] = currentSpin[0];
                        //    lattice[i][j][k][1] = currentSpin[2];
                    }
                } //sw3
                } //sw2
                } //sw1
            } //sw
        } //End skips
        //ENERGY
        //MAGNETIZATION AND SUSCEPTIBILITY
            //Energy sweep - calculate energy of lattice in single sweep
            //Reset e_storage magnetization calc variables
            e_storage = 0;
            Mxsum = 0.0;
            Mysum = 0.0;
            Mzsum = 0.0;
            //Start calculations loops
            for (sw1 = 0; sw1 < LENGTH; sw1++)
            {
            for (sw2 = 0; sw2 < LENGTH; sw2++)
            {
            for (sw3 = 0; sw3 < LENGTH; sw3++)
            {
                //Select random spin's location
                i = sw1;
                j = sw2;
                k = sw3;
                //Six nearest neighbors
                ip = (i+LENGTH-1) % LENGTH;
                im = (i+1) % LENGTH;
                jp = (j+LENGTH-1) % LENGTH;
                jm = (j+1) % LENGTH;
                kp = (k+LENGTH-1) % LENGTH;
                km = (k+1) % LENGTH;
                //Current angles
                theta = lattice[i][j][k][0];
                phi = lattice[i][j][k][1];
                //printf("Theta, phi: %f,%f\n",theta,phi);
                //ENERGY
                //Dot product calculation
                //E = -J (Si_current)*(S1 + S2+ S3+ S4+ S5+ S6)
                //NOTE: (Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
                //Left side of dot product
                spinX = sin(theta)*cos(phi);
                spinY = sin(theta)*sin(phi);
                spinZ = cos(theta);
                //Right side dot product
                spin1theta = lattice[im][j][k][0];
                spin1phi   = lattice[im][j][k][1];
                spin2theta = lattice[ip][j][k][0];
                spin2phi   = lattice[ip][j][k][1];
                spin3theta = lattice[i][jm][k][0];
                spin3phi   = lattice[i][jm][k][1];
                spin4theta = lattice[i][jp][k][0];
                spin4phi   = lattice[i][jp][k][1];
                spin5theta = lattice[i][j][km][0];
                spin5phi   = lattice[i][j][km][1];
                spin6theta = lattice[i][j][kp][0];
                spin6phi   = lattice[i][j][kp][1];
                //Putting the right side of the dot product together....
                //                     Spin 1                              Spin 2                            Spin 3                             Spin 4                            Spin 5                            Spin 6
                spinSum[0] = sin(spin1theta)*cos(spin1phi) +    sin(spin2theta)*cos(spin2phi) +     sin(spin3theta)*cos(spin3phi) +   sin(spin4theta)*cos(spin4phi) +   sin(spin5theta)*cos(spin5phi) +   sin(spin6theta)*cos(spin6phi);
                spinSum[1] = sin(spin1theta)*sin(spin1phi) +    sin(spin2theta)*sin(spin2phi) +     sin(spin3theta)*sin(spin3phi) +   sin(spin4theta)*sin(spin4phi) +   sin(spin5theta)*sin(spin5phi) +   sin(spin6theta)*sin(spin6phi);
                spinSum[2] = cos(spin1theta) +                  cos(spin2theta) +                   cos(spin3theta) +                 cos(spin4theta) +                 cos(spin5theta) +                 cos(spin6theta);
                //Compute dot product
                dotProduct = spinX*spinSum[0] + spinY*spinSum[1] + spinZ*spinSum[2];
                //Energy summation calculation
                e_storage = e_storage + -J*dotProduct;

                //Correct for double counting
                //MAGNETIZATION
                //m_rms = M = (1/N)sqrt(sum(S_i*S_i))
                //NOTE: (Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)cos(phi), cos(theta))
                //Taken from Landau and Binder p.153, (5.16)
                X = spinX;
                Y = spinY;
                Z = spinZ;
                Mxsum += X;
                Mysum += Y;
                Mzsum += Z;
            } //sw1
            } //sw2
            } //sw3
            //Calculate energy
            e[eq] = e_storage/2.0; //Corrected for double counting
            e2[eq] = e[eq] * e[eq];//E^2
            //printf("e[eq]: %f\n",e[eq]);
            //Calculate magnetization
            Mx[eq] = Mxsum/N_;
            My[eq] = Mysum/N_;
            Mz[eq] = Mzsum/N_;
            m[eq] = sqrt(Mx[eq]*Mx[eq] + My[eq]*My[eq] + Mz[eq]*Mz[eq]);
            m2[eq] = m[eq] * m[eq];
        //Data debug
        if (DATA_DEBUG){printf("\t\tMeasurement runs (0-%d): %d\n",RUNS,eq);}
    } //End measurement phase

    //Array summations
    eSum = 0;
    e2Sum = 0;
    mSum = 0;
    m2Sum = 0;
    SxSum = 0;
    SySum = 0;
    SzSum = 0;
    for (ii = 0; ii < RUNS_; ii++)
    {
        printf("eSum 1: %f\n",eSum);
        eSum += e[ii];
        printf("eSum 2: %f\n",eSum);
        e2Sum += e2[ii];
        mSum += m[ii];
        m2Sum += m2[ii];
        SxSum += Mx[ii];
        SySum += My[ii];
        SzSum += Mz[ii];
    }

    //Spin averages
    //<Si>, i = x, y, or z
    SxAvg = SxSum / RUNS_;
    SyAvg = SySum / RUNS_;
    SzAvg = SzSum / RUNS_;
    //<Si^2> or SiAvg*SiAvg, i = x, y, or z
    SxAvgSqr = (SxSum * SxSum) / RUNS_;
    SyAvgSqr = (SxSum * SxSum) / RUNS_;
    SzAvgSqr = (SxSum * SxSum) / RUNS_;

    //Calculations
    data[t][0] = temperature;
    //<E>
    data[t][1] = eSum / RUNS_;
    //<E^2>
    data[t][2] = e2Sum / RUNS_;
    //<E>^2
    data[t][3] = data[t][1]*data[t][1];
    //<M>
    data[t][4] = mSum / RUNS_;
        //<M^2>
        data[t][5] = m2Sum / RUNS_;
        //<M>^2
        data[t][6] = data[t][4]*data[t][4];
    //<E>/N
    data[t][7] = data[t][1] / N_;
    //<M>/N
    data[t][8] = data[t][4] / N_;
    //C
    data[t][9] = (data[t][2] - data[t][3]) / (temperature*temperature*N_);
    //C/N
    data[t][10] = data[t][9] / N_;
    //X
    //Chi = beta*[<sx^2> + <sy^2> + <sz^2>  -  <sx> - <sy> - <sz> ]
    //where <sx> = 1/N sum_i Sx  and <sx^2> = (1/N sum_i Sx )^2
    //Chi = (1.0/temperature)*(SxAvg*SxAvg + SyAvg*SyAvg + SzAvg*SzAvg - SxAvg - SyAvg - SzAvg);
    data[t][11] = (1.0/temperature)*(SxAvgSqr + SyAvgSqr + SzAvgSqr - SxAvg - SyAvg -SzAvg);
    printf("Chi: %f\n", data[t][11]);
    //X/N
    data[t][12] = data[t][11] / N_;
    //Acceptance Ratio
    printf("\t\tAcceptance Ratio: %f\n",accept/test);
}

//Print data to a file
void print_to_file(double data[][13])
{
     FILE *file; //For writing files
    int i; //Counter

    //Data debug
    if (DATA_DEBUG)
    printf("\t\tTemp\t\t<E/N>\t\tC/N\n");
    {
        for(i = 0; i < MAX_TEMP_STEPS; i++)
            printf("\t\t%f\t\t%f\t\t%f\n",data[i][0],data[i][7],data[i][10]);
    }

    //Average energy
    file = fopen("e_avg.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][1]);
    }
    fclose(file);

    //Average energy per spin
    file = fopen("e_avg_persite.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][7]);
    }
    fclose(file);

    //Average magnetization
    file = fopen("m_avg.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][4]);
    }
    fclose(file);

    //Average magnetization per site
    file = fopen("m_avg_persite.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][8]);
    }
    fclose(file);

    //Heat capacity
    file = fopen("c.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][9]);
    }
    fclose(file);

    //Heat capacity per site
    file = fopen("c_persite.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][10]);
    }
    fclose(file);

    //Magnetic susceptibility
    file = fopen("x.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][11]);
    }
    fclose(file);

    //Magnetic susceptibility per site
    file = fopen("x_persite.dat","w+");
    for(i = 0; i < MAX_TEMP_STEPS; i++)
    {
        fprintf(file, "%f\t%f\n",data[i][0],data[i][12]);
    }
    fclose(file);
}
