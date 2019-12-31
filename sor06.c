/***************************************************************************/
/*                                                                         */
/* The sor06.C main program                                                */
/*  Version: 1.0                                                           */
/*  Date: 30 Dec. 2019                                                     */
/*                                                                         */
/*  Primary Authors: M.S. Olufsen                                          */
/*  Key Contributers: M.U. Qureshi & M.J. Colebank                         */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/***************************************************************************/

#include "sor06.h"
#include "tools.h"
#include "arteries.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determed. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{
    double tstart, tend, finaltime;
    
    double f1, f2, f3;
    double Tper, Period, k, Deltat;
    int total_vessels, total_terminal, total_conn,number_of_points;
    
    //==========adjustment in resistances (WK parameters) for the control========
    
    
    if (argc != 8) //argv[0] is the name of the program, here sor06
    {
        printf("Not enough input arguments: only noargc %d but require %d\n", argc, 8);
        return 1;
    }
    
    f1   = atof(argv[1]);
    f2   = atof(argv[2]);
    f3   = atof(argv[3]);
    Tper = atof(argv[4]);
    
    
    total_vessels    = atoi(argv[5]);
    total_terminal   = atoi(argv[6]);
    number_of_points = atoi(argv[7]);
    
    total_conn         = total_vessels-total_terminal;
    nbrves             = total_vessels;
    
    printf("f1 f2 f3 Tper: %lf %lf %lf %lf,\n vessels, terminal, numpts: %d, %d, %d\n",f1, f2, f3, Tper, total_vessels, total_terminal, number_of_points);
    
    // Originally in sor06.h
    Period = Tper*q/Lr3;            // The dimension-less period.
    k      = Period/tmstps;         // Length of a timestep.
    Deltat = Period/plts;           // Interval between each point plottet.
    
    /* Declare string vectors and files to hold 
     * the output from the model. This can 
     * be expanded to have more vessels
     */
    char nameart1 [20];
    char nameart2 [20];
    char nameart3 [20];
    
    sprintf(nameart1, "art1_1.2d");
    FILE *art1 = fopen (nameart1, "w");
    
    sprintf(nameart2, "art2_1.2d");
    FILE *art2 = fopen (nameart2, "w");
    
    sprintf(nameart3, "art3_1.2d");
    FILE *art3 = fopen (nameart3, "w");

    
    
    
    // Workspace used by bound_bif
    for(int i=0; i<18; i++) fjac[i] = new double[18];
    tstart    = 0.0;            // Starting time.
    
    // The number of vessels in the network is given when the governing array of
    // vessels is declared.
    
// =========================NETWORK =================================================================================
    
    Tube   *Arteries[nbrves];                     // Array of blood vessels.
    int conn_rows = (total_vessels-1)/2;          // Number of rows
    int conn_cols = 3;                            // Number of columns
    int connectivity_matrix[conn_rows][conn_cols];// Matrix with vessels
    double bc_matrix[total_terminal][3];          // WK bound. matrix
    int terminal_vessels[total_terminal];         // ID's for term. ves.
    double dimensions_matrix[total_vessels][3];   // Length and radius
    
    FILE *conn;
    conn = fopen("connectivity.txt","rt");
    int parent, daughter1, daughter2, r_in;
    int conn_id = 0;
    
    // Check to see if we have the connectivity file
   if (conn == NULL)
   {
       fprintf(stdout,"Error: Connectivity File Does Not Exist \n");
       return 1;
   }
    while ((r_in = fscanf(conn, "%d %d %d", &parent, &daughter1, &daughter2)) != EOF)
    {
        connectivity_matrix[conn_id][0] = parent;
        connectivity_matrix[conn_id][1] = daughter1;
        connectivity_matrix[conn_id][2] = daughter2;
        conn_id++;
    }
    fclose(conn);
    
       /////////////// Load in the 3-WK Boundary Conditions Matrix
    FILE *BCs;
    BCs = fopen("Windkessel_Parameters.txt","rt");
   if (BCs==NULL)
   {
       fprintf(stdout,"Error: Boundary Conditions File Does Not Exist");
       return 1;
   }
    for (int i=0; i<total_terminal; i++){
        fscanf(BCs, "%lf %lf %lf",&bc_matrix[i][0],&bc_matrix[i][1],&bc_matrix[i][2]);
    }
    
    fclose(BCs);
    

    /////////////// Load in terminal vessels
    
    FILE *terminal_file;
    terminal_file = fopen("terminal_vessels.txt","rt");
    for (int i=0; i<total_terminal; i++){
        fscanf(terminal_file, "%d", &terminal_vessels[i]);
    }
    fclose(terminal_file);
    
    
    /////////////// Load in Dimensions
    FILE *dim_file;
    dim_file = fopen("Dimensions.txt","rt");
    for (int dim_ID = 0; dim_ID < total_vessels; dim_ID++){
        fscanf(dim_file, "%5lf %5lf %5lf", &dimensions_matrix[dim_ID][0],&dimensions_matrix[dim_ID][1],&dimensions_matrix[dim_ID][2]);
    }
    fclose(dim_file);
    
    
    // Initialization of the Arteries.
    // Definition of Class Tube: (Length, topradius, botradius, *LeftDaughter, *RightDaughter,
    //                              points, init, f1, f2, f3, R1, R2,  CT);
    int curr_d1, curr_d2;
    double R1, R2, CT;
    int term_id = total_terminal-1;
    int bc_id = total_terminal-1;
    conn_id = total_conn-1;
    
///// A hard code for a single vessel
    if (total_vessels == 1){
        R1 = bc_matrix[0][0];
        R2 = bc_matrix[0][1];
        CT = bc_matrix[0][2];
        Arteries[0] = new Tube( dimensions_matrix[0][0], dimensions_matrix[0][1], dimensions_matrix[0][2], 0, 0, number_of_points, 1,f1,f2,f3, R1, R2, CT);
    }
    else if (total_vessels > 1){
        // NOTE: In order to set up the tube class, vessels must be initilized
        // starting from the most termin vessel (the vessel with the largest
        // ID number), and then works back to the root vessel in the network.
        for (int i=total_vessels-1; i>=0; i--) {
            //// Define the connecitivity between the vessels
            if (i == connectivity_matrix[conn_id][0])
            {
                curr_d1     = connectivity_matrix[conn_id][1];
                curr_d2     = connectivity_matrix[conn_id][2];
                conn_id--;
            }
            //first vessels
            if (i==0) {
                printf("Root Vessel: L, Rin, Rout: %lf %lf %lf\n",dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][2]);
                Arteries[i] = new Tube( dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][2], Arteries[ curr_d1], Arteries[ curr_d2], number_of_points, 1, f1,f2,f3, 0, 0, 0);
            }
            else{
                if (i == terminal_vessels[term_id] && term_id>=0){
                    R1 = bc_matrix[bc_id][0];
                    R2 = bc_matrix[bc_id][1];
                    CT = bc_matrix[bc_id][2];
                    Arteries[i] = new Tube( dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][2], 0, 0, number_of_points, 0, f1,f2,f3, R1, R2, CT);
                    term_id--;
                    bc_id--;
                }
                else
                {
                    Arteries[i] = new Tube( dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][2], Arteries[ curr_d1], Arteries[ curr_d2], number_of_points, 0, f1,f2,f3, 0, 0, 0);
                }
            }
        }
    }

    
    // Solves the equations until time equals tend.///////////////////////////
    /* ADDED BY M. J. Colebank
     * Rather than specifying the number of cycles as an input to the function,
     * we want to test to see if the solution has converged. If so, we should exit.*/
    
    int period_counter = 1; // Count the number of periods you have solved for
    double norm_sol = 1e+6;
    double sol_tol  = 1e+1;
    printf("NORM_SOL: %f\n",norm_sol);
    double sol_p1[tmstps],sol_p2[tmstps];
    tend      = Deltat;
    
    
    // SOLVE THE MODEL ONCE
    // Note: Only want to test the pressure at the inlet
    int sol_ID = 0;
    while (tend<=period_counter*Period)
    {
    solver (Arteries, tstart, tend, k, Period); 
    sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0]); // for printing
    sol_p1[sol_ID] *= rho*g*Lr/conv;
    tstart = tend;
    tend   = tend + Deltat; // The current ending time is increased by Deltat.
    sol_ID++;
    }
    
    
    // LOOP FOR CONVERGENCE
    double sse;
    while (norm_sol>=sol_tol)
    {
        sol_ID = 0;
        sse    = 0;
        period_counter++;
        if (period_counter>max_cycles)
        {
            printf("ERROR: TOO MANY CYCLES. EXITING. \n");
            return 1;
        }
        while (tend<=period_counter*Period)
        {
            solver (Arteries, tstart, tend, k, Period);
            sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0]); // for printing
            sol_p2[sol_ID] *= rho*g*Lr/conv;
            sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
            tstart = tend;
            tend   = tend + Deltat; // The current ending time is increased by Deltat.
            sol_ID++;
        }
        norm_sol = sse;
        memcpy (sol_p1, sol_p2, sizeof(sol_p2));
        printf("NORM OF SOLN DIFF:%f\n",norm_sol);
    }
    printf("final num cylces:%d\n",period_counter);
    

    //////////////////////////////////////////////
  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration.
  
  period_counter++;
  finaltime = (period_counter+(cycles-1))*Period;
    while (tend <= finaltime)
    {
        for (int j=0; j<nbrves; j++)
        {
            int ArtjN = Arteries[j]->N;
            for (int i=0; i<ArtjN; i++)
            {
                Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
                Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
            }
        }
        
        // Solves the equations until time equals tend.
        solver (Arteries, tstart, tend, k, Period);
        
        // A 2D plot of P(x_fixed,t) is made. The resulting 2D graph is due to
        // the while loop, since for each call of printPt only one point is set.
            Arteries[ 0] -> printPxt (art1, tend, 0);
        if (total_vessels>1) {
            Arteries[ connectivity_matrix[0][1]] -> printPxt (art2, tend, 0);
            Arteries[ connectivity_matrix[0][2]] -> printPxt (art3, tend, 0);
        }
       

        // The time within each print is increased.
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
    }
    
    // In order to termate the program correctly the vessel network and hence
    // all the vessels and their workspace are deleted.
    for (int i=0; i<nbrves; i++) delete Arteries[i];
    
    // Matrices and arrays are deleted
    for (int i=0; i<18; i++) delete[] fjac[i];
    
    fclose (art1);
    if (total_vessels>1) {
    fclose (art2);
    fclose (art3);
    }
    return 0;
}
