/***************************************************************************/
/*                                                                         */
/* The sor06.C main program                                                */
/*  Version: 1.0                                                           */
/*  Date: 13 Sept. 2019                                                    */
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
    int HB, total_vessels, total_terminal, total_conn,number_of_points;
    
    //==========adjustment in resistances (WK parameters) for the control========
    
    
    if (argc != 7) //argv[0] is the name of the program, here sor06
    {
        printf("Not enough input arguments, noargc %d and they are %s\n", argc, argv[0]);
        return 1;
    }
    
    f1 = atof(argv[1]);
    f2 = atof(argv[2]);
    f3 = atof(argv[3]);
    HB = atoi(argv[4]);
    total_vessels  = atoi(argv[5]);
    total_terminal = atoi(argv[6]);
    total_conn     = total_vessels-total_terminal;
    nbrves             = total_vessels;
    number_of_points  = 8;
    
    
    
    /* Declare string vectors and files to hold 
     * the output from the model               */
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
    finaltime = HB*Period;      // Final end-time during a simulation.
    tend      = (HB-1)*Period;  // Timestep before the first plot-point
                                // is reached.
    
    // The number of vessels in the network is given when the governing array of
    // vessels is declared.
    
// =========================NETWORK =================================================================================
    
    Tube   *Arteries[nbrves];                    // Array of blood vessels.
    int conn_rows = (total_vessels-1)/2;
    int conn_cols = 3;
    int connectivity_matrix[conn_rows][conn_cols];
    double bc_matrix[total_terminal][3];
    int terminal_vessels[total_terminal];
    double dimensions_matrix[total_vessels][3];
    FILE *conn;
    conn = fopen("connectivity.txt","rt");
    int parent, daughter1, daughter2, r_in;
    int conn_id = 0;
    
    // Check to see if we have the connecitivity file
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
        //printf("BCs:   %lf    %lf    %lf     \n",bc_matrix[i][0],bc_matrix[i][1],bc_matrix[i][2]);
    }
    
    fclose(BCs);
    

    /////////////// Load in terminal vessels
    
    FILE *terminal_file;
    terminal_file = fopen("terminal_vessels.txt","rt");
    for (int i=0; i<total_terminal; i++){
        fscanf(terminal_file, "%d", &terminal_vessels[i]);
        //printf("terminal vessels: %d\n",terminal_vessels[i]);
    }
    fclose(terminal_file);
    
    
    /////////////// Load in Dimensions
    FILE *dim_file;
    dim_file = fopen("Dimensions.txt","rt");
    for (int dim_ID = 0; dim_ID < total_vessels; dim_ID++){
        fscanf(dim_file, "%5lf %5lf %5lf", &dimensions_matrix[dim_ID][0],&dimensions_matrix[dim_ID][1],&dimensions_matrix[dim_ID][2]);
        printf("Length: %lf    Rin: %lf   Rout: %lf\n", dimensions_matrix[dim_ID][0],dimensions_matrix[dim_ID][1],dimensions_matrix[dim_ID][2]);
    }
    fclose(dim_file);
    
    
    // Initialization of the Arteries.
    // Definition of Class Tube: (Length, topradius, botradius, *LeftDaughter, *RightDaughter, rmin, points, init, K, f1, f2, f3, R1, R2,  CT, LT);
    /////////////// Create a dynamic network (Arteries[1] = new Tube( L1, R1, R1, Arteries[ 2], Arteries[ 3], rm, 40, 0, 0,f1,f2,f3, 0, 0, 0, 0);)
    /*                 FORWARD VERSION             */
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

    
    // Solves the equations until time equals tend.
    solver (Arteries, tstart, tend, k);
    tstart = tend;
    tend = tend + Deltat;
    
    // The loop is continued until the final time
    // is reached. If one wants to make a plot of
    // the solution versus x, tend is set to final-
    // time in the above declaration.
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
        solver (Arteries, tstart, tend, k);
        
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
