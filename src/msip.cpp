#include <cstdio>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <cmath>

using namespace std;

enum MSIP_N
{
    MSIP5,
    MSIP9
};

struct Params
{
    int imax;
    int jmax;
    double tol;
    double psi;
    double omega;
    int maxiter;
};


void msip5(double** A, double psi, int N, int jmax)
{
    int k=0;
    A[k][3] = A[k][3]/A[k][2];
    A[k][4] = A[k][4]/A[k][2];
    for (k=1;k<jmax;k++)
    {
        A[k][1] = A[k][1]/(1+psi*A[k-1][4]);
        A[k][2] = A[k][2] - A[k][1]*A[k-1][3] + psi*A[k][1]*A[k-1][4];
        A[k][3] = A[k][3]/A[k][2];
        A[k][4] = (A[k][4] - psi*A[k][1]*A[k-1][4])/A[k][2];           
    }
    for (k=jmax;k<N-jmax;k++)
    {
        A[k][0] = A[k][0]/(1+psi*A[k-jmax][3]);
        A[k][1] = A[k][1]/(1+psi*A[k-1][4]);
        A[k][2] = A[k][2] - A[k][1]*A[k-1][3] - A[k][0]*A[k-jmax][4] + psi*(A[k][1]*A[k-1][4] + A[k][0]*A[k-jmax][3]);
        A[k][3] = (A[k][3] - psi*A[k][0]*A[k-jmax][3])/A[k][2];
        A[k][4] = (A[k][4] - psi*A[k][1]*A[k-1][4])/A[k][2];
    }    
    for (k=N-jmax;k<N-1;k++)
    {
        A[k][0] = A[k][0]/(1+psi*A[k-jmax][3]);
        A[k][1] = A[k][1]/(1+psi*A[k-1][4]);
        A[k][2] = A[k][2] - A[k][1]*A[k-1][3] - A[k][0]*A[k-jmax][4] + psi*(A[k][1]*A[k-1][4] + A[k][0]*A[k-jmax][3]);
        A[k][3] = (A[k][3] - psi*A[k][0]*A[k-jmax][3])/A[k][2];    
    }
    A[k][0] = A[k][0]/(1+psi*A[k-jmax][3]);
    A[k][1] = A[k][1]/(1+psi*A[k-1][4]);
    A[k][2] = A[k][2] - A[k][1]*A[k-1][3] - A[k][0]*A[k-jmax][4] + psi*(A[k][1]*A[k-1][4] + A[k][0]*A[k-jmax][3]);
}


void msip9(double** LU, double** A,double psi, int N, int jmax) {
    int k=0;
    LU[k][4] = A[k][4];
    LU[k][5] = A[k][5]/LU[k][4];
    LU[k][6] = A[k][6]/LU[k][4];
    LU[k][7] = A[k][7]/LU[k][4];
    LU[k][8] = A[k][8]/LU[k][4];
    for (k=1;k<jmax-1;k++) {
        LU[k][3] = A[k][3]/(1+2*psi*LU[k-1][6]);
        LU[k][4] = A[k][4] - LU[k][3]*LU[k-1][5] + 2*psi*LU[k][3]*LU[k-1][6];
        LU[k][5] = A[k][5]/LU[k][4];
        LU[k][6] = (A[k][6] - LU[k][3]*LU[k-1][7])/LU[k][4];
        LU[k][7] = (A[k][7] - psi*LU[k][3]*LU[k-1][6] - LU[k][3]*LU[k-1][8])/LU[k][4];
        LU[k][8] = A[k][8]/LU[k][4];
    }
    LU[k][2] = A[k][2];
    LU[k][3] = A[k][3]/(1+2*psi*LU[k-1][6]);
    LU[k][4] = A[k][4] - LU[k][3]*LU[k-1][5] - LU[k][2]*LU[k-jmax+1][6] + 2*psi*(LU[k][2]*LU[k-jmax+1][5] + LU[k][3]*LU[k-1][6]) + psi*LU[k][2]*LU[k-jmax+1][8];
    LU[k][5] = (A[k][5] - 2*psi*LU[k][2]*(LU[k-jmax+1][5] + LU[k-jmax+1][8]) - LU[k][2]*LU[k-jmax+1][7])/LU[k][4];
    LU[k][6] = (A[k][6] - LU[k][3]*LU[k-1][7])/LU[k][4];
    LU[k][7] = (A[k][7] - psi*LU[k][3]*LU[k-1][6] - LU[k][3]*LU[k-1][8])/LU[k][4];
    LU[k][8] = A[k][8]/LU[k][4];
    k = jmax;
    LU[k][1] = (A[k][1] - psi*LU[k-jmax+1][5]*A[k][2])/(1-psi*LU[k-jmax][5]*LU[k-jmax+1][5]);
    LU[k][2] = A[k][2] - LU[k][1]*LU[k-jmax][5]; 
    LU[k][3] = (A[k][3] - LU[k][1]*LU[k-jmax][6])/(1+2*psi*LU[k-1][6]);
    LU[k][4] = A[k][4] - LU[k][2]*LU[k-jmax+1][6] - LU[k][3]*LU[k-1][5] - LU[k][1]*LU[k-jmax][7] + 2*psi*(LU[k][2]*LU[k-jmax+1][5] + LU[k][3]*LU[k-1][6]) + psi*LU[k][2]*LU[k-jmax+1][8];
    LU[k][5] = (A[k][5] - 2*psi*LU[k][2]*(LU[k-jmax+1][5] + LU[k-jmax+1][8]) - LU[k][2]*LU[k-jmax+1][7] - LU[k][1]*LU[k-jmax][8])/LU[k][4];
    LU[k][6] = (A[k][6] - LU[k][3]*LU[k-1][7])/LU[k][4]; 
    LU[k][7] = (A[k][7] - psi*LU[k][3]*LU[k-1][6] - LU[k][3]*LU[k-1][8])/LU[k][4];
    LU[k][8] = A[k][8]/LU[k][4];
    for (k=jmax+1;k<N-jmax-1;k++) {
        LU[k][0] = A[k][0];
        LU[k][1] = (A[k][1] - psi*LU[k-jmax+1][5]*A[k][2] - LU[k][0]*LU[k-jmax-1][5])/(1-psi*LU[k-jmax][5]*LU[k-jmax+1][5]);
        LU[k][2] = A[k][2] - LU[k][1]*LU[k-jmax][5];
        LU[k][3] = (A[k][3] - 2*psi*LU[k][0]*LU[k-jmax-1][6] - LU[k][0]*LU[k-jmax-1][7] - LU[k][1]*LU[k-jmax][6])/(1+2*psi*LU[k-1][6]);
        LU[k][4] = A[k][4] - LU[k][0]*LU[k-jmax-1][8] - LU[k][1]*LU[k-jmax][7] - LU[k][2]*LU[k-jmax+1][6] - LU[k][3]*LU[k-1][5] + 2*psi*(LU[k][2]*LU[k-jmax+1][5] + LU[k][3]*LU[k-1][6]) + psi*(LU[k][0]*LU[k-jmax-1][6] + LU[k][2]*LU[k-jmax+1][8]);
        LU[k][5] = (A[k][5] - 2*psi*LU[k][2]*(LU[k-jmax+1][5] + LU[k-jmax+1][8]) - LU[k][2]*LU[k-jmax+1][7] - LU[k][1]*LU[k-jmax][8])/LU[k][4];
        LU[k][6] = (A[k][6] - LU[k][3]*LU[k-1][7])/LU[k][4]; 
        LU[k][7] = (A[k][7] - psi*LU[k][3]*LU[k-1][6] - LU[k][3]*LU[k-1][8])/LU[k][4];
        LU[k][8] = A[k][8]/LU[k][4]; 
    }
    LU[k][0] = A[k][0];
    LU[k][1] = (A[k][1] - psi*LU[k-jmax+1][5]*A[k][2] - LU[k][0]*LU[k-jmax-1][5])/(1-psi*LU[k-jmax][5]*LU[k-jmax+1][5]);
    LU[k][2] = A[k][2] - LU[k][1]*LU[k-jmax][5];
    LU[k][3] = (A[k][3] - 2*psi*LU[k][0]*LU[k-jmax-1][6] - LU[k][0]*LU[k-jmax-1][7] - LU[k][1]*LU[k-jmax][6])/(1+2*psi*LU[k-1][6]);
    LU[k][4] = A[k][4] - LU[k][0]*LU[k-jmax-1][8] - LU[k][1]*LU[k-jmax][7] - LU[k][2]*LU[k-jmax+1][6] - LU[k][3]*LU[k-1][5] + 2*psi*(LU[k][2]*LU[k-jmax+1][5] + LU[k][3]*LU[k-1][6]) + psi*(LU[k][0]*LU[k-jmax-1][6] + LU[k][2]*LU[k-jmax+1][8]);
    LU[k][5] = (A[k][5] - 2*psi*LU[k][2]*(LU[k-jmax+1][5] + LU[k-jmax+1][8]) - LU[k][2]*LU[k-jmax+1][7] - LU[k][1]*LU[k-jmax][8])/LU[k][4];
    LU[k][6] = (A[k][6] - LU[k][3]*LU[k-1][7])/LU[k][4]; 
    LU[k][7] = (A[k][7] - psi*LU[k][3]*LU[k-1][6] - LU[k][3]*LU[k-1][8])/LU[k][4];
    k = N-jmax;
    LU[k][0] = A[k][0];
    LU[k][1] = (A[k][1] - psi*LU[k-jmax+1][5]*A[k][2] - LU[k][0]*LU[k-jmax-1][5])/(1-psi*LU[k-jmax][5]*LU[k-jmax+1][5]);
    LU[k][2] = A[k][2] - LU[k][1]*LU[k-jmax][5];
    LU[k][3] = (A[k][3] - 2*psi*LU[k][0]*LU[k-jmax-1][6] - LU[k][0]*LU[k-jmax-1][7] - LU[k][1]*LU[k-jmax][6])/(1+2*psi*LU[k-1][6]);
    LU[k][4] = A[k][4] - LU[k][0]*LU[k-jmax-1][8] - LU[k][1]*LU[k-jmax][7] - LU[k][2]*LU[k-jmax+1][6] - LU[k][3]*LU[k-1][5] + 2*psi*(LU[k][2]*LU[k-jmax+1][5] + LU[k][3]*LU[k-1][6]) + psi*(LU[k][0]*LU[k-jmax-1][6] + LU[k][2]*LU[k-jmax+1][8]);
    LU[k][5] = (A[k][5] - 2*psi*LU[k][2]*(LU[k-jmax+1][5] + LU[k-jmax+1][8]) - LU[k][2]*LU[k-jmax+1][7] - LU[k][1]*LU[k-jmax][8])/LU[k][4];
    LU[k][6] = (A[k][6] - LU[k][3]*LU[k-1][7])/LU[k][4]; 
    for (k=N-jmax+1;k<N-1;k++) {
        LU[k][0] = A[k][0];
        LU[k][1] = (A[k][1] - psi*LU[k-jmax+1][5]*A[k][2] - LU[k][0]*LU[k-jmax-1][5])/(1-psi*LU[k-jmax][5]*LU[k-jmax+1][5]);
        LU[k][2] = A[k][2] - LU[k][1]*LU[k-jmax][5];
        LU[k][3] = (A[k][3] - 2*psi*LU[k][0]*LU[k-jmax-1][6] - LU[k][0]*LU[k-jmax-1][7] - LU[k][1]*LU[k-jmax][6])/(1+2*psi*LU[k-1][6]);
        LU[k][4] = A[k][4] - LU[k][0]*LU[k-jmax-1][8] - LU[k][1]*LU[k-jmax][7] - LU[k][2]*LU[k-jmax+1][6] - LU[k][3]*LU[k-1][5] + 2*psi*(LU[k][2]*LU[k-jmax+1][5] + LU[k][3]*LU[k-1][6]) + psi*(LU[k][0]*LU[k-jmax-1][6] + LU[k][2]*LU[k-jmax+1][8]);
        LU[k][5] = (A[k][5] - 2*psi*LU[k][2]*(LU[k-jmax+1][5] + LU[k-jmax+1][8]) - LU[k][2]*LU[k-jmax+1][7] - LU[k][1]*LU[k-jmax][8])/LU[k][4];
    }
    LU[k][0] = A[k][0];
    LU[k][1] = (A[k][1] - psi*LU[k-jmax+1][5]*A[k][2] - LU[k][0]*LU[k-jmax-1][5])/(1-psi*LU[k-jmax][5]*LU[k-jmax+1][5]);
    LU[k][2] = A[k][2] - LU[k][1]*LU[k-jmax][5];
    LU[k][3] = (A[k][3] - 2*psi*LU[k][0]*LU[k-jmax-1][6] - LU[k][0]*LU[k-jmax-1][7] - LU[k][1]*LU[k-jmax][6])/(1+2*psi*LU[k-1][6]);
    LU[k][4] = A[k][4] - LU[k][0]*LU[k-jmax-1][8] - LU[k][1]*LU[k-jmax][7] - LU[k][2]*LU[k-jmax+1][6] - LU[k][3]*LU[k-1][5] + 2*psi*(LU[k][2]*LU[k-jmax+1][5] + LU[k][3]*LU[k-1][6]) + psi*(LU[k][0]*LU[k-jmax-1][6] + LU[k][2]*LU[k-jmax+1][8]);
}

void rhs_constructor (double* rhs, double** A, double* b, double* u, int N, int jmax, int stencil) {
    int k=0;
    switch (stencil) {
        case 5:
            rhs[k] = b[k] - A[k][2]*u[k] - A[k][3]*u[k+1] - A[k][4]*u[k+jmax];
            for (k=1;k<jmax;k++)
                rhs[k] = b[k] - A[k][1]*u[k-1] - A[k][2]*u[k] - A[k][3]*u[k+1] - A[k][4]*u[k+jmax];
            for (k=jmax;k<N-jmax;k++)
                rhs[k] = b[k] - A[k][0]*u[k-jmax] - A[k][1]*u[k-1] - A[k][2]*u[k] - A[k][3]*u[k+1] - A[k][4]*u[k+jmax];        
            for (k=N-jmax;k<N-1;k++)
                rhs[k] = b[k] - A[k][0]*u[k-jmax] - A[k][1]*u[k-1] - A[k][2]*u[k] - A[k][3]*u[k+1];
            rhs[k] = b[k] - A[k][0]*u[k-jmax] - A[k][1]*u[k-1] - A[k][2]*u[k];
            break;

        default:
            rhs[k] = b[k] - A[k][4]*u[k] - A[k][5]*u[k+1] - A[k][6]*u[k+jmax-1] - A[k][7]*u[k+jmax] - A[k][8]*u[k+jmax+1];
            for (k=1;k<jmax-1;k++)
                rhs[k] = b[k] - A[k][3]*u[k-1] - A[k][4]*u[k] - A[k][5]*u[k+1] - A[k][6]*u[k+jmax-1] - A[k][7]*u[k+jmax] - A[k][8]*u[k+jmax+1]; 
            rhs[k] = b[k] - A[k][2]*u[k-jmax+1] - A[k][3]*u[k-1] - A[k][4]*u[k] - A[k][5]*u[k+1] - A[k][6]*u[k+jmax-1] - A[k][7]*u[k+jmax] - A[k][8]*u[k+jmax+1];    
            k=jmax;
            rhs[k] = b[k] - A[k][1]*u[k-jmax] - A[k][2]*u[k-jmax+1] - A[k][3]*u[k-1] - A[k][4]*u[k] - A[k][5]*u[k+1] - A[k][6]*u[k+jmax-1] - A[k][7]*u[k+jmax] - A[k][8]*u[k+jmax+1];
            for (k=jmax+1;k<N-jmax-1;k++)
                rhs[k] = b[k] - A[k][0]*u[k-jmax-1] - A[k][1]*u[k-jmax] - A[k][2]*u[k-jmax+1] - A[k][3]*u[k-1] - A[k][4]*u[k] - A[k][5]*u[k+1] - A[k][6]*u[k+jmax-1] - A[k][7]*u[k+jmax] - A[k][8]*u[k+jmax+1];
            rhs[k] = b[k] - A[k][0]*u[k-jmax-1] - A[k][1]*u[k-jmax] - A[k][2]*u[k-jmax+1] - A[k][3]*u[k-1] - A[k][4]*u[k] - A[k][5]*u[k+1] - A[k][6]*u[k+jmax-1] - A[k][7]*u[k+jmax];
            k=N-jmax;
            rhs[k] = b[k] - A[k][0]*u[k-jmax-1] - A[k][1]*u[k-jmax] - A[k][2]*u[k-jmax+1] - A[k][3]*u[k-1] - A[k][4]*u[k] - A[k][5]*u[k+1] - A[k][6]*u[k+jmax-1];
            for (k=N-jmax+1;k<N-1;k++)
                rhs[k] = b[k] - A[k][0]*u[k-jmax-1] - A[k][1]*u[k-jmax] - A[k][2]*u[k-jmax+1] - A[k][3]*u[k-1] - A[k][4]*u[k] - A[k][5]*u[k+1];
            rhs[k] = b[k] - A[k][0]*u[k-jmax-1] - A[k][1]*u[k-jmax] - A[k][2]*u[k-jmax+1] - A[k][3]*u[k-1] - A[k][4]*u[k];    
            break;
	}
}



class Solution
{
    public:
        int N;
        double** A;
        double** LU;
        double* b;
        double* du;
        double* u;
        double* u_exact;
        double* w;
        double* rhs;

        vector<double> err_values;
        vector<double> res_values;


    Solution(int N, MSIP_N msip) {
        // var stencil = msip switch {MSIP_N.MSIP5 => 5, MSIP_N.MSIP9 => 9};
        int stencil = 0x00;
        switch (msip) {
            case MSIP5: stencil = 5; break;
            case MSIP9: stencil = 9; break; 
        }
        this->N = N;
        b = new double[N];
        du = new double[N];
        u_exact = new double[N];
        u = new double[N];
        w = new double[N];
        rhs = new double[N];
        
        A = new double*[N];
        LU = new double*[N];

        for (int i=0;i<N;i++)
            LU[i] = new double[stencil];   

        for (int i=0;i<N;i++)
            A[i] = new double[stencil];        
            
        //----------Construction of system matrices and vectors----------         
        //Initiation of matrices 
        for (int i=0;i<N;i++) {
            for (int j=0;j<stencil;j++) {
                A[i][j] = 0;
                LU[i][j] = 0;
            }
            b[i] = 0;
            u_exact[i] = 0;    
            u[i] = 0;
            du[i] = 0;
            w[i] = 0;
            rhs[i] = 0;
        }
    }

    ~Solution() {
        delete b;
        delete du;
        delete u_exact;
        delete u;
        delete w;
        delete rhs;
        
        for (int i = 0; i < N; i++) 
            delete A[i];

        for (int i = 0; i < N; i++) 
            delete LU[i];
        
        delete A;
        delete LU;

        err_values.clear();
        res_values.clear();
    }
};


class Solver {
    int solver,imax,jmax,N,maxiter,stencil,i,j,k,write,iter,selection,extra_iter;
    double psi,tol,omega,dx,dy,res_init,res,err,err_init;

    // double[][] A, LU;
    // double[] b, du, u_exact, u, w, rhs;
    Solution* solution;

    bool is_constructed, is_initialized, is_decomposed, is_solved;
    
    public: Solver(Params pars, MSIP_N msip) {
        i = 0;
        j = 0;
        k = 0;
        tol = pars.tol;
        psi = pars.psi;
        omega = pars.omega;
        maxiter = pars.maxiter;

        imax = pars.imax;
        jmax = pars.jmax;
        
        if (psi<0)
            psi=0;
        if (psi>=1)
            psi = 0.9; 
        if (omega<=0)
            omega = 1;
        if (omega>2)
            omega = 1;

        imax++;             //number of nodes in the x-direction = number of intervals in the x-direction +1
        jmax++;             //number of nodes in the y-direction = number of intervals in the y-direction +1
        N = imax*jmax;      //total number of nodes 
        dx = (double)1/(imax-1);    //length of the interval in the x-direction
        dy = (double)1/(jmax-1);    //length of the interval in the y-direction
        solution = new Solution(N, msip);
        // solver = 1;
        // if (solver==1)
        // 	stencil = 5;
        // else
        // 	stencil = 9;
        int stencil = 0x00;
        switch (msip) {
            case MSIP5: stencil = 5; break;
            case MSIP9: stencil = 9; break; 
        }
    }

    ~Solver() {
        delete solution;
    }

    private:
        void decomposition() {
            double** A = solution->A;
            double** LU = solution->LU;
        
            //Incomplete LU decomposition
            switch (stencil) {
            	case 5:
            	    for (i=0;i<N;i++) 
            	        for (j=0;j<stencil;j++)
                	        LU[i][j] = A[i][j];            	         
            	    
            	    msip5(LU,psi,N,jmax); 
            	    break;

            	default:
            		for (i=0;i<N;i++)
        		        for (j=0;j<stencil;j++)
            		        LU[i][j] = 0;            		         
            		    
            		msip9(LU,A,psi,N,jmax);
            	    break;
            }           
        }

        void initialSolution() {
            double* u = solution->u;
            double* du = solution->du;
            double* w = solution->w;
            double* rhs = solution->rhs;
        
            //Initial solution 
            for (i=0;i<N;i++) {
                u[i] = 0;
                du[i] = 0;
                w[i] = 0;
                rhs[i] =0;
            }        
        }
        
        void solve() {
            double** A = solution->A;
            double** LU = solution->LU;
            double* b = solution->b;
            double* u = solution->u;
            double* w = solution->w;
            double* du = solution->du;
            double* rhs = solution->rhs;
            double* u_exact = solution->u_exact;
        
            //----------Solution of the system using the Delta Formulation----------
            rhs_constructor(rhs,A,b,u,N,jmax,stencil);
    
            res_init = 0;
            err_init = 0;
            for (i=0;i<N;i++)
            {
                res_init += pow(rhs[i],2);
                err_init += pow(u_exact[i]-u[i],2);
            }
            res_init = sqrt(res_init);      //initial (absolute) residual norm
            err_init = sqrt(err_init);      //initial (absolute) error norm
            iter = 0;                       //iteration counter
    	    switch (stencil)
            {
            case 5:
                while(true)
                {   
                    //calculation of the new rhs and the relative residual
                    rhs_constructor(rhs,A,b,u,N,jmax,stencil);
                    res = 0;
                    err = 0;
                    for (i=0;i<N;i++)
                    {
                        res += pow(rhs[i],2);
                        err += pow(u_exact[i]-u[i],2);
                    }
                    res = sqrt(res)/res_init;       //relative residual norm
                    err = sqrt(err)/err_init;       //relative error norm
                    solution->err_values.push_back(err);
                    solution->res_values.push_back(res);

                    //check convergence 
                    if(iter==maxiter)
                    {
                        cout << "Maximum number of iterations reached\n";
                        break;
                    }
                    if (res<tol || iter > maxiter)
                        break;

                    //count the next iteration
                    iter ++;
            
                    //solution of the forward system L*w=rhs
                    k=0;
                    w[k] = rhs[k]/LU[k][2];
                    for (k=1;k<jmax;k++)
                        w[k] = (rhs[k]-LU[k][1]*w[k-1])/LU[k][2];
                    for (k=jmax;k<N;k++)
                        w[k] = (rhs[k]-LU[k][0]*w[k-jmax]-LU[k][1]*w[k-1])/LU[k][2];   
                    //solution of the backward system U*du=w
                    k=N-1;
                    du[k] = w[k];
                    for (k=N-2;k>=N-jmax-2;k--)
                        du[k] = w[k] - LU[k][3]*du[k+1];
                    for (k=N-jmax-1;k>=0;k--)
                        du[k] = w[k] - LU[k][3]*du[k+1] - LU[k][4]*du[k+jmax];

                    //new u 
                    for (k=1;k<N;k++)
                        u[k] = u[k] + omega*du[k];   
                }
                break;
    
            default:
                while(true)
                {   
                    //calculation of the new rhs and the relative residual
                    rhs_constructor(rhs,A,b,u,N,jmax,stencil);
                    res = 0;
                    err = 0;
                    for (i=0;i<N;i++)
                    {
                        res += pow(rhs[i],2);
                        err += pow(u_exact[i]-u[i],2);
                    }
                    res = sqrt(res)/res_init;       //relative residual norm
                    err = sqrt(err)/err_init;       //relative error norm
                    //save iteration , residual and error 
                    solution->err_values.push_back(err);
                    solution->res_values.push_back(res);

                    //check convergence 
                    if(iter==maxiter)
                    {
                        cout << "Maximum number of iterations reached\n";
                        break;
                    }
                    if (res<tol || iter > maxiter)
                        break;

                    //count the next iteration
                    iter++;
            
                    //solution of the forward system L*w=rhs
                    k=0;
                    w[k] = rhs[k]/LU[k][4];
                    for (k=1;k<jmax-1;k++)
                        w[k] = (rhs[k] - LU[k][3]*w[k-1])/LU[k][4];
                    w[k] = (rhs[k] - LU[k][2]*w[k-jmax+1] - LU[k][3]*w[k-1])/LU[k][4];    
                    k = jmax;
                    w[k] = (rhs[k] - LU[k][1]*w[k-jmax] - LU[k][2]*w[k-jmax+1] - LU[k][3]*w[k-1])/LU[k][4];   
                    for (k=jmax+1;k<N;k++)
                        w[k] = (rhs[k] - LU[k][0]*w[k-jmax-1] - LU[k][1]*w[k-jmax] - LU[k][2]*w[k-jmax+1] - LU[k][3]*w[k-1])/LU[k][4];
                    //solution of the backward system U*du=w
                    k=N-1;
                    du[k] = w[k];
                    for (k=N-2;k>=N-jmax+1;k--)
                        du[k] = w[k] - LU[k][5]*du[k+1];
                    k = N-jmax;
                    du[k] = w[k] - LU[k][5]*du[k+1] - LU[k][6]*du[k+jmax-1];
                    k=N-jmax-1;
                    du[k] = w[k] - LU[k][5]*du[k+1] - LU[k][6]*du[k+jmax-1] - LU[k][7]*du[k+jmax];
                    for (k=N-jmax-2;k>=0;k--)
                        du[k] = w[k] - LU[k][5]*du[k+1] - LU[k][6]*du[k+jmax-1] - LU[k][7]*du[k+jmax] - LU[k][8]*du[k+jmax+1];

                    //new u 
                    for (k=1;k<N;k++)
                        u[k] = u[k] + omega*du[k];  
                }    
                break;
            }
            //----------Results--------------
        	// Console.WriteLine($"System solved\nFinal relative residual norm: {res}\nFinal relative error norm: {err}\nNumber of iterations: {iter}");
        	cout << "System solved\nFinal relative residual norm: " << res <<"\nFinal relative error norm: " << err << "\nNumber of iterations: " << iter << "\n";
        }

    public: 
        /// passes the construction of the A, b, u matrix and vectors to a delegate
        // public void Construction(Action<Solution> action)
        void Construction(void (*action)(Solution*)) {
            if (is_constructed) return;
            action(solution);
            is_constructed = true;
        }

        void Decomposition() {
            // if (!is_initialized) throw new Exception("matrices are not initialized");
            if (is_decomposed) return;
            decomposition();
            is_decomposed = true;
        }
        
        void InitialSolution() {
            // if (!is_constructed) throw new Exception("matrices are not constructed");
            if (is_initialized) return;
            initialSolution();
            is_initialized = true;
        }

        void InitialSolution(void (*action)(Solution*)) {
            // if (!is_constructed) throw new Exception("matrices are not constructed");
            if (is_initialized) return;        
            action(solution);
            is_initialized = true;
        }

        void Solve() {
            // if (!is_decomposed) throw new Exception("matrices are not decomposed");
            if (is_solved) return;
            solve();
            is_solved = true;
        }        

        
        void SerializeU() {
            ofstream residuals;
            residuals.open("Residuals.dat");
            for (int i = 0; i < solution->res_values.size(); i++) {
                auto res = to_string(solution->res_values[i]) + "\n";                
                residuals.write((char*)&res, sizeof(string));
            }
            residuals.close();

            auto u = solution->u;
            auto u_exact = solution->u_exact;
            ofstream xyu;
            xyu.open("u.dat");
            for (j=jmax-1;j>=0;j--)
            {
                for (i=0;i<imax;i++)
                {
                    k = i*jmax+j;
        			auto str = to_string(i*dx) + "  " + to_string(j*dy) + "  " + to_string(u[k]) + "  " + to_string(u_exact[k]) + "\n";
        			xyu.write((char*)&str, sizeof(string));
                }
            }
            xyu.close();        
        }
};
