#include <iostream>
#include <vector>
#include <math.h>

#include "mesh.h"
#include "msip.cpp"

using namespace std;

const double pi = M_PI;
const double K = 2.;

double** A_n1;
double** A_n2;

const int imax = 100;
const int jmax = 100;
const double Dx = 1. / imax;
const double Dy = 1. / jmax;
const int N = (imax + 1) * (jmax + 1);
const MSIP_N scheme = MSIP5;
const int stencil = 5;
// int stencil = 0x00;
// switch (scheme) {
//     case MSIP5: stencil = 5; break;
//     case MSIP9: stencil = 9; break; 
// }


inline double u (double x, double y) {
    return sin(pi*K*x) * sin(pi*K*y);
}

inline double g (double x, double y) {
    return -2. * pow(pi,2) * sin(pi*K*x) * sin(pi*K*y) - pow(sin(pi*K*x),3.) * pow(sin(pi*K*y),3.);
}

void m_copy (double** n, double** n_old) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < stencil; j++) {
             n_old[i][j] = n[i][j];       
        }
    }            
}

void m_set (double** n, double dx, double dy, double (*f)(double, double)) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < stencil; j++) {
            double x = i * dx;
            double y = j * dy;
            n[i][j] = f(x,y);       
        }
    }            
}

void initial_solution (Solution* s) {
    for (int i = 0; i < s->N; i++) {
        s->u[i] = 0;
        s->du[i] = 0;
        s->w[i] = 0;
        s->rhs[i] = 0;
    }
}

void construction (Solution* s) {
    auto A = s->A;
    auto b = s->b;
    auto u_exact = s->u_exact;
    auto dx = Dx;
    auto dy = Dy;
    int i = 0;
    int j = 0;
    int k = 0;

    m_set (A_n1, dx, dx, u);
    m_set (A_n2, dx, dx, u);

    switch (scheme) {
        case MSIP5:
            //boundary conditions , west and east boundaries 
            for (i=0; i <= imax-1; i++) {
                double x = i * dx;
                double y = j * dy;
                j=0;
                k=i*jmax+j;
                A[k][2]= 1.0;
                b[k] = 0.0;
                u_exact[k] = u(x,y);       //exact solution
                j=jmax-1;
                k=i*jmax+j;
                A[k][2] = 1.0;
                b[k] = 0.0;
                u_exact[k] = u(x,y);
            }
            //boundary conditions , north and south boundaries
            for (j=1; j<=jmax-2; j++) {
                double x = i*dx;
                double y = j*dy;
                i=0;
                k=i*jmax+j;
                A[k][2] = 1.0;
                b[k] = 0.0;
                u_exact[k] = u(x,y);
                i=imax-1;
                k=i*jmax+j;
                A[k][2] = 1.0;
                b[k] = 0.0;
                u_exact[k] = u(x,y);
            }
            //internal nodes
            for (i=1; i<=imax-2; i++) {
                for (j=1; j <=jmax-2; j++) {
                    double x = i*dx;
                    double y = j*dy;
                    k = i*jmax+j;
                    A[k][0] = 1./pow(dx,2);
                    A[k][1] = 1./pow(dy,2);
                    A[k][2] = -2./pow(dx,2) - 2./pow(dy,2);
                    A[k][3] = 1./pow(dy,2);
                    A[k][4] = 1./pow(dx,2);
                    b[k] = g(x,y) + A_n1[k][2] * A_n2[k][2] * A_n2[k][2];
                    u_exact[k] = u(x,y);
                }
            }

            m_copy (A, A_n1);
            m_copy (A_n1, A_n2);
        break;
    case MSIP9: 
        //boundary conditions , west and east boundaries 
        for (i=0; i<=imax-1; i++) {
            double x = i*dx;
            double y = j*dy;
            j=0;
            k=i*jmax+j;
            A[k][4] = 1.0;
            b[k] = 0.0;
            u_exact[k] = u(x,y);  // exact solution;
            j=jmax-1;
            k=i*jmax+j;
            A[k][4] = 1.0;
            b[k] = 0.0;
            u_exact[k] = u(x,y);
        }
        //boundary conditions , north and south boundaries
        for (j=1; j<=jmax-2; j++) {
            double x = i*dx;
            double y = j*dy;
            i=0;
            k=i*jmax+j;
            A[k][4] = 1.0;
            b[k] = 0.0;
            u_exact[k] = u(x,y);
            i=imax-1;
            k=i*jmax+j;
            A[k][4] = 1.0;
            b[k] = 0.0;
            u_exact[k] = u(x,y);
        }
        //internal nodes
        for (i=1; i <= imax-2; i++) {
            for (j=1; j<=jmax-2; j++) {
                double x = i*dx;
                double y = j*dy;
                k = i*jmax+j;
                A[k][1] = 1./pow(dx,2);
                A[k][3] = 1./pow(dy,2);
                A[k][4] = -2./pow(dx,2) - 2./pow(dy,2);
                A[k][5] = 1./pow(dy,2);
                A[k][7] = 1./pow(dx,2);
                b[k] = g(x,y) + A_n1[k][4] * A_n2[k][4] * A_n2[k][4];
                u_exact[k] = u(x,y);
            }
        }
        m_copy (A, A_n1);
        m_copy (A_n1, A_n2);
        break;
    }    
}

int main ()
{
    auto mesh = Mesh(1000, 3000);
    
    auto pars = Params{
        .imax = imax,
        .jmax = jmax,
        .tol = 1e-6,
        .psi = 1,
        .omega = 1.0,
        .maxiter = 1000,
    };

    auto solver = new Solver(pars, scheme);

    A_n1 = new double*[N];
    A_n2 = new double*[N];
    for (int i = 0; i < N; i++) {
        A_n1[i] = new double[stencil];
        A_n2[i] = new double[stencil];

        for (int j = 0; j < stencil; j++) {
            A_n1[i][j] = 0;
            A_n2[i][j] = 0;
        }
    }
    
    solver->Construction(construction);
    solver->InitialSolution(initial_solution);
    solver->Decomposition();
    solver->Solve();
    solver->SerializeU();

    return 0; 
}

