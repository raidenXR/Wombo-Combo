namespace MSIP;
using System;
using System.IO;

public enum MSIP_N
{
    MSIP5,
    MSIP9
}

public struct Params
{
    public int imax;
    public int jmax;
    public double tol;
    public double psi;
    public double omega;
    public int maxiter;
}

public class Solution
{
    public double[][] A {get;set;}
    public double[][] LU {get;set;}
    public double[] b {get;set;}
    public double[] du {get;set;}
    public double[] u {get;set;}
    public double[] u_exact {get;set;}
    public double[] w {get;set;}
    public double[] rhs {get;set;}
    public readonly List<double> err_values  = new(1000);
    public readonly List<double> res_values = new(1000);

    public int N {get;private set;}

    public Solution(int N, MSIP_N msip)
    {
        var stencil = msip switch {MSIP_N.MSIP5 => 5, MSIP_N.MSIP9 => 9};
        this.N = N;
        b = new double [N];
        du = new double [N];
        u_exact = new double [N];
        u = new double [N];
        w = new double [N];
        rhs = new double [N];
        A = new double[N][];
        LU = new double[N][];
        for (int i=0;i<N;i++)
            LU[i] = new double[stencil];   

        for (int i=0;i<N;i++)
            A[i] = new double[stencil];        
            
        //----------Construction of system matrices and vectors----------         
        //Initiation of matrices 
        for (int i=0;i<N;i++)
        {
            for (int j=0;j<stencil;j++)
            {
                A[i][j] = 0;
            }
            b[i] = 0;
            u_exact[i] = 0;    
        }
    }
}

public class Solver
{
    int solver,imax,jmax,N,maxiter,stencil,i,j,k,write,iter,selection,extra_iter;
    double psi,tol,omega,dx,dy,res_init,res,err,err_init;

    // double[][] A, LU;
    // double[] b, du, u_exact, u, w, rhs;
    Solution solution;

    bool is_constructed, is_initialized, is_decomposed, is_solved;


    double pow(double x, int y)
    {
    	return Math.Pow(x,y);
    }

    double sqrt (double x)
    {
    	return Math.Sqrt(x);
    }

    public Solver(Params pars, MSIP_N msip)
    {

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
        stencil = msip switch {MSIP_N.MSIP5 => 5, MSIP_N.MSIP9 => 9};

    }

    //Construction of matrix A and vector b - CHANGE IT IF YOU WANT TO SOLVE A DIFFERENT PROBLEM
    void construction()
    {
        var A = solution.A;
        var b = solution.b;
        var u_exact = solution.u_exact;
        
        switch (stencil)
        {
        	case 5:
        	    //boundary conditions , west and east boundaries 
        	    for (i=0;i<imax;i++)
        	    {
        	        j=0;
        	        k=i*jmax+j;
        	        A[k][2] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);       //exact solution
        	        j=jmax-1;
        	        k=i*jmax+j;
        	        A[k][2] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	    }
        	    //boundary conditions , north and south boundaries
        	    for (j=1;j<jmax-1;j++)
        	    {
        	        i=0;
        	        k=i*jmax+j;
        	        A[k][2] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	        i=imax-1;
        	        k=i*jmax+j;
        	        A[k][2] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	    }
        	    //internal nodes
        	    for (i=1;i<imax-1;i++)
        	    {
        	        for (j=1;j<jmax-1;j++)
        	        {
        	            k = i*jmax+j;
        	            A[k][0] = 1/pow(dx,2);
        	            A[k][1] = 1/pow(dy,2);
        	            A[k][2] = -2/pow(dx,2)-2/pow(dy,2);
        	            A[k][3] = 1/pow(dy,2);
        	            A[k][4] = 1/pow(dx,2);
        	            b[k] = 4;
        	            u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	        }
        	    }
        	    break;

        	default:
        	    //boundary conditions , west and east boundaries 
        	    for (i=0;i<imax;i++)
        	    {
        	        j=0;
        	        k=i*jmax+j;
        	        A[k][4] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);       //exact solution
        	        j=jmax-1;
        	        k=i*jmax+j;
        	        A[k][4] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	    }
        	    //boundary conditions , north and south boundaries
        	    for (j=1;j<jmax-1;j++)
        	    {
        	        i=0;
        	        k=i*jmax+j;
        	        A[k][4] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	        i=imax-1;
        	        k=i*jmax+j;
        	        A[k][4] = 1.0;
        	        b[k] = pow(i*dx,2)+pow(j*dy,2);
        	        u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	    }
        	    //internal nodes
        	    for (i=1;i<imax-1;i++)
        	    {
        	        for (j=1;j<jmax-1;j++)
        	        {
        	            k = i*jmax+j;
        	            A[k][1] = 1/pow(dx,2);
        	            A[k][3] = 1/pow(dy,2);
        	            A[k][4] = -2/pow(dx,2)-2/pow(dy,2);
        	            A[k][5] = 1/pow(dy,2);
        	            A[k][7] = 1/pow(dx,2);
        	            b[k] = 4;
        	            u_exact[k] = pow(i*dx,2)+pow(j*dy,2);
        	        }
        	    }
        	    break;
        }
    }

    public int Imax => imax;
    public int Jmax => jmax;
    public double Dx => dx;
    public double Dy => dy;

    /// uses the default constrution of the problem
    public void Construction()
    {
        if (is_constructed) return;
        construction();
        is_constructed = true;
    }

    /// passes the construction of the A, b, u matrix and vectors to a delegate
    public void Construction(Action<Solution> action)
    {
        if (is_constructed) return;
        action(solution);
        is_constructed = true;
    }

    void decomposition()
    {
        var A = solution.A;
        var LU = solution.LU;
        
        //Incomplete LU decomposition
        switch (stencil)
        {
        	case 5:
        	    for (i=0;i<N;i++)
        	    {
        	        for (j=0;j<stencil;j++)
        	        {
        	        LU[i][j] = A[i][j];
        	        } 
        	    }
        	    msip5(LU,psi,N,jmax); 
        	    break;

        	default:
        		for (i=0;i<N;i++)
        		    {
        		        for (j=0;j<stencil;j++)
        		        {
        		        LU[i][j] = 0;
        		        } 
        		    }
        		msip9(LU,A,psi,N,jmax);
        	    break;
        }           
    }

    public void Decomposition()
    {
        if (!is_initialized) throw new Exception("matrices are not initialized");
        if (is_decomposed) return;
        decomposition();
        is_decomposed = true;
    }

    void initialSolution()
    {
        var u = solution.u;
        var du = solution.du;
        var w = solution.w;
        var rhs = solution.rhs;
        
        //Initial solution 
        for (i=0;i<N;i++)
        {
            u[i] = 0;
            du[i] = 0;
            w[i] = 0;
            rhs[i] =0;
        }        
    }

    public void InitialSolution()
    {
        if (!is_constructed) throw new Exception("matrices are not constructed");
        if (is_initialized) return;
        initialSolution();
        is_initialized = true;
    }

    public void InitialSolution(Action<Solution> action)
    {
        if (!is_constructed) throw new Exception("matrices are not constructed");
        if (is_initialized) return;        
        action(solution);
        is_initialized = true;
    }

    void solve()
    {
        var A = solution.A;
        var LU = solution.LU;
        var b = solution.b;
        var u = solution.u;
        var w = solution.w;
        var du = solution.du;
        var rhs = solution.rhs;
        var u_exact = solution.u_exact;
        
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
                solution.err_values.Add(err);
                solution.res_values.Add(res);

                //check convergence 
                if(iter==maxiter)
                {
                    Console.WriteLine("Maximum number of iterations reached");
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
                solution.err_values.Add(err);
                solution.res_values.Add(res);

                //check convergence 
                if(iter==maxiter)
                {
                    Console.WriteLine("Maximum number of iterations reached");
                    break;
                }
                if (res<tol || iter > maxiter)
                    break;

                //count the next iteration
                iter ++;
            
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
    	Console.WriteLine($"System solved\nFinal relative residual norm: {res}\nFinal relative error norm: {err}\nNumber of iterations: {iter}");
    }

    public void SaveOutput()
    {
        var u = solution.u;
        var u_exact = solution.u_exact;
        var A = solution.A;
        var LU = solution.LU;
        var b = solution.b;
    
    	if (!Directory.Exists("results")) Directory.CreateDirectory("results");
    	var fs = File.CreateText("results/A.dat");

    	for (i=0; i<N; i++) {
    		for (j=0; j<stencil; j++) {
    			fs.Write($"{A[i][j]}".PadLeft(10));
    		}
    		fs.Write("\n");
    	}
    	fs.Close();

    	fs = File.CreateText("results/b.dat");
    	for (i=0; i<N; i++) {
    		fs.WriteLine(b[i]);
    	}
    	fs.Close();

    	fs = File.CreateText("results/LU");
    	for (i=0; i<N; i++) {
    		for (j=0; j<stencil; j++) {
    			fs.Write($"{LU[i][j]}".PadLeft(10));
    		}
    		fs.Write("\n");
    	}
    	fs.Close();

        //----------Output-----------
        var x_values = new List<double>(N);
        var y_values = new List<double>(N);
        var u_values = new List<double>(N);
        var u_exact_values = new List<double>(N);
        
        var X = File.CreateText("results/X.dat");
        for (j=0;j<jmax;j++)
        {
            for (i=0;i<imax;i++)
    			X.Write($"{i*dx}".PadLeft(10));
                x_values.Add(i*dx);
            X.Write("\n");  
        }
        X.Close();

        var Y = File.CreateText("results/Y.dat");
        for (j=jmax-1;j>=0;j--)
        {
        	for (i=0;i<imax;i++)
    			Y.Write($"{j*dy}".PadLeft(10));
                y_values.Add(j*dy);
            Y.Write("\n");  
    	}
    	Y.Close();

        var U = File.CreateText("results/U.dat");
        for (j=jmax-1;j>=0;j--)
        {
            for (i=0;i<imax;i++)
            {
                k = i*jmax+j;
    			U.Write($"{u[k]}".PadLeft(10));
                u_values.Add(u[k]);
            }
            U.Write("\n");  
        }
    	U.Close();

        var U_exact = File.CreateText("results/U_exact.dat");
        for (j=jmax-1;j>=0;j--)
        {
            for (i=0;i<imax;i++)
            {
                k = i*jmax+j;
    			U_exact.Write($"{u_exact[k]}".PadLeft(10));
                u_exact_values.Add(u_exact[k]);
            }
            U_exact.Write("\n");  
        }
    	U_exact.Close();

        var residuals = File.CreateText("results/Residuals.dat");
        foreach (var res in solution.res_values)
        {
            residuals.WriteLine(res);
        }
        residuals.Close();
    }

    public void SerializeU()
    {
    	if (!Directory.Exists("results")) Directory.CreateDirectory("results");
        var residuals = File.CreateText("results/Residuals.dat");
        foreach (var res in solution.res_values)
        {
            residuals.WriteLine(res);
        }
        residuals.Close();

        var u = solution.u;
        var u_exact = solution.u_exact;
        var xyu = File.CreateText("results/u.dat");
        for (j=jmax-1;j>=0;j--)
        {
            for (i=0;i<imax;i++)
            {
                k = i*jmax+j;
    			xyu.WriteLine($"{i*dx}   {j*dy}   {u[k]}   {u_exact[k]}");
            }
        }
        xyu.Close();        
    }
    
    [Obsolete]
    public void SaveMatrices()
    {
        if (is_solved) return;
        var A = solution.A;
        var LU = solution.LU;
        var b = solution.b;
        var u = solution.u;
        var w = solution.w;
        var du = solution.du;
        var rhs = solution.rhs;
        var u_exact = solution.u_exact;
    
    	if (!Directory.Exists("results")) Directory.CreateDirectory("results");
    	var fs = File.CreateText("results/A.dat");

    	for (i=0; i<N; i++) {
    		for (j=0; j<stencil; j++) {
    			fs.Write($"{A[i][j]}".PadLeft(10));
    		}
    		fs.Write("\n");
    	}
    	fs.Close();

    	fs = File.CreateText("results/b.dat");
    	for (i=0; i<N; i++) {
    		fs.WriteLine(b[i]);
    	}
    	fs.Close();

    	fs = File.CreateText("results/LU");
    	for (i=0; i<N; i++) {
    		for (j=0; j<stencil; j++) {
    			fs.Write($"{LU[i][j]}".PadLeft(10));
    		}
    		fs.Write("\n");
    	}
    	fs.Close();

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
	
    	var cp = File.CreateText("results/convergence.dat");

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
                //save iteration , residual and error 
    			cp.WriteLine($"iter: {iter}, res: {res}, err: {err}");

                //check convergence 
                if (res<tol)
                    break;

                if(iter==maxiter)
                {
                    Console.WriteLine("Maximum number of iterations reached");
                    break;
                //     cout<<"Maximum number of iterations reached"<<endl<<"Current relative residual norm: "<<res<<endl;
                //     cout<<"Select one of the following: "<<endl<<"1. More iterations"<<endl<<"2. End"<<endl<<"Your selection: ";
                //     cin>>selection;
                //     if (selection == 1)
                //     {
                //         cout<<"Select the number of extra iterations: ";
                //         cin>>extra_iter;
                //         maxiter += extra_iter;
                //     }
                //     else 
                //         break;
                }    

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
    			cp.WriteLine($"iter: {iter}, res: {res}, err: {err}");
                //check convergence 
                if (res<tol)
                    break;

                if(iter==maxiter)
                {
                    Console.WriteLine("Maximum number of iterations reached");
                    break;
                //     cout<<"Maximum number of iterations reached"<<endl<<"Current relative residual norm: "<<res<<endl;
                //     cout<<"Select one of the following: "<<endl<<"1. More iterations"<<endl<<"2. End"<<endl<<"Your selection: ";
                //     cin>>selection;
                //     if (selection == 1)
                //     {
                //         cout<<"Select the number of extra iterations: ";
                //         cin>>extra_iter;
                //         maxiter += extra_iter;
                //     }
                //     else 
                //         break;
                }    

                //count the next iteration
                iter ++;
            
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
        cp.Close();

        //----------Results--------------
    	Console.WriteLine($"System solved\nFinal relative residual norm: {res}\nFinal relative error norm: {err}\nNumber of iterations: {iter}");

        //----------Output-----------
        var X = File.CreateText("results/X.dat");
        for (j=0;j<jmax;j++)
        {
            for (i=0;i<imax;i++)
    			X.Write($"{i*dx}".PadLeft(10));
            X.Write("\n");  
        }
        X.Close();

        var Y = File.CreateText("results/Y.dat");
        for (j=jmax-1;j>=0;j--)
        {
        	for (i=0;i<imax;i++)
    			Y.Write($"{j*dy}".PadLeft(10));
            Y.Write("\n");  
    	}
    	Y.Close();

        var U = File.CreateText("results/U.dat");
        for (j=jmax-1;j>=0;j--)
        {
            for (i=0;i<imax;i++)
            {
                k = i*jmax+j;
    			U.Write($"{u[k]}".PadLeft(10));
            }
            U.Write("\n");  
        }
    	U.Close();

        var U_exact = File.CreateText("results/U_exact.dat");
        for (j=jmax-1;j>=0;j--)
        {
            for (i=0;i<imax;i++)
            {
                k = i*jmax+j;
    			U_exact.Write($"{u_exact[k]}".PadLeft(10));
            }
            U_exact.Write("\n");  
        }
    	U_exact.Close();
        is_solved = true;
    }

    public void Solve()
    {
        if (!is_decomposed) throw new Exception("matrices are not decomposed");
        if (is_solved) return;
        solve();
        is_solved = true;
    }

    public Solution Solution => solution;
    
    void msip5(double[][]  A, double psi, int N, int jmax)
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

    
    void msip9(double[][] LU, double[][] A,double psi, int N, int jmax)
    {
        int k=0;
        LU[k][4] = A[k][4];
        LU[k][5] = A[k][5]/LU[k][4];
        LU[k][6] = A[k][6]/LU[k][4];
        LU[k][7] = A[k][7]/LU[k][4];
        LU[k][8] = A[k][8]/LU[k][4];
        for (k=1;k<jmax-1;k++)
        {
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
        for (k=jmax+1;k<N-jmax-1;k++)
        {
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
        for (k=N-jmax+1;k<N-1;k++)
        {
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

    void rhs_constructor(double[] rhs, double[][] A, double[] b, double[] u, int N, int jmax, int stencil)
    {
        int k=0;
        switch (stencil)
        {
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
}








