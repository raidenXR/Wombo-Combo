#r "bin/Debug/net9.0/msip.dll"

open MSIP
open System
open System.IO

// analytical solutions
let pi = System.Math.PI
let K = 2.
let inline u'(x,y) = sin(pi*K*x) * sin(pi*K*y)
let inline g'(x,y) = -2. * pi**2 * sin(pi*K*x) * sin(pi*K*y) - sin(pi*K*x)**3. * sin(pi*K*y)**3.
let pow(x,y) = pown x y


// set the parameters for the MSIP algorithm
let pars = Params(
    imax = 100,
    jmax = 100,
    tol = 1e-6,
    psi = 1,
    omega = 1.0,
    maxiter = 1000
)

// instanciate an object for the MSIP algorith with 9-scheme
let scheme = MSIP_N.MSIP5
let solver = Solver(pars, scheme)
let stencil = match scheme with | MSIP_N.MSIP5 -> 5 | MSIP_N.MSIP9 -> 9
let N = (pars.imax + 1) * (pars.jmax + 1)
let A_n1 = Array.init N (fun x -> Array.zeroCreate<double> stencil)  // to cache previous iterations
let A_n2 = Array.init N (fun x -> Array.zeroCreate<double> stencil)

// copy current iteration to older
let copy (n:array<array<double>>) (n_old:array<array<float>>) =
    for i in 0..n.Length-1 do
        for j in 0..n[i].Length-1 do
            n_old[i][j] <- n[i][j]

let set (n:array<array<double>>) dx dy f =
    let dx = solver.Dx
    let dy = solver.Dy
    for i in 0..n.Length-1 do
        for j in 0..n[i].Length-1 do
            let x = float(i)*dx
            let y = float(j)*dy
            n[i][j] <- f(x,y)
            

// construct the matrices for the solution 
solver.Construction((fun s ->
    let A = s.A
    let b = s.b
    let u_exact = s.u_exact
    let jmax = solver.Jmax
    let imax = solver.Imax
    let dx = solver.Dx
    let dy = solver.Dy
    let mutable i = 0
    let mutable j = 0
    let mutable k = 0

    set A_n1 dx dx u'
    set A_n2 dx dx u'

    match scheme with
    | MSIP_N.MSIP5 ->
        //boundary conditions , west and east boundaries 
        for i=0 to imax-1 do
            let x = float(i)*dx
            let y = float(j)*dy
            j<-0
            k<-i*jmax+j
            A[k][2] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)       //exact solution
            j<-jmax-1
            k<-i*jmax+j
            A[k][2] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)
        //boundary conditions , north and south boundaries
        for j=1 to jmax-2 do
            let x = float(i)*dx
            let y = float(j)*dy
            i<-0
            k<-i*jmax+j
            A[k][2] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)
            i<-imax-1
            k<-i*jmax+j
            A[k][2] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)
        //internal nodes
        for i=1 to imax-2 do
            for j=1 to jmax-2 do
                let x = float(i)*dx
                let y = float(j)*dy
                k <- i*jmax+j
                A[k][0] <- 1./dx**2
                A[k][1] <- 1./dy**2
                A[k][2] <- -2./dx**2 - 2./dy**2
                A[k][3] <- 1./dy**2
                A[k][4] <- 1./dx**2
                b[k] <- g'(x,y) + A_n1[k][2] * A_n2[k][2] * A_n2[k][2]
                u_exact[k] <- u'(x,y)

        copy A A_n1
        copy A_n1 A_n2
    | MSIP_N.MSIP9 -> 
        //boundary conditions , west and east boundaries 
        for i=0 to imax-1 do
            let x = float(i)*dx
            let y = float(j)*dy
            j<-0
            k<-i*jmax+j
            A[k][4] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)  // exact solution
            j<-jmax-1
            k<-i*jmax+j
            A[k][4] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)
        //boundary conditions , north and south boundaries
        for j=1 to jmax-2 do
            let x = float(i)*dx
            let y = float(j)*dy
            i<-0
            k<-i*jmax+j
            A[k][4] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)
            i<-imax-1
            k<-i*jmax+j
            A[k][4] <- 1.0
            b[k] <- 0.0
            u_exact[k] <- u'(x,y)
            //internal nodes
        for i=1 to imax-2 do
            for j=1 to jmax-2 do
                let x = float(i)*dx
                let y = float(j)*dy
                k <- i*jmax+j
                A[k][1] <- 1./dx**2
                A[k][3] <- 1./dy**2
                A[k][4] <- -2./dx**2 - 2./dy**2
                A[k][5] <- 1./dy**2
                A[k][7] <- 1./dx**2
                b[k] <- g'(x,y) + A_n1[k][4] * A_n2[k][4] * A_n2[k][4]
                u_exact[k] <- u'(x,y)

        copy A A_n1
        copy A_n1 A_n2
))

// set an initial estimation
solver.InitialSolution((fun s ->
    let u = s.u
    let du = s.du
    let w = s.w
    let rhs = s.rhs
    let N = s.N

    for i=0 to N-1 do
        u[i] <- 0.
        du[i] <- 0.
        w[i] <- 0.
        rhs[i] <- 0.
))

// A -> LIU decomposition
solver.Decomposition()
solver.Solve()
solver.SerializeU()
