#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open System
open GridGeneration
open Plotting

let out_nodes = "../dat-files/correct_grid_4.dat"
let out_nodes_2 = "../dat-files/correct_grid_compare_4.dat"
let out_boundaries = "../dat-files/boundaries_4.dat"

let N = 30
let mutable D = Vec2(20., 50.)
let mutable A = Vec2(10., 10.)
let mutable C = Vec2(60., 40.)
let mutable B = Vec2(80., 0.)

let x_range (sides_a:Sides) (sides_b:Sides) N =
    let mutable x_min = sides_a.U[0].X
    let mutable x_max = sides_a.U[0].X
    
    for (a,b) in [(sides_a.D, sides_b.D); (sides_a.R, sides_b.R); (sides_a.L, sides_b.L); (sides_a.U, sides_b.U)] do
        for i in 0..N-1 do
            x_min <- if a[i].X < x_min then a[i].X else x_min
            x_min <- if b[i].X < x_min then b[i].X else x_min
            x_max <- if a[i].X > x_max then a[i].X else x_max
            x_max <- if b[i].X > x_max then b[i].X else x_max
    $"set xrange[{x_min}:{x_max}]"

let y_range (sides_a:Sides) (sides_b:Sides) N =
    let mutable y_min = sides_a.U[0].Y
    let mutable y_max = sides_a.U[0].Y
    
    for (a,b) in [(sides_a.D, sides_b.D); (sides_a.R, sides_b.R); (sides_a.L, sides_b.L); (sides_a.U, sides_b.U)] do
        for i in 0..N-1 do
            y_min <- if a[i].Y < y_min then a[i].Y else y_min
            y_min <- if b[i].Y < y_min then b[i].Y else y_min
            y_max <- if a[i].Y > y_max then a[i].Y else y_max
            y_max <- if b[i].Y > y_max then b[i].Y else y_max
    $"set yrange[{y_min}:{y_max}]"

let create_map (mesh:Vec2[,]) N (table:double[,]) (tag:string) =
    let sb = new System.Text.StringBuilder(1000)
    sb.AppendLine("$" + tag + " << EOD") |> ignore

    for i in 0..N-1 do
        for j in 0..N-1 do
            let d = table[i,j]
            let x = mesh[i,j].X
            let y = mesh[i,j].Y
            sb.AppendLine($"{x} {y} {d}").Append(" ") |> ignore
        sb.AppendLine("\n") |> ignore
    sb.AppendLine("EOD").ToString()
    

let compute_signal_intensity (mesh:Vec2[,]) (antennas_pos:Vec2[]) (signal:double[,]) = 
    // clear previous signal buffer
    for i in 0..N-1 do
        for j in 0..N-1 do
            signal[i,j] <- 0.0

    // compute signal intensity
    for k in 0..5 do
        let dv_x = (B-A).X
        let dv_y = (D-A).Y
        let x0 = antennas_pos[k].X
        let y0 = antennas_pos[k].Y

        for i in 0..N-1 do
            for j in 0..N-1 do
                let A  = 1.      // amplitude
                let x = mesh[i,j].X
                let y = mesh[i,j].Y
                let sx = 0.15 * dv_x
                let sy = 0.15 * dv_y
                signal[i,j] <- signal[i,j] + Signal.Gaussian(A, x, x0, y, y0, sx, sy)

let get_closest_antenna (p:Vec2) (antennas_pos:Vec2[]) =
    let mutable n = 0
    let mutable l = (p - antennas_pos[0]).Length()
    for i in 1..antennas_pos.Length-1 do
        let a = antennas_pos[i]
        if l > (p - a).Length() then
            n <- i
            l <- (p - a).Length()
    n

let _sides =
    let l = Array.zeroCreate<Vec2> N
    let r = Array.zeroCreate<Vec2> N
    let d = Array.zeroCreate<Vec2> N
    let u = Array.zeroCreate<Vec2> N
    for i in 0..N-1 do
        d[i] <- A + ((B-A)/double (N-1)) * (double i)
        u[i] <- D + ((C-D)/double (N-1)) * (double i)

        l[i] <- A + ((D-A)/double (N-1)) * (double i)
        r[i] <- B + ((C-B)/double (N-1)) * (double i)
    GeoRandomizer.RandomSides(l, r, d, u)

let sides_noise =
    let l = Array.zeroCreate<Vec2> N
    let r = Array.zeroCreate<Vec2> N
    let d = Array.zeroCreate<Vec2> N
    let u = Array.zeroCreate<Vec2> N
    for i in 0..N-1 do
        d[i] <- A + ((B-A)/double (N-1)) * (double i)
        u[i] <- D + ((C-D)/double (N-1)) * (double i)

        l[i] <- A + ((D-A)/double (N-1)) * (double i)
        r[i] <- B + ((C-B)/double (N-1)) * (double i)
    GeoRandomizer.RandomSidesWithNoise(l, r, d, u)

// fill the boundaries with non-linear geometry. make it curvy!
let sides_curve =
    let l = Array.zeroCreate<Vec2> N
    let r = Array.zeroCreate<Vec2> N
    let d = Array.zeroCreate<Vec2> N
    let u = Array.zeroCreate<Vec2> N
    for i in 0..N-1 do
        l[i] <- A + ((D-A)/double (N-1)) * (double i)

    for i in 0..N-1 do
        r[i] <- B + ((C-B)/double (N-1)) * (double i)
        r[i].X <- r[i].X + r[i].Y**2 * 0.02
        
    for i in 0..N-1 do
        d[i] <- A + ((B-A)/double (N-1)) * (double i)

    C.X <- r[N-1].X
    for i in 0..N-1 do
        u[i] <- D + ((C-D)/double (N-1)) * (double i)
        u[i].Y <- u[i].Y - 2.0 * cos ((C.Y - D.Y)/double i)

    GeoRandomizer.RandomSides(l, r, d, u)

let sides_curve_noise =
    let l = Array.zeroCreate<Vec2> N
    let r = Array.zeroCreate<Vec2> N
    let d = Array.zeroCreate<Vec2> N
    let u = Array.zeroCreate<Vec2> N
    for i in 0..N-1 do
        l[i] <- A + ((D-A)/double (N-1)) * (double i)

    for i in 0..N-1 do
        r[i] <- B + ((C-B)/double (N-1)) * (double i)
        r[i].X <- r[i].X + r[i].Y**2 * 0.02
        
    for i in 0..N-1 do
        d[i] <- A + ((B-A)/double (N-1)) * (double i)

    C.X <- r[N-1].X
    for i in 0..N-1 do
        u[i] <- D + ((C-D)/double (N-1)) * (double i)
        u[i].Y <- u[i].Y - 4.0 * cos ((C.Y - D.Y)/double i)
    GeoRandomizer.RandomSidesWithNoise(l, r, d, u)


let grids = [
    "sides", _sides
    "sides_noise", sides_noise
    "sides_curve", sides_curve
    "sides_curve_noise", sides_curve_noise
]

let sides_constraint (p:Vec2) (sides:Sides) N =
    let mutable b = true
    let mutable i = 0

    while i < N && b do
        b <- p.X > sides.L[i].X
        i <- i + 1

    i <- 0
    while i < N && b do
        b <- p.X < sides.R[i].X
        i <- i + 1

    i <- 0
    while i < N && b do
        b <- p.Y > sides.D[i].Y
        i <- i + 1
        
    i <- 0
    while i < N && b do
        b <- p.Y < sides.U[i].Y
        i <- i + 1
    b

let distance_constraint (p:int*int) (positions:seq<(int*int)>) =
    let mutable b = Seq.length positions > 0
    let pi = fst p
    let pj = snd p

    for (i,j) in positions do
        b <- b && if (pi <> i || pj <> j) then abs((pi + pj) - (i + j)) >= 8 else b
    b


for (name,sides) in grids do
    let mesh = GeoRandomizer.MeshFromRandomizer(sides, N, N)
    let signal = Array2D.zeroCreate<double> N N
    
    let (antennas_pos, antennas_idx) = 
        let antennas = ResizeArray<Vec2>()
        let mutable i = Random.Shared.Next(1,N-2)
        let mutable j = Random.Shared.Next(1,N-2)
        let mutable cached = [(i,j)]
        antennas.Add(mesh[i,j])

        for n in 0..5 do
            while ((List.contains (i,j) cached) || antennas.Count < 6) do
                i <- Random.Shared.Next(1,N-2)
                j <- Random.Shared.Next(1,N-2)
                if (distance_constraint (i,j) cached) then
                    antennas.Add(mesh[i,j])
                    cached <- (i,j)::cached
                    printfn "antenna added: %d" (antennas.Count)
                
        (Array.ofSeq antennas, Array.ofList cached)

    let gif = System.IO.File.CreateText("output/animation_gif.dat")
    for a in antennas_pos do
        // printfn "%s" (string a)
        gif.WriteLine($"{a.X}  {a.Y}")
    gif.WriteLine("\n")
            
    // Optimization step
    let max_iters = 100
    for n in 1..max_iters do
        compute_signal_intensity mesh antennas_pos signal
        let mutable s = 0.
        // get mean value
        for i in 0..N-1 do
            for j in 0..N-1 do
                s <- s + signal[i,j]
        s <- s / double(N*N)                

        // ************************************
        // ADD CONSTRAINTS
        // ************************************
        for i in 1..N-2 do
            for j in 1..N-2 do
                let p = mesh[i,j]
                let k = get_closest_antenna p antennas_pos 
                if signal[i,j] >  (s * 0.7) then
                    let I,J =
                        match (GeoRandomizer.GetGradient(mesh, signal, i,j).SteepestDescend()) with
                        | GradDirection.LL -> (i-1,j-1)
                        | GradDirection.LC -> (i-1,j+0)
                        | GradDirection.LR -> (i-1,j+1)
                        | GradDirection.UL -> (i+1,j-1)
                        | GradDirection.UC -> (i+1,j+0)
                        | GradDirection.UR -> (i+1,j+1)
                        | GradDirection.CL -> (i,j-1)
                        | GradDirection.CR -> (i,j+1)
                        | _ -> (i,j)
                    let (_i, _j) = (Math.Clamp(I, 1, N-2), Math.Clamp(J,1,N-2))
                    let _p = mesh[_i, _j]
                    if distance_constraint (I,J) antennas_idx && sides_constraint _p sides N then
                        antennas_idx[k] <- (_i,_j)
                        antennas_pos[k] <- mesh[fst antennas_idx[k], snd antennas_idx[k]]
                    
                if signal[i,j] < (s * 0.2) then
                    let I,J = 
                        match (GeoRandomizer.GetGradient(mesh, signal, i,j).SteepestAscend()) with
                        | GradDirection.LL -> (i-1,j-1)
                        | GradDirection.LC -> (i-1,j+0)
                        | GradDirection.LR -> (i-1,j+1)
                        | GradDirection.UL -> (i+1,j-1)
                        | GradDirection.UC -> (i+1,j+0)
                        | GradDirection.UR -> (i+1,j+1)
                        | GradDirection.CL -> (i,j-1)
                        | GradDirection.CR -> (i,j+1)
                        | _ -> (i,j)
                    let (_i, _j) = (Math.Clamp(I, 1, N-2), Math.Clamp(J,1,N-2))
                    let _p = mesh[_i, _j]
                    if distance_constraint (I,J) antennas_idx && sides_constraint _p sides N then
                        antennas_idx[k] <- (_i,_j)
                        antennas_pos[k] <- mesh[fst antennas_idx[k], snd antennas_idx[k]]

        for a in antennas_pos do
            // printfn "%s" (string a)
            gif.WriteLine($"{a.X}  {a.Y}")
        gif.WriteLine("\n")

    let _x_range = if (name.Contains("curve")) then x_range sides_curve sides_curve_noise N else x_range sides sides_noise N
    let _y_range = if (name.Contains("curve")) then y_range sides_curve sides_curve_noise N else y_range sides sides_noise N

    GeoRandomizer.CheckUnitialized(mesh, N, N)
    GeoRandomizer.WriteBoundaryNodes(sides, out_boundaries, N)
    GeoRandomizer.WriteCorrectGrid(out_nodes, mesh, N, N)
        
    let boundary_nodes = GeoRandomizer.SidesToBoundaryNodes(sides, N, N)

    Gnuplot()
    |> Gnuplot.datablockXY (boundary_nodes |> Array.map (fun v -> v.X)) (boundary_nodes |> Array.map (fun v -> v.Y)) "bounds"
    |>> "set terminal png size 840,580"
    |>> $"set output 'output/{name}.png'"
    |>> _x_range
    |>> _y_range
    |>> "unset key"
    |>> $"plot '$bounds' with lines lc rgb 'black' lw 2, \\"
    |>> $"'{out_nodes}' with lines lc rgb 'black'"
    |> Gnuplot.run
    |> ignore

    Gnuplot()
    |> Gnuplot.datablockXY (boundary_nodes |> Array.map (fun v -> v.X)) (boundary_nodes |> Array.map (fun v -> v.Y)) "bounds"
    |>> "set terminal gif animate delay 30"
    |>> $"set output 'output/{name}.gif'"
    |>> _x_range
    |>> _y_range
    |>> "stats 'animation_gif.dat' nooutput"
    |>> "do for [i=1:int(STATS_blocks)] {"
    |>> $"     plot 'output/animation_gif.dat' index (i-1) with circles lw 2, \\"
    |>> "         '$bounds' with lines lc rgb 'black' lw 2, \\"
    |>> $"        '{out_nodes}' with lines lc rgb 'black', \\"
    |>> "}"
    |> Gnuplot.run
    |> ignore

    let _signal =
        let _s = ResizeArray<double>(1000)
        for i in 0..N-1 do
            for j in 0..N-1 do
                _s.Add(signal[i,j])
        _s.ToArray()
    
    Gnuplot()
    |>> "set terminal png size 840,580"
    |>> $"set output 'output/{name}_signal.png'"
    |>> "set palette rgbformulae 33,13,10"
    |>> "set isosample 250,250"
    |>> "set grid"
    |>> "set contour base"
    // |>> "unset surface"
    |>> "set view map"
    |>> "set dgrid3d 30,30 qnorm 4"
    |>> "set cntrparam"
    |> Gnuplot.dataTable (_signal) N N (0.,1.) (0.,1.) "signal" 
    |>> "splot '$signal' w pm3d title 'Signal Intensity'"
    // plt
    // |> Gnuplot.datablockXY (antennas_pos |> Array.map (fun v -> v.X)) (antennas_pos |> Array.map (fun v -> v.Y)) "antennas"
    // |>> "unset key"
    // |>> "set view map"
    // // |>> "set pm3d map impl"
    // // |>> "unset surface"
    // |>> "set contour"
    // |>> "set dgrid3d"
    // |>> "set cntrparam cubicspline"
    // |>> "splot '$signal' with pm3d notitle, '$antennas' using 1:2:(0) with points lc rgb 'black' lw 3"
    |> Gnuplot.run
    |> ignore


    do
        use fs = System.IO.File.CreateText($"output/{name}.csv")
        for a in antennas_pos do
            fs.WriteLine(String.Format("{0:0.000}", a.X) + "," + String.Format("{0:0.000}", a.Y))

