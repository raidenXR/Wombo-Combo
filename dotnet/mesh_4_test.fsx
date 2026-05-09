#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open System
open GridGeneration
open Plotting

let out_nodes = "../dat-files/correct_grid_4.dat"
let out_nodes_2 = "../dat-files/correct_grid_compare_4.dat"
let out_boundaries = "../dat-files/boundaries_4.dat"

let N = 30
let D = Vec2(20., 50.)
let A = Vec2(10., 10.)
let C = Vec2(60., 40.)
let B = Vec2(80., 0.)


let l = Array.zeroCreate<Vec2> N
let r = Array.zeroCreate<Vec2> N
let d = Array.zeroCreate<Vec2> N
let u = Array.zeroCreate<Vec2> N

// fill the boundaries with non-linear geometry. make it curvy!
// for i in 0..N-1 do
//     d[i] <- A + ((B-A)/double (N-1)) * (double (i) + 14.3)
//     u[i] <- D + ((C-D)/double (N-1)) * (double i * 1.7)

//     l[i] <- A + ((D-A)/double (N-1)) * (double i)
//     r[i] <- B + ((C-B)/double (N-1)) * (double i * double i * 0.4)

for i in 0..N-1 do
    d[i] <- A + ((B-A)/double (N-1)) * (double i)
    u[i] <- D + ((C-D)/double (N-1)) * (double i)

    l[i] <- A + ((D-A)/double (N-1)) * (double i)
    r[i] <- B + ((C-B)/double (N-1)) * (double i)

// let sides = GeoRandomizer.RandomSides(l, r, d, u)
let sides = GeoRandomizer.RandomSidesWithNoise(l, r, d, u)
let mesh = GeoRandomizer.MeshFromRandomizer(sides, N, N)
let signal = Array2D.zeroCreate<double> N N

let (antennas_pos, antennas_idx) = 
    let antennas = ResizeArray<Vec2>()
    let mutable i = Random.Shared.Next(1,N-2)
    let mutable j = Random.Shared.Next(1,N-2)
    let mutable cached = [(i,j)]
    antennas.Add(mesh[i,j])

    for n in 0..5 do
        while ((List.contains (i,j) cached) && antennas.Count < 6) do
            i <- Random.Shared.Next(1,N-2)
            j <- Random.Shared.Next(1,N-2)
        antennas.Add(mesh[i,j])
        cached <- (i,j)::cached
                
    (Array.ofSeq antennas, Array.ofList cached)


let create_map N (table:double[,]) (tag:string) =
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
    

let compute_signal_intensity () = 
    // compute signal intensity
    for k in 0..5 do
        let dv_x = (B-A).X
        let dv_y = (C-A).Y
        let x0 = antennas_pos[k].X
        let y0 = antennas_pos[k].Y

        for i in 0..N-1 do
            for j in 0..N-1 do
                let A  = 1.      // amplitude
                let x = mesh[i,j].X
                let y = mesh[i,j].Y
                let sx = 0.15 * dv_x
                let sy = 0.15 * dv_y
                signal[i,j] <- Signal.Gaussian(A, x, x0, y, y0, sx, sy)


// Optimization step
let max_iters = 30
for n in 1..max_iters do
    compute_signal_intensity()
    let mutable s = 0.
    for k in 0..antennas_idx.Length-1 do
        let (i,j) = antennas_idx[k]
        s <- s + signal[i,j]
        s <- s / 6.

    for k in 0..antennas_idx.Length-1 do
        let (i,j) = antennas_idx[k]
        // printfn "%d, %d" i j
        if signal[i,j] > 1.2 * s then
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
            antennas_idx[k] <- (Math.Clamp(I, 1, N-2), Math.Clamp(J,1,N-2))
            antennas_pos[k] <- mesh[fst antennas_idx[k], snd antennas_idx[k]]
        if signal[i,j] < 0.4 * s then
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
            antennas_idx[k] <- (Math.Clamp(I, 1, N-2), Math.Clamp(J,1,N-2))
            antennas_pos[k] <- mesh[fst antennas_idx[k], snd antennas_idx[k]]
        

for k in 0..5 do
    printfn "%s" (string antennas_pos[k])
            

GeoRandomizer.CheckUnitialized(mesh, N, N)

GeoRandomizer.WriteBoundaryNodes(sides, out_boundaries, N);

GeoRandomizer.WriteCorrectGrid(out_nodes, mesh, N, N)

// let bnodes = [|A; B; C; D; A|]
let boundary_nodes = GeoRandomizer.SidesToBoundaryNodes(sides, N, N)
// let bnodes = GeoRandomizer.BoundsFromMesh(mesh, N, N)

Gnuplot()
|> Gnuplot.datablockXY (boundary_nodes |> Array.map (fun v -> v.X)) (boundary_nodes |> Array.map (fun v -> v.Y)) "bounds"
|>> "set terminal png size 840,580"
|>> "set output 'correct_grid_compare_4.png'"
|>> "unset key"
// |>> $"plot '$bounds' with lines lc rgb 'black' lw 2"
|>> $"plot '$bounds' with lines lc rgb 'black' lw 2, \\"
|>> $"'{out_nodes}' with lines lc rgb 'black'"
// |>> $"plot '{out_nodes}' with lines lc rgb 'black'"
|> Gnuplot.run
|> ignore

let tabular_data = create_map N signal "signal"
// printfn "%s" tabular_data

let plt = Gnuplot()
plt.writeln (tabular_data)

plt
|> Gnuplot.datablockXY (antennas_pos |> Array.map (fun v -> v.X)) (antennas_pos |> Array.map (fun v -> v.Y)) "antennas"
|>> "set terminal png size 840,580"
|>> "set output 'correct_contour_compare_4.png'"
|>> "unset key"
// |>> "set isosample 250,250"
|>> "set view map"
|>> "set palette rgbformulae 33,13,10"
// |>> "set pm3d map impl"
// |>> "unset surface"
|>> "set contour"
|>> "set dgrid3d"
|>> "set cntrparam cubicspline"
|>> "splot '$signal' with pm3d notitle, '$antennas' using 1:2:(0) with points lc rgb 'black' lw 3"
|> Gnuplot.run
|> ignore


