#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open System
open GridGeneration
open Plotting

let out_nodes = "../dat-files/correct_grid_4.dat"
let out_nodes_2 = "../dat-files/correct_grid_compare_4.dat"
let out_boundaries = "../dat-files/boundaries_4.dat"

let N = 21
let D = Vec2(20., 50.)
let A = Vec2(10., 10.)
let C = Vec2(30., 40.)
let B = Vec2(30., 20.)


let l = Array.zeroCreate<Vec2> N
let r = Array.zeroCreate<Vec2> N
let d = Array.zeroCreate<Vec2> N
let u = Array.zeroCreate<Vec2> N

// let normal_mesh = Array2D.zeroCreate<Vec2> N N
// for i in 0..N-1 do
//     for j in 0..N-1 do
//         normal_mesh[i,j].X <- (A + ((D-A)/double N) * (double i)).X
//         normal_mesh[i,j] <- (B + ((C-B)/double N) * (double i)).X
//         normal_mesh[i,j] <- (A + ((B-A)/double N) * (double i)).X
//         normal_mesh[i,j] <- (D + ((C-D)/double N) * (double i)).X



for i in 0..N-1 do
    d[i] <- A + ((B-A)/double N) * (double i)
    u[i] <- D + ((C-D)/double N) * (double i)

    l[i] <- A + ((D-A)/double N) * (double i)
    r[i] <- B + ((C-B)/double N) * (double i)

    // printfn "%s, %s" (string r[i]) (string r[i])

let sides = GeoRandomizer.RandomSides(l, r, d, u)
let bounds = GeoRandomizer.SidesToBoundaryNodes(sides, N)
let mesh = GeoRandomizer.MeshFromRandomizer(sides, N, N)

GeoRandomizer.CheckUnitialized(mesh, N-1, N-1)

GeoRandomizer.WriteBoundaryNodes(sides, out_boundaries, N);

GeoRandomizer.WriteCorrectGrid(out_nodes, mesh, N-1, N-1)

// let bnodes = [|A; B; C; D; A|]
let bnodes = GeoRandomizer.SidesToBoundaryNodes(sides, N)

Gnuplot()
|> Gnuplot.datablockXY (bnodes |> Array.map (fun v -> v.X)) (bnodes |> Array.map (fun v -> v.Y)) "bounds"
|>> "set terminal png size 840,580"
|>> "set output 'correct_grid_compare_4.png'"
|>> "unset key"
// |>> $"plot '$bounds' with lines lc rgb 'black' lw 2"
|>> $"plot '$bounds' with lines lc rgb 'black' lw 2, \\"
|>> $"'{out_nodes}' with lines lc rgb 'black'"
|> Gnuplot.run
|> ignore

