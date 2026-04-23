#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open System
open GridGeneration
open Plotting

let path = "../dat-files/coordinates_dense.dat"
let out_nodes = "../dat-files/mesh_vertices.dat"

let N = 70
let boundary_nodes = Geo.ReadFromFile(path)
let bounds   = Geo.GetBounds(boundary_nodes)
let mutable grid     = Geo.StencilFromBoundaryNodes(boundary_nodes, N)

// do
//     Writer.WriteBoundaryStencil("stencil_manual.txt", boundary_nodes, N, '*', ' ')
//     exit 0

grid.stencil_buffer <- Geo.ParseStencilBuffer("stencil_manual.txt", '*', N)

let _        = Geo.FillStencil(grid)
// let _        = Geo.FillStencil_2(grid,boundary_nodes)
let nodes    = Geo.MeshFromStencil(grid)
let antennas = Signal.InitialAntennasPositions(nodes, 6)
let signal_intensity = Array.zeroCreate<double> (N*N)
let I_opt = 8.

// do
//     Writer.WriteGridStencil("stencil_full.txt", grid, '*', ' ')
//     Writer.WriteVertices(out_nodes, grid)
//     Writer.WriteBoundaryStencil("stencil.txt", boundary_nodes, N)

// these are just for plotting
let bx = boundary_nodes |> Seq.map (fun v -> v.X) |> Array.ofSeq
let by = boundary_nodes |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let vx = nodes |> Seq.map (fun v -> v.X) |> Array.ofSeq
let vy = nodes |> Seq.map (fun v -> v.Y) |> Array.ofSeq
    
let centers = Geo.Centers(boundary_nodes);
let cx = centers |> Array.map (fun v -> v.X)
let cy = centers |> Array.map (fun v -> v.Y)

let normals = Geo.Normals(boundary_nodes);
let nx = normals |> Array.map (fun v -> v.X)
let ny = normals |> Array.map (fun v -> v.Y)

let cnx = [|for i in 0..cx.Length-1 -> nx[i] - cx[i]|]
let cny = [|for i in 0..cx.Length-1 -> ny[i] - cy[i]|]
    
// minimize this function
let obj_fn () =
    let dv_x = (Array.max vx) - (Array.min vx)
    let dv_y = (Array.max vy) - (Array.min vy)
    let mutable diff = 0.
    Array.Clear signal_intensity
    
    for i in 0..N-1 do
        for j in 0..N-1 do
            if grid.stencil_buffer[i*N + j] then
                let v = Geo.ToCartesianSystem(grid, i, j)
                let x = v.X
                let y = v.Y
                for k in 0..antennas.Length-1 do
                    let p = antennas[k]
                    let A  = 1.      // amplitude
                    let x0 = p.X      // center x
                    let y0 = p.Y      // center y
                    // let s_x = (Nn*Nn / 8.) * dv_x / double(nodes.Count)   // x spread of the blob
                    // let s_y = (Nn*Nn / 8.) * dv_y / double(nodes.Count)   // y spread of the blob
                    let s_x = 0.15 * dv_x
                    let s_y = 0.15 * dv_y
                    signal_intensity[i*N + j] <- signal_intensity[i*N + j] + Signal.Gaussian(A, x, x0, y, y0, s_x, s_y)
                // diff <- diff + abs(signal_intensity[i*N + j] - I_opt)
                diff <- diff + abs(signal_intensity[i*N + j] - 0.5)
    diff
    

// Initial estimation for Optimization alg - loop
Signal.InitialAntennasPositions(nodes, 6)

let gif = System.IO.File.CreateText("animation_gif.dat")
do
    for i in 0..5 do
        let x = antennas[i].X
        let y = antennas[i].Y
        gif.WriteLine($"{x}  {y}")
    gif.WriteLine("\n")

let n_iter_max = 60
let mutable n_iter = 0
while obj_fn() > 1e-3 && n_iter <= n_iter_max do
    let struct(ci,cj) = Console.GetCursorPosition()
    printfn "iteration: %d" n_iter
    Signal.MoveAntennas(grid, antennas, boundary_nodes, signal_intensity)
    Console.SetCursorPosition(0,cj)
        
    for i in 0..5 do
        let x = antennas[i].X
        let y = antennas[i].Y
        gif.WriteLine($"{x}  {y}")
        // let idx = Geo.ToStencilSystem(N, Vec2(x,y), bounds)
        // gif.WriteLine($"{idx.I}  {idx.J}")
        // printfn "%d, %d" idx.I idx.J
    // printf "\n"
    if (n_iter < n_iter_max) then gif.WriteLine("\n")
    n_iter <- n_iter + 1
gif.Close()

// Gnuplot()
// |> Gnuplot.datablockXY bx by "coordinates"
// |> Gnuplot.datablockXY vx vy "vertices"
// |> Gnuplot.datablockXYZW cx cy cnx cny "normals"
// |> Gnuplot.datablockXY cx cy "centers"
// |>> "set terminal png size 840,580"
// |>> "set output 'shape.png'"
// |>> "unset key"
// |>> "plot $coordinates using 1:2 with lp lc rgb 'black' lw 2, \\"
// |>> "$normals using 1:2:3:4 with vectors lc rgb 'red' lw 1 , \\"
// |>> "$centers using 1:2 with points lc rgb 'green' lw 2 "
// // |>> $"'{out_nodes}' using 1:2 with lines lc rgb 'black'"
// |> Gnuplot.run
// |> ignore

// Gnuplot()
// |> Gnuplot.datablockXY bx by "coordinates"
// |> Gnuplot.datablockXY vx vy "vertices"
// |> Gnuplot.datablockXY nx ny "normals"
// |>> "set terminal png size 840,580"
// |>> "set output 'grid.png'"
// |>> "unset key"
// |>> "plot $coordinates using 1:2 with lines lc rgb 'black' lw 2, \\"
// // |>> "$normals using 1:2 with lp lc rgb 'red' lw 2, \\"
// |>> $"'{out_nodes}' using 1:2 with lines lc rgb 'black'"
// |> Gnuplot.run
// |> ignore


Gnuplot()
|> Gnuplot.datablockXY bx by "coordinates"
|> Gnuplot.datablockXY vx vy "vertices"
|>> "set terminal gif animate delay 30"
|>> "set output 'antennas.gif'"
|>> "stats 'animation_gif.dat' nooutput"
|>> "do for [i=1:int(STATS_blocks)] {"
|>> "     plot 'animation_gif.dat' index (i-1) with circles, \\"
|>> $"        '{out_nodes}' using 1:2 with lines lc rgb 'black', \\"
|>> "         '$coordinates' using 1:2 with lines lw 2 lc rgb 'black'"
|>> "}"
|> Gnuplot.run
|> ignore


let si = [|
    for i in 0..N-1 do
        for j in 0..N-1 do
            if grid.stencil_buffer[i*N + j] then
                yield signal_intensity[i*N + j]
|]

Gnuplot()
|> Gnuplot.datablockXYZ vx vy si "vertices"
|> Gnuplot.datablockXY bx by "coordinates"
|> Gnuplot.datablockXY (antennas |> Array.map (fun v -> v.X)) (antennas |> Array.map (fun v -> v.Y)) "antennas"
|>> "set terminal png size 840,580"
|>> "set output 'antennas.png'"
|>> "PIC = '../images/antenna_transparent.png'"
|>> "set for [i=1:6] pixmap i PIC center at word($antennas[i], 1), word($antennas[i], 2)"
|>> "unset key"
|>> "set isosample 250,250"
|>> "set view map"
|>> "set palette rgbformulae 33,13,10"
|>> "splot $vertices using 1:2:3 with lp palette lw 2 pt 5, \\"
|>> "'$coordinates' using 1:2:(0) with lines lw 2 lc rgb 'black'"
|> Gnuplot.run
|> ignore


