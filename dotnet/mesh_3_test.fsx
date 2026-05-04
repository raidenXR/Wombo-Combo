#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open System
open GridGeneration
open Plotting

let path = "../dat-files/coordinates_dense.dat"
let out_nodes = "../dat-files/correct_grid.dat"
let out_nodes_2 = "../dat-files/correct_grid_compare.dat"

let N = 50
let boundary_nodes = Geo.ReadFromFile(path)
printfn "boundary_nodes.len: %d" (boundary_nodes.Length)
let bounds   = Geo.GetBounds(boundary_nodes)
let mutable grid     = Geo.StencilFromBoundaryNodes(boundary_nodes, N)

// do
//     Writer.WriteBoundaryStencil("correct_stencil_manual.txt", boundary_nodes, N, '*', ' ')
//     exit 0
    
grid.stencil_buffer <- Geo.ParseStencilBuffer("correct_stencil_manual.txt", '*', N)
let _        = Geo.FillStencil(grid)
let nodes    = Geo.MeshFromStencil(grid)

let mesh = Geo.CreateCorrectWhateverGrid(grid, N)

Writer.WriteCorrectGrid(out_nodes, mesh, N)
Writer.WriteVertices(out_nodes_2, grid)

let bx = boundary_nodes |> Seq.map (fun v -> v.X) |> Array.ofSeq
let by = boundary_nodes |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let vx = nodes |> Seq.map (fun v -> v.X) |> Array.ofSeq
let vy = nodes |> Seq.map (fun v -> v.Y) |> Array.ofSeq

Gnuplot()
|> Gnuplot.datablockXY bx by "coordinates"
|>> "set terminal png size 840,580"
|>> "set output 'correct_grid_compare.png'"
|>> "unset key"
// |>> "plot $coordinates using 1:2 with lines lc rgb 'black' lw 2, \\"
// |>> "$normals using 1:2 with lp lc rgb 'red' lw 2, \\"
|>> $"plot '{out_nodes_2}' using 1:2 with lines lc rgb 'black'"
|> Gnuplot.run
|> ignore

Gnuplot()
|> Gnuplot.datablockXY bx by "coordinates"
|>> "set terminal png size 840,580"
|>> "set output 'correct_grid.png'"
|>> "unset key"
|>> "plot $coordinates using 1:2 with lines lc rgb 'black' lw 2, \\"
// |>> "$normals using 1:2 with lp lc rgb 'red' lw 2, \\"
|>> $"'{out_nodes}' using 1:2 with lines lc rgb 'black'"
|> Gnuplot.run
|> ignore
