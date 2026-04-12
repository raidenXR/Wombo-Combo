#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open GridGeneration
open Plotting

let path = "../dat-files/coordinates_dense.dat"
let out_all      = "../dat-files/mesh_all.dat"
let out_vertices = "../dat-files/mesh_vertices.dat"
let out_indices  = "../dat-files/mesh_indices.dat"

// let mesh = Mesh.FromFile(path, 8);
// printfn "mesh.vertices.len: %d" mesh.vertices.Count

let vertices = Geo.ReadFromFile(path);
let centers = Geo.Centers(vertices);
let normals = Geo.Normals(vertices);
let bounds = BoundingBox.GetBounds(vertices)
let dx = bounds.v_max.X - bounds.v_min.X
let dy = bounds.v_max.Y - bounds.v_min.Y
printfn "dx: %g, dy: %g" dx dy
let mesh = Mesh.FromVertices(vertices, 40);

let cx = vertices |> Seq.map (fun v -> v.X) |> Array.ofSeq
let cy = vertices |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let px = centers |> Seq.map (fun v -> v.X) |> Array.ofSeq
let py = centers |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let nx = normals |> Seq.map (fun v -> v.X) |> Array.ofSeq
let ny = normals |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let vx = mesh.vertices |> Seq.map (fun v -> v.X) |> Array.ofSeq
let vy = mesh.vertices |> Seq.map (fun v -> v.Y) |> Array.ofSeq

do
    let fs = System.IO.File.CreateText(out_vertices)
    mesh.WriteVertices(fs)
    fs.Close()

Gnuplot(false, true, None)
|> Gnuplot.datablockXY vx vy "vertices"
|> Gnuplot.datablockXY px py "centers"
|> Gnuplot.datablockXY nx ny "normals"
|> Gnuplot.datablockXY cx cy "coordinates"
|>> $"plot '{out_vertices}' using 1:2 with lines lc rgb 'black',\\"
|>> "'$vertices' using 1:2 with points lc rgb 'black',\\"
|>> "'$centers' using 1:2 with points lc rgb 'red',\\"
|>> "'$normals' using 1:2 with points lc rgb 'blue',\\"
|>> "'$coordinates' using 1:2 with lp lw 2 lc rgb 'black'"
// |>> "plot '$vertices' using 1:2 with points lc rgb 'black',\\"
// |>> "'$coordinates' using 1:2 with lp lw 2 lc rgb 'black'"
|> Gnuplot.run
|> ignore

System.Console.ReadKey()


// do
//     let fs = System.IO.File.CreateText(out_all)
//     mesh.WriteAll(fs)
//     fs.Close()
// do
//     let fs = System.IO.File.CreateText(out_indices)
//     mesh.WriteIndices(fs)
//     fs.Close()
