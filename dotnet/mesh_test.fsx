#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open System
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

let N = 30
let mesh = Mesh.FromVertices(vertices, N);

let cx = vertices |> Seq.map (fun v -> v.X) |> Array.ofSeq
let cy = vertices |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let px = centers |> Seq.map (fun v -> v.X) |> Array.ofSeq
let py = centers |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let nx = normals |> Seq.map (fun v -> v.X) |> Array.ofSeq
let ny = normals |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let vx = mesh.vertices |> Seq.map (fun v -> v.X) |> Array.ofSeq
let vy = mesh.vertices |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let dv_x = (Array.max vx) - (Array.min vx)
let dv_y = (Array.max vy) - (Array.min vy)

let z = Array.zeroCreate<double> (vx.Length)

let antennas_x = Array.zeroCreate<double> 6 
let antennas_y = Array.zeroCreate<double> 6 
let signal = Array.zeroCreate<double> (mesh.vertices.Count)

do
    for i in 1..6 do
    let vs = mesh.vertices
    let a = Random.Shared.Next(vs.Count)
    let mutable b = a
    while a = b do b <- Random.Shared.Next(vs.Count)
    let p = Geo.Center(vs[a], vs[b])
    antennas_x[i-1] <- p.X
    antennas_y[i-1] <- p.Y
    
let gaussian_function A x x0 y y0 s_x s_y =
    let t1 = ((x-x0)**2.) / (2. * s_x**2.)
    let t2 = ((y-y0)**2.) / (2. * s_y**2.)
    A * exp(-(t1 + t2))

for j in 0..mesh.vertices.Count-1 do
    let x  = mesh.vertices[j].X 
    let y  = mesh.vertices[j].Y
    for i in 1..6 do
        let p = Vec2(antennas_x[i-1], antennas_y[i-1])
        let A  = 5.      // amplitude
        let x0 = p.X      // center x
        let y0 = p.Y      // center y
        // let s_x = 50
        // let s_y = 40
        let s_x = 100. * dv_x / double(mesh.vertices.Count)   // x spread of the blob
        let s_y = 100. * dv_y / double(mesh.vertices.Count)   // y spread of the blob
        signal[j] <- signal[j] + gaussian_function A x x0 y y0 s_x s_y
    
printfn "signal-> min: %g, max: %g, mean: %g" (Array.min signal) (Array.max signal) (((Array.max signal)-(Array.min signal))/(double signal.Length))    

do
    let fs = System.IO.File.CreateText(out_vertices)
    mesh.WriteVertices(fs)
    fs.Close()

Gnuplot(false, true, None)
|> Gnuplot.datablockXYZ vx vy signal "vertices"
// |> Gnuplot.datablockXY px py "centers"
// |> Gnuplot.datablockXY nx ny "normals"
|> Gnuplot.datablockXY cx cy "coordinates"
|> Gnuplot.datablockXY antennas_x antennas_y "antennas"
|>> "PIC = '../images/antenna_transparent.png'"
|>> "set for [i=1:6] pixmap i PIC center at word($antennas[i], 1), word($antennas[i], 2)"
|>> "unset key"
|>> "set isosample 250,250"
// |>> "set table 'antennas.dat'"
// |>> "splot '$antennas' using 1:2:(0) with points lw 3 lc rgb 'red'"
// |>> "unset table"
// |>> "set pm3d interpolate 0,0"
// |>> "set pm3d depth"
// |>> "set surface"
|>> "set view map"
// |>> "set dgrid3d"
// |>> "set contour base"
// |>> "set cntrparam cubicspline"
// |>> "unset surface"
|>> "set palette rgbformulae 33,13,10"
// |>> "set table '$vertices'"
// |>> "set dgrid3d 20,10 qnorm 4"
// |>> "set cntrparam level incremental -3,0.5,3"
// |>> "splot '$vertices' with pm3d, 'antennas.dat' with points lw 3 lc rgb 'red'"
|>> "splot $vertices using 1:2:3 with points palette lw 3 pt 3, \\"
// |>> $"'{out_vertices}' using 1:2:(0) with lines lc rgb 'black',\\"
// |>> "'$vertices' using 1:2:(0) with points lc rgb 'black',\\"
|>> "'$coordinates' using 1:2:(0) with lp lw 2 lc rgb 'black', \\"
// |>> "'$antennas' using 1:2:(0) with points lw 4 lc rgb 'red' notitle"
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
