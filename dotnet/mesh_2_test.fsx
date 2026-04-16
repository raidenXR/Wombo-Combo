#r "bin/Debug/net10.0/msip.dll"
#r "net10.0/GnuplotFS.dll"

open System
open GridGeneration
open Plotting

let path = "../dat-files/coordinates_dense.dat"
let out_nodes = "../dat-files/mesh_vertices.dat"

let N = 70
let boundary_nodes = Geo.ReadFromFile(path)
let bounds  = Geo.GetBounds(boundary_nodes)
let stencil = Geo.StencilFromBoundaryNodes(boundary_nodes, N)
Geo.FillStencil(stencil, N, bounds,Marching.Horizontal)
let nodes   = Geo.MeshFromStencil(stencil, N, bounds)

do
    Writer.WriteVertices(out_nodes, stencil, N, bounds)
    Writer.WriteBoundaryStencil("stencil.txt", boundary_nodes, N)


// these are just for plotting
let bx = boundary_nodes |> Seq.map (fun v -> v.X) |> Array.ofSeq
let by = boundary_nodes |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let vx = nodes |> Seq.map (fun v -> v.X) |> Array.ofSeq
let vy = nodes |> Seq.map (fun v -> v.Y) |> Array.ofSeq
let dv_x = (Array.max vx) - (Array.min vx)
let dv_y = (Array.max vy) - (Array.min vy)

let z = Array.zeroCreate<double> (vx.Length)

let antennas_x = Array.zeroCreate<double> 6 
let antennas_y = Array.zeroCreate<double> 6 
let signal = Array.zeroCreate<double> (nodes.Count)
let I_opt = 8.

let move_antennas () =
    let dx = (bounds.v_max.X - bounds.v_min.X) / double(N)
    let dy = (bounds.v_max.Y - bounds.v_min.Y) / double(N)

    for i in 0..5 do
        let p = Vec2(antennas_x[i], antennas_y[i])
        let p_idx = Geo.ToStencilSystem(N, p, bounds)
        let grad  = Geo.GetGradientAtNode(stencil, N, p_idx.I, p_idx.J, bounds, signal)

        let dir =
            if signal[Geo.GetIndex(stencil, N, p_idx.I, p_idx.J)] - I_opt > 0 then grad.SteepestDescend()
            elif signal[Geo.GetIndex(stencil, N, p_idx.I, p_idx.J)] - I_opt < 0 then grad.SteepestAscend()
            else GradDirection.Self
            
        match dir with
        | GradDirection.UL ->
            antennas_x[i] <- antennas_x[i] - dx
            antennas_y[i] <- antennas_y[i] + dy
        | GradDirection.UC ->
            antennas_y[i] <- antennas_y[i] + dy
        | GradDirection.UR ->
            antennas_x[i] <- antennas_x[i] + dx
            antennas_y[i] <- antennas_y[i] + dy
        | GradDirection.CL ->
            antennas_x[i] <- antennas_x[i] - dx
        | GradDirection.CR ->
            antennas_x[i] <- antennas_x[i] + dx
        | GradDirection.LL ->
            antennas_x[i] <- antennas_x[i] - dx
            antennas_y[i] <- antennas_y[i] - dy
        | GradDirection.LC ->
            antennas_y[i] <- antennas_y[i] - dy
        | GradDirection.LR ->
            antennas_x[i] <- antennas_x[i] + dx
            antennas_y[i] <- antennas_y[i] - dy
        | _ -> ()          
            

    
// minimize this function
let obj_fn () =
    let gaussian_function A x x0 y y0 s_x s_y =
        let t1 = ((x-x0)**2.) / (2. * s_x**2.)
        let t2 = ((y-y0)**2.) / (2. * s_y**2.)
        A * exp(-(t1 + t2))

    let mutable diff = 0.
    
    for j in 0..nodes.Count-1 do
        let x  = nodes[j].X 
        let y  = nodes[j].Y
        let Nn = double N
        // clear signals buffer
        signal[j] <- 0.
        for i in 0..5 do
            let p = Vec2(antennas_x[i], antennas_y[i])
            let A  = 5.      // amplitude
            let x0 = p.X      // center x
            let y0 = p.Y      // center y
            let s_x = (Nn*Nn / 8.) * dv_x / double(nodes.Count)   // x spread of the blob
            let s_y = (Nn*Nn / 8.) * dv_y / double(nodes.Count)   // y spread of the blob
            signal[j] <- signal[j] + gaussian_function A x x0 y y0 s_x s_y
        diff <- diff + abs(signal[j] - I_opt)
    // printfn "signal-> min: %g, max: %g, mean: %g" (Array.min signal) (Array.max signal) (((Array.max signal)-(Array.min signal))/(double signal.Length))    
    diff
    

// Initial estimation for Optimization alg - loop
do
    for i in 1..6 do
    let a = Random.Shared.Next(nodes.Count)
    let mutable b = a
    while a = b do b <- Random.Shared.Next(nodes.Count)
    let p = Geo.Center(nodes[a], nodes[b])
    antennas_x[i-1] <- p.X
    antennas_y[i-1] <- p.Y

let gif = System.IO.File.CreateText("animation_gif.dat")
let n_iter_max = 30
let mutable n_iter = 0
while obj_fn() > 1e-3 && n_iter <= n_iter_max do
    move_antennas()
    // create an animation of the loop
    
    for i in 0..5 do
        let x = antennas_x[i]
        let y = antennas_y[i]
        gif.WriteLine($"{x}  {y}")
        // let idx = Geo.ToStencilSystem(N, Vec2(x,y), bounds)
        // gif.WriteLine($"{idx.I}  {idx.J}")
        // printfn "%d, %d" idx.I idx.J
    // printf "\n"
    if (n_iter < n_iter_max) then gif.WriteLine("\n")
    n_iter <- n_iter + 1
gif.Close()

Gnuplot(false, true, None)
|> Gnuplot.datablockXY bx by "coordinates"
|> Gnuplot.datablockXY vx vy "vertices"
|>> "set terminal gif animate delay 50"
|>> "set output 'antennas.gif'"
|>> "stats 'animation_gif.dat' nooutput"
// |>> $"set xrange [0:{N}]"
// |>> $"set yrange [0:{N}]"
|>> "do for [i=1:int(STATS_blocks)] {"
|>> "     plot 'animation_gif.dat' index (i-1) with circles, \\"
|>> $"        '{out_nodes}' using 1:2 with lines lc rgb 'black', \\"
|>> "         '$coordinates' using 1:2 with lines lw 2 lc rgb 'black'"
|>> "}"
|> Gnuplot.run
|> ignore


Gnuplot(false, true, None)
|> Gnuplot.datablockXYZ vx vy signal "vertices"
|> Gnuplot.datablockXY bx by "coordinates"
|> Gnuplot.datablockXY antennas_x antennas_y "antennas"
|>> "set terminal png size 840,580"
|>> "set output 'antennas.png'"
|>> "PIC = '../images/antenna_transparent.png'"
|>> "set for [i=1:6] pixmap i PIC center at word($antennas[i], 1), word($antennas[i], 2)"
|>> "unset key"
|>> "set isosample 250,250"
|>> "set view map"
|>> "set palette rgbformulae 33,13,10"
|>> "splot $vertices using 1:2:3 with points palette lw 3 pt 5, \\"
// |>> $"'{out_nodes}' using 1:2:(0) with lines lc rgb 'black',\\"
// |>> "'$vertices' using 1:2:(0) with points lc rgb 'black',\\"
|>> "'$coordinates' using 1:2:(0) with lines lw 2 lc rgb 'black', \\"
// |>> "'$coordinates' using 1:2:(0) with lp lw 2 lc rgb 'black', \\"
// |>> "'$antennas' using 1:2:(0) with points lw 4 lc rgb 'red' notitle"
|> Gnuplot.run
|> ignore

// System.Console.ReadKey()

