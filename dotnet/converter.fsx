#r "net10.0/GnuplotFS.dll"

open System
open Serializers
open System.IO

let table =
    File.ReadAllLines("../coordinates.dat")
    |> Array.map (fun line -> line.Split(',') |> List.ofArray)
    |> List.ofArray

Html.HtmlBuilder()
|> Html.table None [] [] table
|> Html.close("../cordinates_table.html")

