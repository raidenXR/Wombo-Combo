# set term png size 640,480
# set output 'mesh.png'

plot 'mesh_vertices.dat' using 1:2 with points lc rgb 'black', \
    'coordinates_dense_tabbed.dat' using 2:3 with lines lw 2 lc rgb 'black'
    
# plot 'coordinates_dense_tabbed.dat' using 2:3 with lines lw 2 lc rgb 'black'

