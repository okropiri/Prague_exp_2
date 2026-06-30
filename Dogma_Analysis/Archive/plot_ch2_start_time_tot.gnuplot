if (!exists("infile")) infile = "ch2_start_time_tot_hist2d.dat"
if (!exists("outfile")) outfile = "ch2_start_time_tot_hist2d.png"
if (!exists("plot_title")) plot_title = "NCAL ch2 TOT vs raw leading time, no modulo"
if (!exists("xmin")) xmin = -10000
if (!exists("xmax")) xmax = 10000
if (!exists("ymin")) ymin = 0
if (!exists("ymax")) ymax = 128

set terminal pngcairo size 1800,1000 enhanced font "DejaVu Sans,14"
set output outfile
set title plot_title
set xlabel "ch2 leading time [ns]"
set ylabel "ch2 TOT [ns]"
set cblabel "counts / bin"
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set view map
set pm3d map
set palette defined (0 "black", 1 "navy", 2 "blue", 3 "cyan", 4 "yellow", 5 "red", 6 "white")
set logscale cb
set cbrange [1:*]
unset key
set grid xtics ytics lc rgb "#444444"

splot infile using 1:2:($3 > 0 ? $3 : 1/0) with pm3d