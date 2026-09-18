nmax = 240

set xrange [0:3]
set grid
set encoding utf8
set terminal wxt font "Arial,12"
do for [k=0:nmax-1] {
    set multiplot layout 2,1

    # --- верх: Ez ---
    set ylabel 'E_z'
    set yrange [-1:1]
    plot 'vectorized_results/ez_vectorized.dat' index k using 1:2 with lines lc rgb 'red' title sprintf('Ez, кадр %d', k)

    # --- низ: Hy ---
    set ylabel 'H_y'
    set yrange [-3e-3:3e-3]
    plot 'vectorized_results/ez_vectorized.dat' index k using 1:2 with lines lc rgb 'blue' title 'Hy'

    unset multiplot
    pause 0.05
}