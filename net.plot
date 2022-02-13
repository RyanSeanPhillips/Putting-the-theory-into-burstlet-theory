#####################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
set term pdfcairo enhanced size 7.5,4.*1.75 font "Helvetica,12
set out "Fig04_net.pdf"
set multi lay 4,1

set style line 101 lc rgb 'black' lt 1 lw 2
set border 3 front ls 101 lw 4 #(1 = bottom, 2= left, 4=top, 8 = right *ex x1y1 boarder = 3, x1y1y2 boarder = 11 )
set tics nomirror out scale 1

set ylabel"spks/s/N"
set xlabel"Time (s)"
set ytic 20
set xtic 20
set yrange[0:80]
set xrange[0:300]

plot "data/net/1.sp"   using ($1):2 w l lw .25 lc rgb "blue" title""

plot "data/net/2.sp"   using ($1):2 w l lw .25 lc rgb "blue" title""

plot "data/net/3.sp"   using ($1):2 w l lw .25 lc rgb "blue" title""

plot "data/net/4.sp"   using ($1):2 w l lw .25 lc rgb "blue" title""


unset multi










