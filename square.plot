################################################################################################################################
##################################################################################################################################
################################################################################################################################
##################################################################################################################################
set term pdfcairo enhanced size 0.5*8.5,2.*2.5 font "Helvetica,12
set out "Fig01_sq_wave.pdf"
set multi lay 4,2


set style line 101 lc rgb 'black' lt 1 lw 2
set border 3 front ls 101 lw 4 #(1 = bottom, 2= left, 4=top, 8 = right *ex x1y1 boarder = 3, x1y1y2 boarder = 11 )
set tics nomirror out scale 1




set xrange[0:200]
set xtic 10
set y2range[0:0.085]
set yrange[-0.00004:0.00019]
set y2tic .1
set ytic .00005

#unset yrange
unset border
unset xtic
unset ytic
unset y2tic

set size 2*0.48,1./4
set label 1 '50s' at 25,0.000145 center front font "Arial-Bold,12
set arrow 1 from 0,0.000125 to 50,0.000125 nohead lw 6
plot'data/square/1.sp'   using ($1-400):(1*$3) w l lw 4 lc rgb "green"  axis x1y2 title 'Ca2+ Pulse',\
                    ''   using ($1-400):4 w l lw 4 lc rgb "blue" title 'Cai',\
                    ''   using ($1-400):($5/10) w l lw 3 lc rgb "red-orange"  axis x1y1 title 'CaER'
set multiplot next 
set xrange[0:40]
set label 1 '10s' at 15,0.000145 center front font "Arial-Bold,12
set arrow 1 from 10,0.000125 to 20,0.000125 nohead lw 6
set label 2 '0.5{/Symbol m}M' at 40,0.00009+0.00005 center tc rgb "blue" front font "Arial-Bold,12
set label 3 '5.0{/Symbol m}M' at 40,0.00006+0.00005 center tc rgb "red-orange" front font "Arial-Bold,12
set arrow 2 from 38,0.00005+0.00005 to 38,0.0001+0.00005 nohead lw 6     
set size 2*0.48,1./4                
plot'data/square/2.sp'   using ($1-240):(1*$3) w l lw 4 lc rgb "green"  axis x1y2 title '',\
                    ''   using ($1-240):4 w l lw 4 lc rgb "blue" title '',\
                    ''   using ($1-240):($5/10) w l lw 3 lc rgb "red-orange"  axis x1y1 title ''


unset arrow
unset label
set label 1 '10s' at 15,0.000185 center front font "Arial-Bold,12
set arrow 1 from 10,0.000165 to 20,0.000165 nohead lw 6
unset y2tic
set multiplot next           
set size 2*0.48,1./4           
plot'data/square/3.sp'   using ($1-250):(1*$3) w l lw 4 lc rgb "green"  axis x1y2 title '',\
                    ''   using ($1-250):4 w l lw 4 lc rgb "blue" title '',\
                    ''   using ($1-250):($5/10) w l lw 3 lc rgb "red-orange"  axis x1y1 title ''
set multiplot next        
set size 2*0.48,1./4              
plot'data/square/4.sp'   using ($1-250-1):(1*$3) w l lw 4 lc rgb "green"  axis x1y2 title '',\
                    ''   using ($1-250-1):4 w l lw 4 lc rgb "blue" title '',\
                    ''   using ($1-250-1):($5/10) w l lw 3 lc rgb "red-orange"  axis x1y1 title ''

unset arrow
unset label




unset multi








