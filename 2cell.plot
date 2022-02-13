
####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
set term pdfcairo enhanced size 7.5,1.*1.75 font "Helvetica,12
set out "Fig02_2cell.pdf"
set multi lay 2,2


unset label
unset object
set key bottom center
set xrange[0:75]
set yrange[-62:40]
#unset xrange
#unset yrange
set ytic 20
set xtic 10
set cbtic .2
unset border
unset tics
unset xlabel
unset ylabel

set label 1 "100s" at 575,-39 center
set arrow 1 from 525,-50 to 625,-50 front nohead lw 4


set xrange[500:800]
plot'data/2cell/1.sp'   using ($1-0):($2+50) w l lw 1 lc rgb "blue" title '',\
                             ''   using ($1-0):3 w l lw 1 lc rgb "red" title ''


unset arrow
unset label
set label 1 "10s" at 25,-39.0 center
set arrow 1 from 20,-50 to 30,-50 front nohead lw 4
set label 2 "40mV" at 11.,-35
set arrow 2 from 10.5,-55 to 10.5,-15 front nohead lw 4

set xrange[0:35]
plot'data/2cell/2.sp'   using ($1-107.5):($2+50) w l lw 1 lc rgb "blue" title '',\
                             ''   using ($1-107.5):3 w l lw 1 lc rgb "red" title ''

unset label 2
unset arrow 2

set label 1 "10s" at 25,-65.5 center
set arrow 1 from 20,-78 to 30,-78 front nohead lw 4
plot'data/2cell/3.sp'   using ($1-110):($2+50) w l lw 1 lc rgb "blue" title '',\
                             ''   using ($1-110):3 w l lw 1 lc rgb "red" title ''


plot'data/2cell/4.sp'   using ($1-110):($2+50) w l lw 1 lc rgb "blue" title '',\
                             ''   using ($1-110):3 w l lw 1 lc rgb "red" title ''
       
unset multi










