#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 4    last modified 2013-10-02 
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2013
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal wxt 0
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set timefmt z "%d/%m/%y,%H:%M"
set zdata 
set timefmt y "%d/%m/%y,%H:%M"
set ydata 
set timefmt x "%d/%m/%y,%H:%M"
set xdata 
set timefmt cb "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy u,v
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
unset grid
set raxis
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
set style line 3  linetype 3 linecolor rgb "grey50"  linewidth 1.000 pointtype 3 pointsize default pointinterval 0
set style line 4  linetype 4 linecolor rgb "grey70"  linewidth 1.000 pointtype 4 pointsize default pointinterval 0
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
set parametric
unset decimalsign
set view 65, 330, 1, 1
set samples 100, 100
set isosamples 100, 100
set surface
set contour both
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set cntrparam order 4
set cntrparam linear
set cntrparam levels incremental -1000,100,1200
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ 0.00000 : 5.00000 ] noreverse nowriteback
set vrange [ 0.00000 : 6.28319 ] noreverse nowriteback
set xlabel "" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ -5.00000 : 5.00000 ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ -5.00000 : 5.00000 ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ -1000.00 : 1000.00 ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_GB.UTF-8"
set pm3d explicit at s
set pm3d depthorder
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette defined ( 0 1 0 0, 0.499 1 0.6471 0, 0.5 0.302 0.302 0.302,\
     0.501 0 1 0, 1 0 0 1 )
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables
f(u,v) = 48.0 * (a**2 * cos(v)**2 - u**2) * (a**4 * cos(v)**4 - 18.0 * u**2 * a**2 * cos(v)**2 + u**4) / (a**2 * cos(v)**2 + u**2)**6
g(u,v) = (a**2 * cos(v)**2 + u**2)**6
clamp(u,v) = (f(u,v) > 1000.0) ? 1000.0 : (f(u,v) > -1000.0 ? f(u,v) : -1000.0)
GNUTERM = "wxt"
a = 1.0
GPFUN_f = "f(u,v) = 48.0 * (a**2 * cos(v)**2 - u**2) * (a**4 * cos(v)**4 - 14.0 * u**2 * a**2 * cos(v)**2 + u**4) / (a**2 * cos(v)**2 + u**2)**6"
GPFUN_g = "g(u,v) = (a**2 * cos(v)**2 + u**2)**6"
GPFUN_clamp = "clamp(u,v) = (f(u,v) > 1000.0) ? 1000.0 : (f(u,v) > -1000.0 ? f(u,v) : -1000.0)"
r = 10
th = 3.14159265358979
splot u*cos(v),u*sin(v),clamp(u,v) title "Kretschmann Scalar" ls 3 palette, (1.0 + sqrt(1.0 - a * a))*cos(v),(1.0 + sqrt(1.0 - a * a))*sin(v),0.0 title "Horizon"
#    EOF
