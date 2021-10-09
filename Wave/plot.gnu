set view 40, 30, 0.85, 2.
set samples 50,50
set isosamples 50,50
set cbrange [-1:1]
set xrange [-0.3:1.3]
set yrange [-0.3:1.3]
set zrange [-1:1]
set pm3d
t=0
dt = 0.01
do for [i=0:20] {
	t = t + dt
	filename(i) = sprintf(".u.dat.%03d", i)
	splot filename(i) w l title ''
	set title sprintf("T = %4.2f", t)
}