# plot output of test_fresnel.f90, showing the fresnel integrals and
# cornu spiral

set terminal pdfcairo

set xlabel 'x'

set output 'fresnel_integrals_normalized.pdf'
set yrange [0 to 1]
set xrange [0 to 5]
set title 'Fresnel Integrals (Normalized)'
plot '../output/fresnel.dat' u 1:2 with lines linewidth 1.5 title 'C(x)', \
	'../output/fresnel.dat' u 1:3 with lines linewidth 1.5 title 'S(x)', \
	0.5 with lines linestyle 0 linewidth 1.5 title ''
set output

set output 'cornu_spiral.pdf'
set grid
set yrange [-1 to 1]
set xrange [-1 to 1]
set title 'Cornu Spiral'
plot '../output/fresnel.dat' u 2:3 with lines linestyle -1 linewidth 1.5 title ''
unset grid
set output
