set terminal pdfcairo font "Times-Roman,7"
set key right bottom
set termoption dashed
set pointsize 0.5
set style line 1 lt 1 linecolor rgb "#C00000" lw 3 pt 5
set style line 2 lt 2 linecolor rgb "#00A000" lw 3 pt 5
set style line 3 lt 3 linecolor rgb "#4060D0" lw 3 pt 5
set style line 4 lt 4 linecolor rgb "#825900" lw 3 pt 5
set style line 5 lt 5 linecolor rgb "#9932CC" lw 3 pt 5
set style line 6 lt 6 linecolor rgb "#DAA520" lw 3 pt 5
set output "graphs/conf_acc.pdf"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/conf_acc_3_100" using 1:2 w linespoints title "3 workers 100 tasks" ls 1, "plotdata/conf_acc_3_300" using 1:2 w linespoints title "3 workers 300 tasks" ls 2, "plotdata/conf_acc_7_100" using 1:2 w linespoints title "7 workers 100 tasks" ls 3, "plotdata/conf_acc_7_300" using 1:2 w linespoints title "7 workers 300 tasks" ls 4
set output "graphs/dens_unc.pdf"
set key right top
set xlabel "Density"
set ylabel "Size of Interval"
plot "plotdata/dens_unc_3_300" using 1:2 w linespoints title "3 workers, 300 tasks" ls 1, "plotdata/dens_unc_7_100" using 1:2 w linespoints title "7 workers, 100 tasks" ls 2, "plotdata/dens_unc_7_300" using 1:2 w linespoints title "7 workers, 300 tasks" ls 3
set output "graphs/conf_acc_real.pdf"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/conf_acc_imcomp" using 1:2 w linespoints title "Image Comparison" ls 1, "plotdata/conf_acc_rte" using 1:2 w linespoints title "RTE" ls 2, "plotdata/conf_acc_temporal" using 1:2 w linespoints title "Temporal" ls 3
set output "graphs/conf_acc_real_improved.pdf"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/conf_acc_imcomp_improved" using 1:2 w linespoints title "Image Comparison" ls 1, "plotdata/conf_acc_rte_improved" using 1:2 w linespoints title "RTE" ls 2, "plotdata/conf_acc_temporal_improved" using 1:2 w linespoints title "Temporal" ls 3
set output "graphs/conf_unc.pdf"
set key left top
set xlabel "Confidence Level"
set ylabel "Size of Interval"
plot "plotdata/conf_unc_3_100" using 1:2 w linespoints title "new technique, 3 workers, 100 tasks" ls 1,  "plotdata/conf_unc_3_100" using 1:3 w linespoints title "old technique, 3 workers, 100 tasks" ls 2,  "plotdata/conf_unc_7_100" using 1:2 w linespoints title "new technique, 7 workers, 100 tasks" ls 3,  "plotdata/conf_unc_7_100" using 1:3 w linespoints title "old technique, 7 workers, 100 tasks" ls 4
set output "graphs/conf_unc_optimization.pdf"
set key left top
set xlabel "Confidence Level"
set ylabel "Size of Interval"
plot "plotdata/conf_unc_noopt_7_100" using 1:3 w linespoints title "No Optimization" ls 1,  "plotdata/conf_unc_opt_7_100" using 1:3 w linespoints title "With Optimization" ls 2
set output "graphs/conf_unc_realdata_optimization.pdf"
set key left top
set xlabel "Confidence Level"
set ylabel "Size of Interval"
plot "plotdata/conf_acc_imcomp_noopt" using 1:3 w linespoints title "Image Comparison without Optimization" ls 1, "plotdata/conf_acc_rte_noopt" using 1:3 w linespoints title "RTE Comparison without Optimization" ls 2, "plotdata/conf_acc_temporal_noopt" using 1:3 w linespoints title "Temporal without Optimization" ls 3, "plotdata/conf_acc_imcomp" using 1:3 w linespoints title "Image Comparison" ls 4, "plotdata/conf_acc_rte" using 1:3 w linespoints title "RTE" ls 5, "plotdata/conf_acc_temporal" using 1:3 w linespoints title "Temporal"  ls 6 

#### kary graphs begin below

set output "graphs/k_conf_acc.pdf"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/k_conf_acc_2_100" using 1:2 w linespoints title "arity 2, 100 tasks", "plotdata/k_conf_acc_3_100" using 1:2 w linespoints title "arity 3, 100 tasks" ls 1, "plotdata/k_conf_acc_4_100" using 1:2 w linespoints title "arity 4, 100 tasks" ls 2, "plotdata/k_conf_acc_2_1000" using 1:2 w linespoints title "arity 2, 1000 tasks" ls 3, "plotdata/k_conf_acc_3_1000" using 1:2 w linespoints title "arity 3, 1000 tasks" ls 4, "plotdata/k_conf_acc_4_1000" using 1:2 w linespoints title "arity 4, 1000 tasks" ls 5 
set output "graphs/k_dens_unc.pdf"
set key right top
set xlabel "Density"
set ylabel "Average Size of Interval" 
plot "plotdata/k_dens_unc_2_0.8_500" using 1:3 w linespoints title "Arity 2" ls 1, "plotdata/k_dens_unc_3_0.8_500" using 1:3 w linespoints title "Arity 3" ls 2, "plotdata/k_dens_unc_4_0.8_500" using 1:3 w linespoints title "Arity 4" ls 3 
set output "graphs/k_conf_acc_real.pdf"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/k_conf_mooc_3" using 1:2 w linespoints title "MOOC arity 3" ls 1, "plotdata/k_conf_wsd_2" using 1:2 w linespoints title "WSD arity 2" ls 2, "plotdata/k_conf_wordsim_2" using 1:2 w linespoints title "Wordsim arity 2" ls 3
