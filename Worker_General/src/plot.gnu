set terminal postscript enhanced
set output "graphs/conf_acc.ps"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/conf_acc_3_100" using 1:2 w linespoints title "3 workers 100 tasks", "plotdata/conf_acc_3_300" using 1:2 w linespoints title "3 workers 300 tasks", "plotdata/conf_acc_7_100" using 1:2 w linespoints title "7 workers 100 tasks", "plotdata/conf_acc_7_300" using 1:2 w linespoints title "7 workers 300 tasks"
set output "graphs/dens_unc.ps"
set key right top
set xlabel "Density"
set ylabel "Size of Interval"
plot "plotdata/dens_unc_3_300" using 1:2 w linespoints title "3 workers, 300 tasks", "plotdata/dens_unc_7_100" using 1:2 w linespoints title "7 workers, 100 tasks", "plotdata/dens_unc_7_300" using 1:2 w linespoints title "7 workers, 300 tasks"
set output "graphs/conf_acc_real.ps"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/conf_acc_imcomp" using 1:2 w linespoints title "Image Comparison", "plotdata/conf_acc_rte" using 1:2 w linespoints title "RTE", "plotdata/conf_acc_temporal" using 1:2 w linespoints title "Temporal"
set output "graphs/conf_acc_real_improved.ps"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/conf_acc_imcomp_improved" using 1:2 w linespoints title "Image Comparison", "plotdata/conf_acc_rte_improved" using 1:2 w linespoints title "RTE", "plotdata/conf_acc_temporal_improved" using 1:2 w linespoints title "Temporal"
set output "graphs/conf_unc.ps"
set key left top
set xlabel "Confidence Level"
set ylabel "Size of Interval"
plot "plotdata/conf_unc_3_100" using 1:2 w linespoints title "new technique, 3 workers, 100 tasks",  "plotdata/conf_unc_3_100" using 1:3 w linespoints title "old technique, 3 workers, 100 tasks",  "plotdata/conf_unc_7_100" using 1:2 w linespoints title "new technique, 7 workers, 100 tasks",  "plotdata/conf_unc_7_100" using 1:3 w linespoints title "old technique, 7 workers, 100 tasks"
set output "graphs/conf_unc_optimization.ps"
set key left top
set xlabel "Confidence Level"
set ylabel "Size of Interval"
plot "plotdata/conf_unc_noopt_7_100" using 1:3 w linespoints title "No Optimization",  "plotdata/conf_unc_opt_7_100" using 1:3 w linespoints title "With Optimization"
set output "graphs/conf_unc_realdata_optimization.ps"
set key left top
set xlabel "Confidence Level"
set ylabel "Size of Interval"
plot "plotdata/conf_acc_imcomp_noopt" using 1:3 w linespoints title "Image Comparison without Optimization", "plotdata/conf_acc_rte_noopt" using 1:3 w linespoints title "RTE Comparison without Optimization", "plotdata/conf_acc_temporal_noopt" using 1:3 w linespoints title "Temporal without Optimization", "plotdata/conf_acc_imcomp" using 1:3 w linespoints title "Image Comparison", "plotdata/conf_acc_rte" using 1:3 w linespoints title "RTE", "plotdata/conf_acc_temporal" using 1:3 w linespoints title "Temporal"  

#### kary graphs begin below

set output "graphs/k_conf_acc.ps"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/k_conf_acc_2_100" using 1:2 w linespoints title "arity 2, 100 tasks", "plotdata/k_conf_acc_3_100" using 1:2 w linespoints title "arity 3, 100 tasks", "plotdata/k_conf_acc_4_100" using 1:2 w linespoints title "arity 4, 100 tasks", "plotdata/k_conf_acc_2_1000" using 1:2 w linespoints title "arity 2, 1000 tasks", "plotdata/k_conf_acc_3_1000" using 1:2 w linespoints title "arity 3, 1000 tasks", "plotdata/k_conf_acc_4_1000" using 1:2 w linespoints title "arity 4, 1000 tasks" 
set output "graphs/k_dens_unc.ps"
set key right top
set xlabel "Density"
set ylabel "Average Size of Interval" 
plot "plotdata/k_dens_unc_2_0.8_500" using 1:3 w linespoints title "Arity 2", "plotdata/k_dens_unc_3_0.8_500" using 1:3 w linespoints title "Arity 3", "plotdata/k_dens_unc_4_0.8_500" using 1:3 w linespoints title "Arity 4" 
set output "graphs/k_conf_acc_real.ps"
set key right bottom
set xlabel "Confidence Level"
set ylabel "Accuracy"
plot x title "Ideal interval-accuracy", "plotdata/k_conf_mooc_3" using 1:2 w linespoints title "MOOC arity 3", "plotdata/k_conf_wsd_2" using 1:2 w linespoints title "WSD arity 2", "plotdata/k_conf_wordsim_2" using 1:2 w linespoints title "Wordsim arity 2"
