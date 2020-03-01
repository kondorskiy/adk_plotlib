
/*====================================================================

  C++ subroutines to plot 2D and 3D real functions using GNUplot
  and the test & example program.

  Alexey D. Kondorskiy,
  P.N.Lebedev Physical Institute of the Russian Academy of Science.
  E-mail: kondorskiy@lebedev.ru, kondorskiy@gmail.com.

====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>


/*--------------------------------------------------------------------
  Minimal possible value to plot in logarythmic scale.
--------------------------------------------------------------------*/
const double LOG_THRESH = 1.0e-30;


/*--------------------------------------------------------------------
  Width of the figures in pixels.
--------------------------------------------------------------------*/
const int X_SIZE = 1280;


/*--------------------------------------------------------------------
  Print tabulated 2D function "func" to the file.
--------------------------------------------------------------------*/
void print2Ddata(
  double (*func)(double),         // Function to plot.
  const double &x_min,            // Minimal value of argument.
  const double &x_max,            // Maximal value of argument.
  const int &x_num,               // Number of argument points.
  const std::string &name,        // Data identification name.
  double &y_min,                  // Result minimal value of function.
  double &y_max)                  // Result maximal value of function.
{
  std::string file_name = name + ".dat";
  std::ofstream fout(file_name.c_str(), std::ios::out);
  for(int i = 0; i < x_num; ++i) {
    double x = x_min + i*(x_max - x_min)/(x_num - 1);

    double y = (*func)(x);
    fout << x << " " << y << "\n";

    if(i == 0){
      y_min = y;
      y_max = y;
    } else {
      if(y_min > y) y_min = y;
      if(y_max < y) y_max = y;
    }
  }
  fout.close();
}


/*--------------------------------------------------------------------
  Print tabulated 3D polar surface for function "func(phi, a)"
  with a scan parameter, a, varying along the cylinder axis.
  ARGUMENT phi = 0..360. SHOULD BE IN DEGREES! a - is the parameter.
--------------------------------------------------------------------*/
void print2DdataScanXYZ(
  double (*func)(double, double), // Function to plot.
  const int &x_num,             // Number of argument points.
  const double &a_min,          // Minimal value of the parameter.
  const double &a_max,          // Maximal value of the parameter.
  const int &a_num,             // Number of points for the parameter.
  const std::string &name,      // Data identification name.
  double &r_min,                // Result minimal value of function.
  double &r_max)                // Result maximal value of function.
{
  std::string file_name = name + ".dat";
  std::ofstream fout(file_name.c_str(), std::ios::out);

  double d2r = M_PI/180.0;

  for(int iz = 0; iz < a_num; ++iz) {     // Parameter.
    double z = a_min + iz*(a_max - a_min)/(a_num - 1);

    for(int ia = 0; ia < x_num; ++ia) {    // Azimuth.
      double phi = ia*360.0/(x_num - 1);

      double r = (*func)(phi, z);

      double x = r*cos(phi*d2r);
      if( (ia == 0) && (iz == 0) ) {
        r_min = x;
        r_max = x;
      } else {
        if(x < r_min) r_min = x;
        if(x > r_max) r_max = x;
      }

      double y = r*sin(phi*d2r);
      if(y < r_min) r_min = y;
      if(y > r_max) r_max = y;

      fout << x << " " << y << " " << z << "\n";
    }                                       // Azimuth.
    fout << "\n";
  }                                       // Parameter.
  fout.close();
}


/*--------------------------------------------------------------------
  Print tabulated 3D function "func" to the file.
--------------------------------------------------------------------*/
void print3Ddata(
  double (*func)(double, double), // Function to plot.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const int &x_num,               // Number of X-argument points.
  const double &y_min,            // Minimal value of Y-argument.
  const double &y_max,            // Maximal value of Y-argument.
  const int &y_num,               // Number of Y-argument points.
  const std::string &name,        // Data identification name.
  const bool &is_log_scale,       // Flag to plot data in LOG scale.
  double &z_min,                  // Result minimal value of function.
  double &z_max)                  // Result maximal value of function.
{
  std::string file_name = name + ".dat";
  std::ofstream fout(file_name.c_str(), std::ios::out);
  for(int i = 0; i < x_num; ++i) {
    double x = x_min + i*(x_max - x_min)/(x_num - 1);

    for(int j = 0; j < y_num; ++j) {
      double y = y_min + j*(y_max - y_min)/(y_num - 1);

      fout << x << " " << y << " ";
      double z = (*func)(x, y);
      if( is_log_scale && (z < LOG_THRESH) )
        fout << LOG_THRESH << "\n";
      else
        fout << z << "\n";

      if(i == 0){
        z_min = z;
        z_max = z;
      } else {
        if(z_min > z) z_min = z;
        if(z_max < z) z_max = z;
      }

    }
    if(i < x_num) fout << "\n";
  }
  fout.close();
}


/*--------------------------------------------------------------------
  Print tabulated 3D polar surface for function "func(theta, phi)".
  ANGLES SHOULD BE IN DEGREES! theta = 0..180; phi = 0..360.
--------------------------------------------------------------------*/
void print3DpolarXYZ(
  double (*func)(double, double), // Function to plot.
  const int &aa_num,              // Number of points for azimuth.
  const int &pa_num,              // Number of points for polar angle.
  const std::string &name,        // Data identification name.
  double &r_min,                  // Result minimal value of function.
  double &r_max)                  // Result maximal value of function.
{
  std::string file_name = name + ".dat";
  std::ofstream fout(file_name.c_str(), std::ios::out);

  double d2r = M_PI/180.0;

  for(int ia = 0; ia < aa_num; ++ia) {      // Azimuth.
    double phi = ia*360.0/(aa_num - 1);
    for(int ip = 0; ip < pa_num; ++ip) {      // Polar angle.
      double theta = ip*180.0/(pa_num - 1);

      double r = (*func)(theta, phi);

      double x = r*sin(theta*d2r)*cos(phi*d2r);
      if( (ia == 0) && (ip == 0) ) {
        r_min = x;
        r_max = x;
      } else {
        if(x < r_min) r_min = x;
        if(x > r_max) r_max = x;
      }

      double y = r*sin(theta*d2r)*sin(phi*d2r);
      if(y < r_min) r_min = y;
      if(y > r_max) r_max = y;

      double z = r*cos(theta*d2r);
      if(z < r_min) r_min = z;
      if(z > r_max) r_max = z;

      fout << x << " " << y << " " << z << "\n";
    }                 // Polar angle.
    fout << "\n";
  }                 // Azimuth.
  fout.close();
}


/*--------------------------------------------------------------------
  Run GNUplot script to plot 2D data form the file.
--------------------------------------------------------------------*/
void GNUplot2Dsimple(
  const std::string &short_name,  // Data file name without extension.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const double &y_min,            // Minimal value of function.
  const double &y_max,            // Maximal value of function.
  const std::string &title,       // Title of the figure.
  const std::string &x_title,     // Title of X-axis.
  const std::string &y_title,     // Title of Y-axis.
  const bool &is_eps)             // Flag to plot .EPS figure.
{
  double ymn = y_min;
  double ymx = y_max;
  if( (ymn > 0.0) && (ymn < 0.1*y_max) )
    ymn = 0.0;
  if( (ymx < 0.0) && (ymx > 0.1*y_min) )
    ymx = 0.0;

  std::string file_name = "plot.plt";
  std::ofstream fout(file_name.c_str(), std::ios::out);

  if(is_eps) {
    fout << "set terminal postscript eps enhanced font \",24\"\n";
    fout << "set output \"" << short_name << ".eps\"\n";
  } else {
    fout << "set terminal pngcairo enhanced size "
      << X_SIZE << "," << int(X_SIZE*0.75) << " font \",24\"\n";
    fout << "set output \"" << short_name << ".png\"\n";
  }

  fout << "set title \"" << title << "\"\n";
  fout << "set xrange[" << x_min << ":" << x_max << "]\n";
  fout << "set yrange[" << ymn << ":" << ymx << "]\n";
  fout << "set xlabel \"" << x_title << "\"\n";
  fout << "set ylabel \"" << y_title << "\"\n";
  fout << "set nokey\n";
  fout << "set bmargin 4\n";
  fout << "set mxtics 2\n";
  fout << "set mytics 2\n";
  fout << "set grid xtics ytics mxtics mytics\n";
  fout << "set format y \"%2.1t✕10^{%L}\"\n";
  fout << "set style line 1 lt 1 lc 1 lw 4\n";
  fout << "plot \"" << short_name
    << ".dat\" u 1:2 w l smooth csplines ls 1\n";
  fout.close();

  std::string command = "gnuplot " + file_name;
  system(command.c_str());
  command = "rm " + file_name;
  system(command.c_str());
}


/*--------------------------------------------------------------------
  Run GNUplot script to 2D polar plot form the data file.
--------------------------------------------------------------------*/
void GNUplot2Dpolar(
  const std::string &short_name,  // Data file name without extension.
  const double &f_max,            // Maximal value of function.
  const std::string &title,       // Title of the figure.
  const bool &is_eps)             // Flag to plot .EPS figure.
{
  std::string file_name = "plot.plt";
  std::ofstream fout(file_name.c_str(), std::ios::out);

  if(is_eps) {
    fout << "set terminal postscript eps enhanced font \",24\"\n";
    fout << "set output \"" << short_name << ".eps\"\n";
  } else {
    fout << "set terminal pngcairo size "
      << X_SIZE << "," << X_SIZE << " font \",24\"\n";
    fout << "set output \"" << short_name << ".png\"\n";
  }

  fout << "set title \"" << title << "\"\n";
  fout << "set polar\n";
  fout << "set angles degrees\n";
  fout << "set style line 10 lt 1 lc 0 lw 0.3\n";
  fout << "set grid polar 45\n";
  fout << "set grid ls 10\n";
  fout << "set xrange [" << -f_max << ":" << f_max << "]\n";
  fout << "set yrange [" << -f_max << ":" << f_max << "]\n";
  fout << "set xtics axis\n";
  fout << "set ytics axis\n";
  fout << "unset xtics\n";
  fout << "unset ytics\n";
  fout << "set size square\n";
  fout << "unset border\n";
  fout << "set_label(x, text) = sprintf";
  fout << "(\"set label \'%s\' at (" << f_max << "*cos(%f)), ";
  fout << "(" << f_max << "*sin(%f)) center\", text, x, x)\n";
  fout << "eval set_label(0, \"0\")\n";
  fout << "eval set_label(45, \"45\")\n";
  fout << "eval set_label(90, \"90\")\n";
  fout << "eval set_label(135, \"135\")\n";
  fout << "eval set_label(180, \"180\")\n";
  fout << "eval set_label(225, \"225\")\n";
  fout << "eval set_label(270, \"270\")\n";
  fout << "eval set_label(315, \"315\")\n";
  fout << "unset key\n";
  fout << "set style line 1 lt 1 lc 1 lw 4\n";
  fout << "plot \"" << short_name
    << ".dat\" u 1:2 w l ls 1\n";
  fout.close();

  std::string command = "gnuplot " + file_name;
  system(command.c_str());
  command = "rm " + file_name;
  system(command.c_str());
}


/*--------------------------------------------------------------------
  Run GNUplot script to plot 3D cylinder of 2D polar plots of some
  function "func(phi, a)" with a scan parameter, a, varying along
  the cylinder axis.
--------------------------------------------------------------------*/
void GNUplot2DpolarScan(
  const std::string &short_name,  // Data file name without extension.
  const double &r_min,            // Result minimal value of function.
  const double &r_max,            // Result maximal value of function.
  const double &a_min,            // Minimal value of the parameter.
  const double &a_max,            // Maximal value of the parameter.
  const std::string &title)       // Title of the figure.
{
  std::string file_name = "plot.plt";
  std::ofstream fout(file_name.c_str(), std::ios::out);

  fout << "set terminal pngcairo size " 
    << X_SIZE << "," << X_SIZE << " font \",24\"\n";

  fout << "set title \"" << title << "\"\n";
  fout << "set angles degrees\n";
  fout << "set style line 1 lc rgb \'#606060\' lt 1 lw 1\n";
  fout << "set pm3d depthorder hidden3d 1\n";
  fout << "set hidden3d\n";
  fout << "set style fill transparent solid 0.75\n";
  // fout << "set palette rgb 9,9,3\n";
  fout << "set palette rgbformulae 22,13,-31\n";
  fout << "unset colorbox\n";
  fout << "unset key\n";
  // fout << "set border 4095\n";
  fout << "unset border\n";
  fout << "unset tics\n";
  fout << "set ticslevel 0\n";
  fout << "set isosamples 50,100\n";

  fout << "set style line 2 lc rgb \'#404040\' lt 1 lw 3\n";
  fout << "set style arrow 1 nohead ls 2\n";

  double scl = 1.3;
  fout << "set xrange [" << scl*r_min << ":" << scl*r_max << "]\n";
  fout << "set yrange [" << scl*r_min << ":" << scl*r_max << "]\n";
  fout << "set zrange [" << a_min - 0.2*(a_max + a_min) << ":" 
    << a_max + 0.2*(a_max + a_min) << "]\n";

  // fout << "set xyplane at " << 0.5*(a_max + a_min) << "\n";

  fout << "set_label(x, text, z, s) = sprintf";
  fout << "(\"set label \'%s\' at (" << 1.2*r_max << "*cos(%f)), ";
  fout << "(" << 1.2*r_max 
    << "*sin(%f)), %f center %s\", text, x, x, z, s)\n";

  fout << "set style line 3 lc rgb \'#660000\' lt 1 lw 8\n";
  fout << "set style arrow 2 nohead ls 3\n";
  fout << "set arrow from 0, 0, " << a_min - 0.4*(a_max + a_min);
  fout << " to 0, 0, " << a_max << " as 2\n";
  fout << "set arrow from 0, 0, " << a_max << " to 0, 0, ";
  fout << a_max + 0.4*(a_max + a_min) << " as 2 front\n";

  fout << "set_line(x, z, s) = sprintf";
  fout << "(\"set arrow from (" << -1.1*r_max << "*cos(%f)), ";
  fout << "(" << -1.1*r_max << "*sin(%f)), %f to (" << 1.1*r_max 
    << "*cos(%f)), " << "(" << 1.1*r_max;
  fout << "*sin(%f)), %f as 1 %s\", x, x, z, x, x, z, s)\n";
  double z = a_min;
  fout << "eval set_line(0, " << z << ", \"\")\n";
  fout << "eval set_line(45, " << z << ", \"\")\n";
  fout << "eval set_line(90, " << z << ", \"\")\n";
  fout << "eval set_line(135, " << z << ", \"\")\n";
  fout << "eval set_label(0, \"0\", " << z << ", \"\")\n";
  fout << "eval set_label(45, \"45\", " << z << ", \"\")\n";
  fout << "eval set_label(90, \"90\", " << z << ", \"\")\n";
  fout << "eval set_label(135, \"135\", " << z << ", \"\")\n";
  fout << "eval set_label(180, \"180\", " << z << ", \"\")\n";
  fout << "eval set_label(225, \"225\", " << z << ", \"\")\n";
  fout << "eval set_label(270, \"270\", " << z << ", \"\")\n";
  fout << "eval set_label(315, \"315\", " << z << ", \"\")\n";
  z = a_max;
  fout << "eval set_line(0, " << z << ", \"front\")\n";
  fout << "eval set_line(45, " << z << ", \"front\")\n";
  fout << "eval set_line(90, " << z << ", \"front\")\n";
  fout << "eval set_line(135, " << z << ", \"front\")\n";
  fout << "eval set_label(0, \"0\", " << z << ", \"front\")\n";
  fout << "eval set_label(45, \"45\", " << z << ", \"front\")\n";
  fout << "eval set_label(90, \"90\", " << z << ", \"front\")\n";
  fout << "eval set_label(135, \"135\", " << z << ", \"front\")\n";
  fout << "eval set_label(180, \"180\", " << z << ", \"front\")\n";
  fout << "eval set_label(225, \"225\", " << z << ", \"front\")\n";
  fout << "eval set_label(270, \"270\", " << z << ", \"front\")\n";
  fout << "eval set_label(315, \"315\", " << z << ", \"front\")\n";

  /* fout << "set output \"" << short_name << "-rp0.png\"\n";
  fout << "set view 60,0,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n"; */

  fout << "set output \"" << short_name << "-rp1.png\"\n";
  fout << "set view 60,72,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout << "set output \"" << short_name << "-rp2.png\"\n";
  fout << "set view 60,144,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout << "set output \"" << short_name << "-rp3.png\"\n";
  fout << "set view 60,216,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout << "set output \"" << short_name << "-rp4.png\"\n";
  fout << "set view 60,288,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout.close();

  std::string command = "gnuplot " + file_name;
  system(command.c_str());
  command = "rm " + file_name;
  system(command.c_str());
}


/*--------------------------------------------------------------------
  Run GNUplot script to plot 3D colormap for the data form the file.
--------------------------------------------------------------------*/
void GNUplot3Dcolormap(
  const std::string &short_name,  // Data file name without extension.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const double &y_min,            // Minimal value of Y-argument.
  const double &y_max,            // Maximal value of Y-argument.
  const double &z_min,            // Minimal value of function.
  const double &z_max,            // Maximal value of function.
  const std::string &title,       // Title of the figure.
  const std::string &x_title,     // Title of X-axis.
  const std::string &y_title,     // Title of Y-axis.
  const double &ratio,            // Aspect ratio.
  const bool &is_log_scale)       // Flag to plot data in LOG scale.
{
  std::string file_name = "plot.plt";
  std::ofstream fout(file_name.c_str(), std::ios::out);
  fout << "set terminal pngcairo size " 
    << X_SIZE << "," << int(X_SIZE*ratio) << " font \",22\"\n";
  fout << "set size ratio " << ratio << "\n";
  fout << "set pm3d map\n"; // fout << "set view map\n";
  fout << "set palette rgbformulae 22,13,-31\n";
  fout << "unset clabel\n";
  fout << "unset key\n";
  fout << "set output \"" << short_name << ".png\"\n";
  if(is_log_scale)
    fout << "set title \"Log. scale colors. " << title << "\"\n";
  else
    fout << "set title \"" << title << "\"\n";
  fout << "set xlabel \"" << x_title << "\"\n";
  fout << "set ylabel \"" << y_title << "\"\n";
  fout << "set xrange[" << x_min << ":" << x_max << "]\n";
  fout << "set yrange[" << y_min << ":" << y_max << "]\n";
  if(is_log_scale)
    fout << "set logscale cb \n";
  fout << "splot \"" << short_name << ".dat\" with pm3d\n";
  fout.close();

  std::string command = "gnuplot " + file_name;
  system(command.c_str());
  command = "rm " + file_name;
  system(command.c_str());
}


/*--------------------------------------------------------------------
  Run GNUplot script to plot 3D surface for the data form the file.
--------------------------------------------------------------------*/
void GNUplot3Dsurface(
  const std::string &short_name,  // Data file name without extension.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const double &y_min,            // Minimal value of Y-argument.
  const double &y_max,            // Maximal value of Y-argument.
  const double &z_min,            // Minimal value of function.
  const double &z_max,            // Maximal value of function.
  const std::string &title,       // Title of the figure.
  const std::string &x_title,     // Title of X-axis.
  const std::string &y_title,     // Title of Y-axis.
  const double &ratio,            // Aspect ratio.
  const bool &is_log_scale)       // Flag to plot data in LOG scale.
{
  // int ntics = 3;    // Number of tics

  std::string file_name = "plot.plt";
  std::ofstream fout(file_name.c_str(), std::ios::out);
  fout << "set terminal pngcairo size " 
    << X_SIZE << "," << int(X_SIZE*ratio) << " font \",22\"\n";
  fout << "set size ratio " << ratio << "\n";
  // fout << "set view 20,340 #HBB: ,1,2\n";
  fout << "set view 20,340,1,2\n";
  fout << "set samples 100, 100\n";
  fout << "set isosamples 100, 100\n";
  fout << "set style line 1 lc rgb \'#157545\' "
    << "lt 1 lw 0.25 # --- green\n";
  fout << "set pm3d depthorder hidden3d 1\n";
  fout << "set hidden3d\n";
  fout << "set style fill transparent solid 0.75\n";
  fout << "set palette rgb 9,9,3\n";
  fout << "unset colorbox\n";
  fout << "unset key\n";
  fout << "unset clabel\n";
  fout << "set output \"" << short_name << ".png\"\n";
  if(is_log_scale)
    fout << "set title \"Log. scale colors. " << title << "\"\n";
  else
    fout << "set title \"" << title << "\"\n";
  fout << "set xlabel \"" << x_title << "\"\n";
  fout << "set ylabel \"" << y_title << "\"\n";
  fout << "set ticslevel 0\n";
  // fout << "set ztics " << (z_max - z_min)/ntics << "\n";
  fout << "set xrange[" << x_min << ":" << x_max << "]\n";
  fout << "set yrange[" << y_min << ":" << y_max << "]\n";
  if(is_log_scale)
    fout << "set logscale cb \n";
  fout << "splot \"" << short_name << ".dat\" with pm3d\n";
  fout.close();

  std::string command = "gnuplot " + file_name;
  system(command.c_str());
  command = "rm " + file_name;
  system(command.c_str());
}


/*--------------------------------------------------------------------
  Run GNUplot script to plot 3D colored surface
  for the data form the file.
--------------------------------------------------------------------*/
void GNUplot3DsurfaceColor(
  const std::string &short_name,  // Data file name without extension.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const double &y_min,            // Minimal value of Y-argument.
  const double &y_max,            // Maximal value of Y-argument.
  const double &z_min,            // Minimal value of function.
  const double &z_max,            // Maximal value of function.
  const std::string &title,       // Title of the figure.
  const std::string &x_title,     // Title of X-axis.
  const std::string &y_title,     // Title of Y-axis.
  const double &ratio,            // Aspect ratio.
  const bool &is_log_scale)       // Flag to plot data in LOG scale.
{
  // int ntics = 3;    // Number of tics

  std::string file_name = "plot.plt";
  std::ofstream fout(file_name.c_str(), std::ios::out);
  fout << "set terminal pngcairo size " 
    << X_SIZE << "," << int(X_SIZE*ratio) << " font \",22\"\n";
  fout << "set size ratio " << ratio << "\n";
  // fout << "set view 20,340 #HBB: ,1,2\n";
  fout << "set view 20,340,1,2\n";
  fout << "set samples 100, 100\n";
  fout << "set isosamples 100, 100\n";
  fout << "set palette rgbformulae 22,13,-31\n";
  fout << "unset clabel\n";
  fout << "unset key\n";
  fout << "set output \"" << short_name << ".png\"\n";
  if(is_log_scale)
    fout << "set title \"Log. scale colors. " << title << "\"\n";
  else
    fout << "set title \"" << title << "\"\n";
  fout << "set xlabel \"" << x_title << "\"\n";
  fout << "set ylabel \"" << y_title << "\"\n";
  fout << "set ticslevel 0\n";
  // fout << "set ztics " << (z_max - z_min)/ntics << "\n";
  fout << "set xrange[" << x_min << ":" << x_max << "]\n";
  fout << "set yrange[" << y_min << ":" << y_max << "]\n";
  if(is_log_scale)
    fout << "set logscale cb \n";
  fout << "splot \"" << short_name << ".dat\" with pm3d\n";
  fout.close();

  std::string command = "gnuplot " + file_name;
  system(command.c_str());
  command = "rm " + file_name;
  system(command.c_str());
}


/*--------------------------------------------------------------------
  Run GNUplot script to 3D polar plot form the data file.
--------------------------------------------------------------------*/
void GNUplot3Dpolar(
  const std::string &short_name,  // Data file name without extension.
  const double &r_min,            // Result minimal value of function.
  const double &r_max,            // Result maximal value of function.
  const std::string &title)       // Title of the figure.
{
  std::string file_name = "plot.plt";
  std::ofstream fout(file_name.c_str(), std::ios::out);

  fout << "set terminal pngcairo size " 
    << X_SIZE << "," << X_SIZE << " font \",24\"\n";

  fout << "set style line 1 lc rgb \'#157545\' "
    << "lt 1 lw 1 # --- green\n";
  fout << "set pm3d depthorder hidden3d 1\n";
  fout << "set hidden3d\n";
  fout << "set style fill transparent solid 0.75\n";
  fout << "set palette rgb 9,9,3\n";
  fout << "unset colorbox\n";
  fout << "unset key\n";
  // fout << "set border 4095\n";
  // fout << "unset border\n";
  // fout << "unset tics\n";
  fout << "set ticslevel 0\n";
  fout << "set isosamples 50,100\n";

  double scl = 1.1;
  fout << "set xrange [" << scl*r_min << ":" << scl*r_max << "]\n";
  fout << "set yrange [" << scl*r_min << ":" << scl*r_max << "]\n";
  fout << "set zrange [" << scl*r_min << ":" << scl*r_max << "]\n";

  fout << "set xyplane at 0\n";
  fout << "set zeroaxis linetype 1 linewidth 2.5\n";
  fout << "set xlabel \"X\"\n";
  fout << "set ylabel \"Y\"\n";
  fout << "set zlabel \"Z\"\n";

  /* fout << "set output \"" << short_name << "-rp0.png\"\n";
  fout << "set view 60,0,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n"; */

  fout << "set output \"" << short_name << "-rp1.png\"\n";
  fout << "set view 60,72,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout << "set output \"" << short_name << "-rp2.png\"\n";
  fout << "set view 60,144,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout << "set output \"" << short_name << "-rp3.png\"\n";
  fout << "set view 60,216,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout << "set output \"" << short_name << "-rp4.png\"\n";
  fout << "set view 60,288,1.5,1\n";
  fout << "splot \"" << short_name << ".dat\" w pm3d\n";

  fout.close();

  std::string command = "gnuplot " + file_name;
  system(command.c_str());
  command = "rm " + file_name;
  system(command.c_str());
}


/*--------------------------------------------------------------------
  Plot 2D graph of some function "func".
--------------------------------------------------------------------*/
void plot2Dsimple(
  double (*func)(double),       // Function to plot.
  const double &x_min,          // Minimal value of argument.
  const double &x_max,          // Maximal value of argument.
  const int &x_num,             // Number of argument points.
  const std::string &name,      // Data identification name.
  const std::string &title,     // Title of the figure.
  const std::string &x_title,   // Title of X-axis.
  const std::string &y_title,   // Title of Y-axis.
  const bool &is_eps)           // Flag to plot additional EPS figure.
{
  double y_min, y_max;

  // Print the data to plot.
  print2Ddata(func, x_min, x_max, x_num, name, y_min, y_max);

  // Plot using GNUplot script.
  if(is_eps)
    GNUplot2Dsimple(name, x_min, x_max, y_min, y_max,
      title, x_title, y_title, true);
  GNUplot2Dsimple(name, x_min, x_max, y_min, y_max,
    title, x_title, y_title, false);
}


/*--------------------------------------------------------------------
  Plot 2D polar plot of some function "func".
  ARGUMENT SHOULD BE IN DEGREES!
--------------------------------------------------------------------*/
void plot2Dpolar(
  double (*func)(double),       // Function to plot.
  const int &x_num,             // Number of argument points.
  const std::string &name,      // Data identification name.
  const std::string &title,     // Title of the figure.
  const bool &is_eps)           // Flag to plot additional EPS figure.
{
  double f_min, f_max;

  // Print the data to plot.
  print2Ddata(func, 0.0, 360.0, x_num, name, f_min, f_max);

  // Plot using GNUplot script.
  if(is_eps)
    GNUplot2Dpolar(name, f_max, title, true);
  GNUplot2Dpolar(name, f_max, title, false);
}


/*--------------------------------------------------------------------
  Plot cylinder of 2D polar plots of some function "func(phi, a)"
  with a scan parameter, a, varying along the cylinder axis.
  ARGUMENT phi SHOULD BE IN DEGREES! a - is the parameter.
--------------------------------------------------------------------*/
void plot2DpolarScan(
  double (*func)(double, double), // Function to plot.
  const int &x_num,             // Number of argument points.
  const double &a_min,          // Minimal value of the parameter.
  const double &a_max,          // Maximal value of the parameter.
  const int &a_num,             // Number of points for the parameter.
  const std::string &name,      // Data identification name.
  const std::string &title)     // Title of the figure.
{
  double r_min, r_max;

  // Print the data to plot.
  print2DdataScanXYZ(func, x_num,
    a_min, a_max, a_num, name, r_min, r_max);

  // Plot using GNUplot script.
  GNUplot2DpolarScan(name, r_min, r_max,
    a_min, a_max, title);
}


/*--------------------------------------------------------------------
  Plot 3D colormap for some function "func".
--------------------------------------------------------------------*/
void plot3Dcolormap(
  double (*func)(double, double), // Function to plot.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const int &x_num,               // Number of X-argument points.
  const double &y_min,            // Minimal value of Y-argument.
  const double &y_max,            // Maximal value of Y-argument.
  const int &y_num,               // Number of Y-argument points.
  const std::string &name,        // Data identification name.
  const std::string &title,       // Title of the figure.
  const std::string &x_title,     // Title of X-axis.
  const std::string &y_title,     // Title of Y-axis.
  const bool &is_log_scale)       // Flag to plot data in LOG scale.
{
  double z_min, z_max;

  // Print the data to plot.
  print3Ddata(func, x_min, x_max, x_num,
    y_min, y_max, y_num, name, is_log_scale, z_min, z_max);

  // Plot using GNUplot script.
  double ratio = fabs( (y_max - y_min)/(x_max - x_min) );
  GNUplot3Dcolormap(name, x_min, x_max, y_min, y_max, z_min, z_max,
    title, x_title, y_title, ratio, is_log_scale);
}


/*--------------------------------------------------------------------
  Plot 3D surface for some function "func".
--------------------------------------------------------------------*/
void plot3Dsurface(
  double (*func)(double, double), // Function to plot.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const int &x_num,               // Number of X-argument points.
  const double &y_min,            // Minimal value of Y-argument.
  const double &y_max,            // Maximal value of Y-argument.
  const int &y_num,               // Number of Y-argument points.
  const std::string &name,        // Data identification name.
  const std::string &title,       // Title of the figure.
  const std::string &x_title,     // Title of X-axis.
  const std::string &y_title,     // Title of Y-axis.
  const bool &is_log_scale)       // Flag to plot data in LOG scale.
{
  double z_min, z_max;

  // Print the data to plot.
  print3Ddata(func, x_min, x_max, x_num,
    y_min, y_max, y_num, name, is_log_scale, z_min, z_max);

  // Plot using GNUplot script.
  double ratio = fabs( (y_max - y_min)/(x_max - x_min) );
  GNUplot3Dsurface(name, x_min, x_max, y_min, y_max, z_min, z_max,
    title, x_title, y_title, ratio, is_log_scale);
}


/*--------------------------------------------------------------------
  Plot 3D colored surface for some function "func".
--------------------------------------------------------------------*/
void plot3DsurfaceColor(
  double (*func)(double, double), // Function to plot.
  const double &x_min,            // Minimal value of X-argument.
  const double &x_max,            // Maximal value of X-argument.
  const int &x_num,               // Number of X-argument points.
  const double &y_min,            // Minimal value of Y-argument.
  const double &y_max,            // Maximal value of Y-argument.
  const int &y_num,               // Number of Y-argument points.
  const std::string &name,        // Data identification name.
  const std::string &title,       // Title of the figure.
  const std::string &x_title,     // Title of X-axis.
  const std::string &y_title,     // Title of Y-axis.
  const bool &is_log_scale)       // Flag to plot data in LOG scale.
{
  double z_min, z_max;

  // Print the data to plot.
  print3Ddata(func, x_min, x_max, x_num,
    y_min, y_max, y_num, name, is_log_scale, z_min, z_max);

  // Plot using GNUplot script.
  double ratio = fabs( (y_max - y_min)/(x_max - x_min) );
  GNUplot3DsurfaceColor(name, x_min, x_max, y_min, y_max, z_min, z_max,
    title, x_title, y_title, ratio, is_log_scale);
}


/*--------------------------------------------------------------------
  Plot 3D polar surface for function "func(theta, phi)".
  ANGLES SHOULD BE IN DEGREES! theta = 0..180; phi = 0..360.
--------------------------------------------------------------------*/
void plot3Dpolar(
  double (*func)(double, double), // Function to plot.
  const int &aa_num,              // Number of points for azimuth.
  const int &pa_num,              // Number of points for polar angle.
  const std::string &name,        // Data identification name.
  const std::string &title)       // Title of the figure.
{
  double r_min, r_max;

  // Print the data to plot.
  print3DpolarXYZ(func, aa_num, pa_num, name, r_min, r_max);

  // Plot using GNUplot script.
  GNUplot3Dpolar(name, r_min, r_max, title);
}


/*--------------------------------------------------------------------
  Test 2D function.
--------------------------------------------------------------------*/
double my_func2(double x)
{
  double d2r = M_PI/180.0;
  return pow(cos(x*d2r), 2);
}


/*--------------------------------------------------------------------
  Test 2D function with parameter.
--------------------------------------------------------------------*/
double my_func2scan(double x, double a)
{
  double d2r = M_PI/180.0;
  return pow(cos(x*d2r + a), 2);
}


/*--------------------------------------------------------------------
  Test 3D function.
--------------------------------------------------------------------*/
double my_func3(double x, double y)
{
  double d2r = M_PI/180.0;
  return 1.0 + 2.0*cos(2.0*x*d2r);
}


/*********************************************************************
  MAIN
*********************************************************************/
int main(int argc, char **argv)
{
  // Plot 2D graph of some function "func".
  plot2Dsimple(my_func2, 0.0, 360.0, 180,
    "Test2Dsimple", "Test for simple 2D plot ",
    "Theta [Degrees]", "Test function", true);

  /* Plot 2D ploar plot of some function "func".
     ARGUMENT SHOULD BE IN DEGREES!   */
  plot2Dpolar(my_func2, 180, "Test2Dpolar",
   "Test for 2D polar plot", true);

  /* Plot cylinder of 2D polar plots of some function "func(phi, a)"
     with a scan parameter, a, varying along the cylinder axis.
     ARGUMENT phi SHOULD BE IN DEGREES! a - is the parameter. */
  plot2DpolarScan(my_func2scan, 90, 0.0, 1.57, 16, "Test2DpolarScan",
    "Test for cylinder of 2D polar plots with a scan parameter");

  // Plot 3D colormap for some function "func".
  plot3Dcolormap(my_func3, 0.0, 360.0, 180, 0.0, 180.0, 90,
    "Test3Dcolormap", "Test plot for 3D",
    "Phi [Degrees]", "Theta [Degrees]", false);

  // Plot 3D surface for some function "func".
  plot3Dsurface(my_func3, 0.0, 360.0, 180, 0.0, 180.0, 90,
    "Test3Dsurface", "Test plot for 3D",
    "Phi [Degrees]", "Theta [Degrees]", false);

  // Plot 3D colored surface for some function "func".
  plot3DsurfaceColor(my_func3, 0.0, 360.0, 180, 0.0, 180.0, 90,
    "Test3DsurfaceColor", "Test plot for 3D",
    "Phi [Degrees]", "Theta [Degrees]", false);

  /* Plot 3D polar surface for function "func(theta, phi)".
     ANGLES SHOULD BE IN DEGREES! theta = 0..180; phi = 0..360. */
  plot3Dpolar(my_func3, 90, 45,
    "Test3Dpolar", "Test for 3D polar plot");

  return 0;
}


//====================================================================
