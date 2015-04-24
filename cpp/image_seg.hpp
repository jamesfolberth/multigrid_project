// image_seg.hpp
//
// Spring 2015
// James Folberth, Nathan Heavner, Rachel Tutmaher

#ifndef _IMAGE_SEG_HPP_
#define _IMAGE_SEG_HPP_ 

#include <algorithm>
#include <cmath>
#include <climits>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <list>
#include <valarray>
#include <vector>
#include <string>
#include <chrono> // for algorithm runtime things
#include <fstream>

#include <opencv2/opencv.hpp> // We use OpenCV for image I/O and cv::Mat 

#include "matrix_coo.hpp"
#include "matrix_crs.hpp"
#include "classical_solvers.hpp"
#include "multigrid.hpp"
#include "utils.hpp"

using namespace std;

////////////////////////
// Image Segmentation //
////////////////////////
// {{{

// Class to store all image segmentation variables for a specific level
// The idea here is to create an STL vector of levels
class image_level {
   public:
      matrix_crs<double> A; // coupling matrix
      matrix_crs<double> S; // variance matrix
      matrix_crs<double> L; // weighted boundary length matrix
      matrix_crs<double> V; // area matrix
      matrix_crs<double> G; // boundary length matrix

      valarray<double> Gamma; // saliency vector

      valarray<double> I;   // intensity vector

      // construction
      image_level() = default;

      // destruction, move, copy
      image_level(const image_level& ) = default;
      image_level& operator=(const image_level& ) = default;
      image_level(image_level&& ) = default;
      image_level& operator=(image_level&& ) = default;
      ~image_level() = default;

};


// Class to store the image segmentation parameters
// These parameters are global, in the sense that the shouldn't vary across 
// the levels
class seg_params {
   public:
      unsigned n;       // Size of image: n x n
      
      double alpha;     // top-level intensity scaling factor
      double alpha_til; // coarse-level intensity rescaling factor
      double beta;      // coarse-level variance rescaling factor
      double theta;     // coarsening strength threshold
      double gamma;     // saliency threshold
      double d1;        // sharpening threshold
      unsigned sigma;     // segment detection threshold level
      unsigned rho;       // variance rescaling detection threshold level

      // construction
      seg_params() = default;
      seg_params(double init_alpha, double init_alpha_til, double init_beta,
            double init_theta, double init_gamma, double init_d1,
            double init_sigma, double init_rho) {
         alpha     = init_alpha;
         alpha_til = init_alpha_til;
         beta      = init_beta;
         theta     = init_theta;
         gamma     = init_gamma;
         d1        = init_d1;
         sigma     = init_sigma;
         rho       = init_rho;
      }

      // destruction, move, copy
      seg_params(const seg_params& ) = default;
      seg_params& operator=(const seg_params& ) = default;
      seg_params(seg_params&& ) = default;
      seg_params& operator=(seg_params&& ) = default;
      ~seg_params() = default;

};


// First Level
//////////////
void build_first_level(const valarray<double>& I, const seg_params& params,
      list<image_level>::iterator it);

matrix_crs<double> build_A1(const valarray<double>& I,
      const seg_params& params);


// General Level
////////////////
// Some of these are used on the first level
void build_L(list<image_level>::iterator it, 
      const seg_params& params);

void build_G(list<image_level>::iterator it, 
      const seg_params& params);

matrix_crs<double> build_interp(const matrix_crs<double>& A,
      const vector<unsigned>& C, unsigned M, unsigned M_next);

matrix_crs<double> build_scaled_interp(const matrix_crs<double>& P);

matrix_crs<double> coarse_variance(const matrix_crs<double>& Sf, 
      const valarray<double>& Sc);

void rescale_coarse_coupling(matrix_crs<double>& A, 
      const list<image_level>::iterator it, const unsigned l_next,
      const seg_params& params);

// V-cycle
//////////
matrix_crs<double> image_vcycle(unsigned l, unsigned M, 
      list<image_level>& levels, list<image_level>::iterator it,
      seg_params& params);


// Graph Coarsening
///////////////////
vector<unsigned> coarsen_AMG(const list<image_level>::iterator it, 
      const seg_params& params);

vector<unsigned> strongly_influenced_by_j(const matrix_crs<double>& A_bar,
      const vector<unsigned>& T, const unsigned j);

vector<unsigned> strongly_influenced_by_j_trans(const matrix_crs<double>& A_bar_trans,
      const vector<unsigned>& T, const unsigned j);

vector<unsigned> strongly_influence_k(const matrix_crs<double>& A_bar,
      const vector<unsigned>& T, const unsigned k);


// }}}


///////////////////////
// Image Seg Drivers //
///////////////////////
// {{{

matrix_crs<double> image_seg(const cv::Mat& img, seg_params& params);

// Overloaded
matrix_crs<double> image_seg(const cv::Mat& img);

// Helpers
void set_params(seg_params& params, const double alpha, 
      const double alpha_til, const double beta, const double theta,
      const double gamma, const double d1, const unsigned sigma,
      const unsigned rho);

void assign_uniquely(matrix_crs<double>& U);

vector<valarray<double>> mat_to_vecs(const matrix_crs<double>& U);

// }}}


///////////////////
// OpenCV things //
///////////////////
// {{{

cv::Mat load_image(const string& img_filename, const unsigned mode = 0);

valarray<double> image_to_intensity(cv::Mat img, seg_params& params);

cv::Mat intensity_to_image1(const valarray<double>& I, const unsigned n);

void write_seg_images(const cv::Mat& img, const matrix_crs<double>& U,
      const string& filename_base, const unsigned mode = 1);

// }}}


//////////////////////
// Algorithm Timing //
//////////////////////
// {{{

void runtime_scaling(const cv::Mat& img, seg_params& params, const string& out_filename);

// }}}


#endif
