// utils.cpp
//
// James Folberth
// Spring 2015

#include "utils.hpp"

using namespace std;

//////////////////
// Vector stuff //
//////////////////

// print a vector as a column vector
template<typename T>
void print_vector(valarray<T>& v) {
   for (size_t i = 0; i < v.size(); ++i) {
      cout << setw(_PRINT_VECTOR_WIDTH_) << setfill(' ') 
           << setprecision(_PRINT_VECTOR_PREC_)
           << static_cast<double>(v[i]) << endl;
   }
   cout << endl;
}

template<typename T>
void print_vector(const valarray<T>& v) {
   for (size_t i = 0; i < v.size(); ++i) {
      cout << setw(_PRINT_VECTOR_WIDTH_) << setfill(' ') 
           << setprecision(_PRINT_VECTOR_PREC_)
           << static_cast<double>(v[i]) << endl;
   }
   cout << endl;
}

template<typename T>
void print_vector(vector<T>& v) {
   for (auto it = v.begin(); it != v.end(); ++it) {
      cout << setw(_PRINT_VECTOR_WIDTH_) << setfill(' ') 
           << setprecision(_PRINT_VECTOR_PREC_)
           << static_cast<double>(*it) << endl;
   }
   cout << endl;
}

template<typename T>
void print_vector(const vector<T>& v) {
   for (auto it = v.begin(); it != v.end(); ++it) {
      cout << setw(_PRINT_VECTOR_WIDTH_) << setfill(' ') 
           << setprecision(_PRINT_VECTOR_PREC_)
           << static_cast<double>(*it) << endl;
   }
   cout << endl;
}


// return a vector with uniform[0,1] random entries
template<typename T>
valarray<T> rand_vec(const unsigned m, const T low, const T high) {
   random_device rd;
   std::mt19937 e2(rd());
   uniform_real_distribution<T> val_dist(low,high);

   valarray<T> v(m);

   for (size_t i=0; i<m; ++i) {
      v[i] = val_dist(e2);
   }
   
   return v;
}

// vector norms
// p = 0 is infinity norm
// p = 1 is l^1 norm
// p = 2 is l^2 norm; like BLAS dnrm2
template<typename T>
T norm(const valarray<T>& v, const unsigned p) {
   T res = 0, scale = 0., absvi, ssq = 1., tmp;

   if ( v.size() < 1 ) {
      return static_cast<T>(0.);
   }

   else if ( v.size() == 1 ) {
      return abs(v[0]);
   }

   else {
      switch (p) {
         case 0:
            //for (size_t i = 0; i < v.size(); ++i) {
            //   if ( res < abs(v[i]) )
            //      res = abs(v[i]);
            //}
            //return res;
            for (auto it = begin(v); it != end(v); ++it) {
               if ( res < abs(*it) )
                  res = abs(*it);
            }
            return res;
            //return MAX(abs(v.min()), abs(v.max()));

         case 1:
            for (size_t i = 0; i < v.size(); ++i) {
               res += abs(v[i]);
            }
            return res;
            //for ( auto it = begin(v); it != end(v); ++it)
            //   res += abs(*it);
            //return res;

         case 2:
            for (size_t i = 0; i < v.size(); ++i) {
               if ( v[i] != 0. ) {
                  absvi = abs(v[i]);
                  if ( scale < absvi ) {
                     tmp = scale/absvi;
                     ssq = 1. + ssq*tmp*tmp;
                     scale = absvi;
                  }
                  else {
                     tmp = absvi/scale;
                     ssq += tmp*tmp;
                  }
               }
            }
            res = scale*sqrt(ssq);
            return res;

         default:
            cerr << "error: utils:norm(vector): unsupported p-norm: p = " << p 
                 << endl;
            exit(-1);
      }  
   }
}

// Discrete L^2 norms
template<typename T>
T dl2norm(const valarray<T>& v, const unsigned nx) {
   return pow(static_cast<T>(nx+1), -0.5)*norm(v,2);
}

template<typename T>
T dl2norm(const valarray<T>& v, const unsigned nx, const unsigned ny) {
   return pow(static_cast<T>(nx+1),-0.5) 
      * pow(static_cast<T>(ny+1),-0.5)
      * norm(v,2);
}

//////////
// Misc //
//////////

// Integer pow
int pow(int b, int e) {
   if ( e == 0 ) return 1;
   else if ( e == 1 ) return b;
   else if ( e < 0 ) {
      cerr << "utils.cpp:pow<int>: negative exponent not handled" << endl;
      exit(-1);
   }
   else return pow(b, e-1);
}

unsigned pow(unsigned b, unsigned e) {
   if ( e == 0 ) return 1;
   else if ( e == 1 ) return b;
   else if ( e < 0 ) {
      cerr << "utils.cpp:pow<unsigned>: negative exponent not handled" << endl;
      exit(-1);
   }
   else return pow(b, e-1);
}


/////////////
// Sorting //
/////////////
//void sort_perm(vector<unsigned>& p, vector<unsigned>::iterator v_begin,
//      vector<unsigned>::iterator v_end) {
//
//   // p = {0, 1, 2, ..., num-1}
//   auto num = distance(v_begin, v_end);
//   if ( static_cast<unsigned>(p.size()) < num ) {
//      cerr << "error: utils.cpp:sort_perm: permutation vector p is not large "
//           << "enough!" << endl;
//      exit(-1);
//   }
//   
//   iota(p.begin(), p.begin() + num, 0);
//   sort(p.begin(), p.begin() + num, [&] (unsigned i, unsigned j) -> bool 
//         { return *(v_begin + i) < *(v_begin + j); });
//
//}





// Force instantiation
template void print_vector<double>(valarray<double>& v);
template void print_vector<double>(const valarray<double>& v);
template void print_vector<double>(vector<double>& v);
template void print_vector<double>(const vector<double>& v);

template valarray<double> rand_vec<double>(const unsigned,
                                         const double, const double);
template double norm<double>(const valarray<double>&, const unsigned);
template double dl2norm<double>(const valarray<double>&, const unsigned);
template double dl2norm<double>(const valarray<double>&, const unsigned,
      const unsigned);

template void print_vector<unsigned>(vector<unsigned>& v);
template void print_vector<unsigned>(const vector<unsigned>& v);

template void print_vector<int>(vector<int>& v);
template void print_vector<int>(const vector<int>& v);

