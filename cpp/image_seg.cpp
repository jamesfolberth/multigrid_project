// image_seg.cpp
//
// Spring 2015
// James Folberth, Nathan Heavner, Rachel Tutmaher

#include "image_seg.hpp"

////////////////////////
// Image Segmentation //
////////////////////////
// {{{

// First Level
//////////////
void build_first_level(const valarray<double>& I, const seg_params& params,
      list<image_level>::iterator it) {
// {{{
   // connectivity
   it->A = build_A1(I, params);

   // variance (init to zeros)
   it->S = zeros<double>(I.size(), 1); // TODO: should this be M x M or M x 1?  Ingris 2010 seems contradictory.  Shouldn't matter, either way
   //it->S = zeros<double>(I.size(), I.size());
 
   // weighted boundary length
   build_L(it, params);

   // area 
   it->V = it->A;
   fill(it->V.val.begin(), it->V.val.end(), 1.);

   // boundary length
   build_G(it, params);

   // Build initial saliency vector
   it->Gamma = diag<double>(it->L);
   it->Gamma /= diag<double>(it->G);

   // intensity vector
   it->I = I;

// }}}
}

matrix_crs<double> build_A1(const valarray<double>& I,
      const seg_params& params) {
// {{{
// Build A by rows, since we store in CRS
// The ordering of I for a 5x5 image (e.g. square.png) is
//
// 0   1   2   3   4
// 5   6   7   8   9
// 10  11  12  13  14
// 15  16  17  18  19
// 20  21  22  23  24
//
// NOTE: can call std::vector::insert with an initializer list

   unsigned n = params.n;
   double a = params.alpha;

   vector<unsigned> row_ptr, col_ind;
   vector<double> val;
   row_ptr.resize(I.size()+1);
   col_ind.reserve(4*2+4*(n-2)*3+(n-2)*(n-2)*4);
   val.reserve(4*2+4*(n-2)*3+(n-2)*(n-2)*4);

   if ( _DEBUG_ >= 2) {
      cout << "debug 2: building A1 (nnz = " 
           << 4*2+4*(n-2)*3+(n-2)*(n-2)*4 << ")" << endl;
   }

   row_ptr[0] = 0;
   for (unsigned ind = 0; ind < I.size(); ++ind) {
      row_ptr[ind+1] = row_ptr[ind];

      if ( 0 <= ind && ind < n ) { // on top
         if ( ind % n == 0 ) { // on left side
            // we must have ind == 0, top-left
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " top row, left side "
                    << "I[ind] = " << I[ind] << endl; 
            }

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else if ( (ind + 1) % n == 0 ) { // on right side
            // we must have ind == n-1, top-right
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " top row, right side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else { // in the middle
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " top row, middle side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }
      }

      else if ( (n <= ind) && (ind < (n-1)*n) ) { // in the middle rows
         if ( ind % n == 0 ) { // on left side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " middle row, left side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else if ( (ind + 1) % n == 0 ) { // on right side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " middle row, right side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else { // in the middle
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " middle row, middle side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }
      }

      else { // in bottom row
         if ( ind % n == 0 ) { // on left side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " bottom row, left side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));
         }

         else if ( (ind + 1) % n == 0 ) { // on right side
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " bottom row, right side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));
         }

         else { // in the middle
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " bottom row, middle side "
                    << "I[ind] = " << I[ind] << endl; 
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-n);
            val.push_back(exp(-a*abs(I[ind]-I[ind-n])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+1);
            val.push_back(exp(-a*abs(I[ind]-I[ind+1])));
         }
      }
   }

   return matrix_crs<double>(row_ptr, col_ind, val, I.size(), I.size(), 1);
// }}}
}


// General level
////////////////
void build_L(list<image_level>::iterator it, 
      const seg_params& params) {
// {{{
//     {  -A_{ij}                  i != j
// L = {
//     {   \sum_{k\neq i} A_{ik}   i == j
//
// Equivalently,
// 1) L \gets -A
// 2) L_{ii} \gets -\sum_{k\neq i} L_{ik}
//
// Should be 
//    L = -A; for i in 1:size(L,1); L[i,i] -= sum(L[i,:]); end
//
// or
//    L = -A + diagm(vec(sum(A,2)))
//
// in julia.
//
// JMF 30-03-2015: this was a bottleneck in the code
//      look into std::lower_bound
//      http://www.cplusplus.com/reference/algorithm/lower_bound/
//
//      Is it better to start with empty vectors and build L from scratch?
//      As opposed to inserting things along the way, which is linear in 
//      remaining entries for std::vector.

   it->L = -1.*it->A;
   
   vector<double> val(it->L.m,0.);

   for (unsigned i = 0; i < it->L.m; ++i) {
      for (unsigned jp = it->L.row_ptr[i]; jp < it->L.row_ptr[i+1]; ++jp) {
         val[i] -= it->L.val[jp];
      }
   }

   vector<unsigned> row_ptr(it->L.m + 1), col_ind(it->L.m);
   for (unsigned i = 0; i < it->L.m; ++i) {
      row_ptr[i] = i;
      col_ind[i] = i;
   }
   row_ptr[it->L.m] = it->L.m;

   matrix_crs<double> sdiag(row_ptr, col_ind, val, it->L.m, it->L.n, 1);

   it->L += sdiag;

   //cout << "A = " << endl;
   //it->A.print_full();
   //cout << "L = " << endl;
   //it->L.print_full();
   
// }}}
}

void build_G(list<image_level>::iterator it, 
      const seg_params& params) {
// {{{
//     {  -V_{ij}                  i != j
// G = {
//     {   \sum_{k\neq i} V_{ik}   i == j
//
// Equivalently,
// 1) G \gets -V
// 2) G_{ii} \gets -\sum_{k\neq i} G_{ik}
//
// Should be 
//    G = -V; for i in 1:size(G,1); G[i,i] -= sum(G[i,:]); end
// in julia.
//

   it->G = -1.*it->V;
   
   vector<double> val(it->G.m,0.);

   for (unsigned i = 0; i < it->G.m; ++i) {
      for (unsigned jp = it->G.row_ptr[i]; jp < it->G.row_ptr[i+1]; ++jp) {
         val[i] -= it->G.val[jp];
      }
   }

   vector<unsigned> row_ptr(it->G.m + 1), col_ind(it->G.m);
   for (unsigned i = 0; i < it->G.m; ++i) {
      row_ptr[i] = i;
      col_ind[i] = i;
   }
   row_ptr[it->G.m] = it->G.m;

   matrix_crs<double> sdiag(row_ptr, col_ind, val, it->G.m, it->G.n, 1);

   it->G += sdiag;
 
   //cout << "V = " << endl;
   //it->V.print_full();
   //cout << "G = " << endl;
   //it->G.print_full();
 
// }}}
}

matrix_crs<double> build_interp(const matrix_crs<double>& A,
      const vector<unsigned>& C, unsigned M, unsigned M_next) {
// {{{
// Build the interpolation matrix P based on 
// connectivity matrix A and C-points found from `coarsen_AMG`
//
// It is assumed that C is sorted in increasing order; this should always
// happen due to the construction of C by `coarsen_AMG`.
//

   vector<unsigned> row_ptr(M+1,0), col_ind;
   vector<double> val;
   
   col_ind.reserve(A.col_ind.size());
   val.reserve(A.col_ind.size());

   unsigned C_ind = 0;
   vector<unsigned> A_inds; // vector of indexes for A.col_ind/A.val
   A_inds.reserve( C.size() );
   vector<unsigned> C_inds; // vector of indexes for C
   C_inds.reserve( C.size() );
 
   for (unsigned i = 0; i < M; ++i) {
      row_ptr[i+1] = row_ptr[i];

      if ( C_ind < C.size() && i == C[C_ind] ) { // inject C-point up
         ++row_ptr[i+1];
         col_ind.push_back(C_ind);
         val.push_back(1.);

         ++C_ind;
      }

      else { // interpolate to F-point
         
         // Find all of the C-points where also A_{iC_j} != 0
         // I wasn't able to find a way to use STL's `set_intersection`
         // from <algorithm> for this, since we need the index of A.col_ind 
         // so we can also get A.val.

         // compute ith row sum of A but only over the C-points
         // it is assumed that A's CRS data structure is ordered/sorted and 
         // that the C-points are ordered
         A_inds.clear();
         auto jp_it = A.col_ind.begin() + A.row_ptr[i];
         auto jp_it_end = A.col_ind.begin() + A.row_ptr[i+1];

         C_inds.clear();
         auto j_it = C.begin();
         auto j_it_end = C.end();
 
         double row_sum = 0.;
         while ( jp_it != jp_it_end && j_it != j_it_end ) {
            if ( *jp_it < *j_it ) {
               //++jp_it;
               jp_it = lower_bound(jp_it, jp_it_end, *j_it);
            }

            else if ( *j_it < *jp_it ) {
               //++j_it;
               j_it = lower_bound(j_it, j_it_end, *jp_it); // this search makes a _huge_ difference!
            }

            else {
               row_sum += A.val[jp_it - A.col_ind.begin()];

               A_inds.push_back(jp_it - A.col_ind.begin());
               C_inds.push_back(j_it - C.begin());

               ++jp_it; ++j_it;
            }
         }

         for (unsigned j = 0; j < C_inds.size(); ++j) {
            ++row_ptr[i+1];
            col_ind.push_back(C_inds[j]);
            val.push_back( A.val[A_inds[j]] / row_sum );
         }
 
      }
   }

   return matrix_crs<double>(row_ptr, col_ind, val, M, M_next, 1);
// }}}
}

matrix_crs<double> build_scaled_interp(const matrix_crs<double>& P) {
// {{{
// Build the scaled interpolation matrix Ptil based on the normal
// interpolation matrix P
//
// In Julia,
//  col_sums = sum(P,1); Pt = P; for j in 1:size(Pt,2); Pt[:,j] /= col_sums[j]; end
//

   matrix_crs<double> Ptil = P;
   vector<double> col_sums(P.n,0.);

   //cout << "P = " << endl;
   //P.print_full();
   //cout << P << endl;
   //cout << "Ptil = " << endl;
   //Ptil.print_full();
   
   // compute column sums
   for (unsigned i = 0; i < P.m; ++i) {
      for (unsigned jp = P.row_ptr[i]; jp < P.row_ptr[i+1]; ++jp) {
         col_sums[P.col_ind[jp]] += P.val[jp];
      }
   }

   // scale columns of P by col_sums to form Ptil
   for (unsigned i = 0; i < Ptil.m; ++i) {
      for (unsigned jp = Ptil.row_ptr[i]; jp < Ptil.row_ptr[i+1]; ++jp) {
         Ptil.val[jp] /= col_sums[Ptil.col_ind[jp]];
      }
   }

   //cout << "Ptil = " << endl;
   //Ptil.print_full();
   //cout << Ptil << endl;
 
   return Ptil;
// }}}
}

matrix_crs<double> coarse_variance(const matrix_crs<double>& Sf, 
      const valarray<double>& Sc) {
// {{{
// Build the coarse-level variance matrix S by concatenating the restricted
// fine-level variance matrix Sf and coarse-level variance vector Sc
//
  
   // copy fine-level variance matrix
   matrix_crs<double> S = Sf;

   //cout << "Sf = " << endl;
   //cout << Sf << endl;

   //cout << "Sc = " << endl;
   //print_vector(Sc);

   // Append the contents of Sc to the right-side of S
   unsigned num_added = 0;
   for (unsigned i = 0; i < S.m; ++i) {
      // We don't want to append a bunch of zeros (Sc usually has some, 
      // at least on the first level)
      if ( Sc[i] != 0. ) {
         S.col_ind.insert(S.col_ind.begin() + S.row_ptr[i+1] + num_added, S.n);
         S.val.insert(S.val.begin() + S.row_ptr[i+1] + num_added, Sc[i]);
         ++num_added;
         
      }
      
      S.row_ptr[i+1] += num_added;
   }

   // Adjust size of S
   ++S.n;

   //cout << "S = " << endl;
   //cout << S << endl;

   return S;
// }}}
}

void rescale_coarse_coupling(matrix_crs<double>& A, 
      const list<image_level>::iterator it, const unsigned l_next,
      const seg_params& params) {
// {{{
// Rescale the input matrix A (= A^[r+1]) according to coarse-level intensity
// and also, if l_next >= rho, the multi-level variance

   double a_til = params.alpha_til;

   if ( l_next >= params.rho ) { // rescale using intensity and variance 
   
      double b = params.beta;
   
      // Store the address of S, which we'll use repeatedly later
      const matrix_crs<double> *S = &( next(it)->S );
      
      //cout << "S = " << endl;
      //next(it)->S.print_full();

      //cout << "S->row_ptr = " << endl;
      //print_vector(S->row_ptr);

      for (unsigned i = 0; i < A.m; ++i) {
         
         double Ii = next(it)->I[i];
         unsigned sip_end = S->row_ptr[i+1];

         for (unsigned jp = A.row_ptr[i]; jp < A.row_ptr[i+1]; ++jp) {
            unsigned j = A.col_ind[jp];

            // rescale using coarse-level intensity
            A.val[jp] *= exp(-a_til*abs(Ii - next(it)->I[j]));

            // rescale using multi-level variance
            // Computing the 2-norms is probably slow
            // We're also not utilizing the symmetry of A; that's a factor of 2
            unsigned sip = S->row_ptr[i]; // need to reset sip each time, since we increment it
            unsigned sjp = S->row_ptr[j];
            unsigned sjp_end = S->row_ptr[j+1];

            double nrm = 0.;
            while ( sip < sip_end || sjp < sjp_end ) {
               if ( sjp >= sjp_end || (sip < sip_end && S->col_ind[sip] < S->col_ind[sjp]) ) {
                  nrm += pow(S->val[sip], 2.);
                  ++sip;
               }
               else if ( sip >= sip_end || (sjp < sjp_end && S->col_ind[sjp] < S->col_ind[sip]) ) {
                  nrm += pow(S->val[sjp], 2.);
                  ++sjp;
               }
               else {
                  nrm += pow(S->val[sip]-S->val[sjp], 2.);
                  ++sip; ++sjp;
               }
            }
            nrm = sqrt(nrm);
            
            A.val[jp] *= exp(-b*nrm);

            //cout << "(i,j) = " << "(" << i << "," << j << ")  nrm = " << nrm << endl;
         } 
      }
   }

   else { // don't rescale using multi-level variance
      for (unsigned i = 0; i < A.m; ++i) {
         
         double Ii = next(it)->I[i];
         for (unsigned jp = A.row_ptr[i]; jp < A.row_ptr[i+1]; ++jp) {
            unsigned j = A.col_ind[jp];
            A.val[jp] *= exp(-a_til*abs(Ii - next(it)->I[j]));
         } 
      }
   } 
// }}}
}


// V-cycle
//////////
matrix_crs<double> image_vcycle(unsigned l, unsigned M, 
      list<image_level>& levels, list<image_level>::iterator it,
      seg_params& params) {
// {{{

   if ( l <= params.sigma ) {
      fill(begin(it->Gamma), end(it->Gamma), static_cast<double>(INFINITY));
   }

   // coarsen the graph
   vector<unsigned> C = coarsen_AMG(it, params);

   unsigned M_next = static_cast<unsigned>(C.size());

   cout << "l = " << l << "  M_next = " << M_next << endl;

   // if no further coarsening is obtained, every node represents a salient 
   // segment.
   if ( M == M_next ) {
      return eye_crs<double>(M,M);
   }

   // We're going down another level, so "make room" in levels
   unsigned l_next = l + 1;

   image_level next_level;
   levels.push_back(next_level);

   // build interpolation matrix P
   matrix_crs<double> P = build_interp(it->A, C, M, M_next);
   matrix_crs<double> P_trans = transpose(P);

   // build column-scaled interpolation matrix P
   matrix_crs<double> Ptil = build_scaled_interp(P);
   matrix_crs<double> Ptil_trans = transpose(Ptil);

   // coarse-level intensity vector
   next(it)->I = Ptil_trans*(it->I);

   // coarse-level variance
   valarray<double> I2 = (it->I)*(it->I);
   valarray<double> nextI2 = (next(it)->I)*(next(it)->I);
   next(it)->S = coarse_variance(Ptil_trans*(it->S), Ptil_trans*I2 - nextI2);


   // coarse-level coupling matrix
   next(it)->A = (it->A)*P;
   next(it)->A = P_trans*(next(it)->A);
   //next(it)->A.clean();
   rescale_coarse_coupling(next(it)->A, it, l_next, params);

   // coarse-level area matrices
   // "W = next(it)->A;"
   next(it)->V = P_trans*(it->V)*P;
   //next(it)->V.clean();

   // coarse-level boundary length matrices
   // XXX: I think we only need the diagonal of L,G
   build_L(next(it), params);
   build_G(next(it), params);
   //next(it)->L.clean();
   //next(it)->G.clean();

   // Build new saliency vector
   valarray<double> diagL = diag(next(it)->L);
   valarray<double> diagV = diag(next(it)->V);
   valarray<double> diagG = diag(next(it)->G);
   valarray<double> diagW = diag(next(it)->A);

   next(it)->Gamma.resize(M_next,0.);
   double tmp = 0.;
   for (unsigned i = 0; i < M_next; ++i) {
      tmp = (diagL[i]/diagG[i]) / (diagW[i]/diagV[i]);
      if ( it->Gamma[C[i]] == 0 ) continue;
      else if ( tmp > params.gamma ) next(it)->Gamma[i] = tmp;
      else continue;
   }

   // Recurse
   matrix_crs<double> U = image_vcycle(l_next, M_next, levels, next(it),
         params);

   // Interpolate U
   U = P*U;

   // Sharpen overlapping segments
   for (unsigned ind = 0; ind < U.val.size(); ++ind) {
      if ( U.val[ind] < params.d1 ) U.val[ind] = 0.;
      else if ( U.val[ind] > 1.-params.d1 ) U.val[ind] = 1.;
      else continue;
   }

   // We probably introduced some zeros into U via P*U and also the sharpening.
   U.clean();
   
   return U;
// }}}
}


// Graph Coarsening
///////////////////
vector<unsigned> coarsen_AMG(const list<image_level>::iterator it, 
      const seg_params& params) {
// {{{
// Coarsen the graph using something similar to the "standard" AMG coarsener
//

   unsigned M = it->A.m;

   vector<unsigned> row_ptr(M+1,0), col_ind;
   vector<double> val;

   vector<unsigned> lambda(M);

   // Build A_bar
   // and count the number of nonzeros in each column
   for (unsigned i = 0; i < M; ++i) {
      row_ptr[i+1] = row_ptr[i];
     
      // compute row sum (excluding diagonal element, if present)
      double row_sum = 0.;
      for (unsigned jp = it->A.row_ptr[i]; jp < it->A.row_ptr[i+1]; ++jp) {
         if ( it->A.col_ind[jp] != i ) {
            row_sum += it->A.val[jp];
         }
      }

      // Use max off-diagonal element, instead of row_sum 
      //double max_elem = -static_cast<double>(INFINITY);
      //for (unsigned jp = it->A.row_ptr[i]; jp < it->A.row_ptr[i+1]; ++jp) {
      //   if ( it->A.col_ind[jp] != i && it->A.val[jp] > max_elem ) {
      //      max_elem = it->A.val[jp];
      //   }
      //}

      // XXX: Testing coarsen_AMG via test.cpp:fig8_4
      //double min_elem = *min_element(it->A.row_ptr.begin() + it->A.row_ptr[i],
      //      it->A.row_ptr.begin() + it->A.row_ptr[i+1]);

      // now add elements to A_bar
      for (unsigned jp = it->A.row_ptr[i]; jp < it->A.row_ptr[i+1]; ++jp) {
         unsigned j = it->A.col_ind[jp];
         double Aij = it->A.val[jp];

         if ( j != i && Aij >= params.theta*row_sum ) {
         //if ( j != i && Aij >= params.theta*max_elem ) {
         //if ( j != i && -Aij >= -params.theta*min_elem ) { // XXX: Testing coarsen_AMG via test.cpp:fig8_4
            ++row_ptr[i+1];
            col_ind.push_back(j);
            val.push_back(Aij);

            // count the number of nonzeros in each column
            lambda[j] += 1;
         }
      }
   }

   matrix_crs<double> A_bar(row_ptr, col_ind, val, M, M, 1);
   matrix_crs<double> A_bar_trans = transpose(A_bar);
   row_ptr.resize(0);
   col_ind.resize(0);
   val.resize(0);

   // start assigning C/F nodes
   // T[i] = 0  <- unassigned
   // T[i] = 1  <- C-point
   // T[i] = 2  <- F-point
   vector<unsigned> T(M,0);
   vector<unsigned> Teqz; Teqz.reserve(M); // vector of indexes j where T[j] != 0.  
   // This is used for a performance optimization.  strongly_influence* can 
   // probably be optimized to use Teqz instead of just T, but that's a TODO
   //
   // In an earlier version of this optimization, it was actually slower to use
   // a std::list than a vector.  Is this due to storage of all those 64-bit pointers?

   // if salient, designate as C node
   for (unsigned i = 0; i < M; ++i) {
      if ( it->Gamma[i] < params.gamma ) {
         T[i] = 1;
         lambda[i] = 0;
      }
      else {
         Teqz.push_back(i);
      }
   }

   //// while any of the T[i] are zero
   //// equivalently, while not all T[i] > 0
   //vector<unsigned>::iterator Teqz_it;
   //auto Teqz_it_end = Teqz.end();

   //while ( Teqz.size() > 0 ) {

   //   unsigned ml = 0;
   //   Teqz_it = Teqz.begin(); // in case lambda is a vector of zeros
   //   Teqz_it_end = Teqz.end();
   //   for (auto it = Teqz.begin(); it != Teqz_it_end; ++it) { // TODO This is the bottleneck in the code.  I think this makes the code O(n^2) in num pixels
   //      if ( lambda[*it] > ml ) {
   //         Teqz_it = it;
   //         ml = lambda[*it];
   //      }
   //   }

   //   unsigned jm = *Teqz_it;
   //   T[jm] = 1;
   //   Teqz.erase(Teqz_it);
   //   lambda[jm] = 0;

   //   vector<unsigned> K = strongly_influenced_by_j_trans(A_bar_trans, T, jm);

   //   vector<unsigned> H;
   //   Teqz_it = Teqz.begin();
   //   for (auto k_it = K.begin(); k_it != K.end(); ++k_it) {
   //      
   //      T[*k_it] = 2; // nodes in K become F-points
   //      Teqz_it = lower_bound(Teqz_it, Teqz.end(), *k_it); // this will always find the value
   //      Teqz_it = Teqz.erase(Teqz_it);

   //      lambda[*k_it] = 0;

   //      H = strongly_influence_k(A_bar, T, *k_it);
   //      for (auto h_it = H.begin(); h_it != H.end(); ++h_it) {
   //         ++lambda[*h_it];
   //      }
   //   }
   //}


   // Here's some stuff I SHAMELESSLY stole from pyamg's ruge_stuben.h
   // XXX This implements the standard RS AMG coarsening.  This is slightly different than the coarsening presented by Inglis et al.  However, when I ran against the 60x60 checkered disk image with their parameters (row-sum strengths), I got the results they found.  So, maybe this algorithm will work.
  
   //for each value of lambda, create an interval of nodes with that value
   // ptr - is the first index of the interval
   // count - is the number of indices in that interval
   // index to node - the node located at a given index
   // node to index - the index of a given node
   vector<unsigned> interval_ptr(M+1,0);
   vector<unsigned> interval_count(M+1,0);
   vector<unsigned> index_to_node(M);
   vector<unsigned> node_to_index(M);

   for (unsigned i = 0; i < M; ++i) {
      interval_count[lambda[i]]++;
   }
   for (unsigned i = 0, cumsum = 0; i < M; ++i) {
      interval_ptr[i] = cumsum;
      cumsum += interval_count[i];
      interval_count[i] = 0;
   }
   for (unsigned i = 0; i < M; ++i) {
      unsigned lambda_i = lambda[i];
      unsigned index    = interval_ptr[lambda_i] + interval_count[lambda_i];
      index_to_node[index] = i;
      node_to_index[i]     = index;
      ++interval_count[lambda_i];
   }

   //Now assign nodes to C/F in T, in descending order of lambda
   unsigned jj_end = 0;
   unsigned kk_end = 0;
   for (int top_index = M - 1; top_index != -1; --top_index) {
      unsigned i        = index_to_node[top_index];
      unsigned lambda_i = lambda[i];


      //remove i from its interval
      --interval_count[lambda_i];

      if ( T[i] == 2 )
      {
         continue;
      } 
      else 
      {

         T[i] = 1;

         //For each j in S^T_i /\ U
         jj_end = A_bar_trans.row_ptr[i+1];
         for(unsigned jj = A_bar_trans.row_ptr[i]; jj < jj_end; ++jj){
            unsigned j = A_bar_trans.col_ind[jj];

            if ( T[j] == 0 ) {
               T[j] = 2;

               //For each k in S_j /\ U
               kk_end = A_bar.row_ptr[j+1];
               for (unsigned kk = A_bar.row_ptr[j]; kk < kk_end; ++kk) {
                  unsigned k = A_bar.col_ind[kk];

                  if ( T[k] == 0 ) {         
                     //move k to the end of its current interval
                     if ( lambda[k] >= M - 1 ) continue;

                     unsigned lambda_k = lambda[k];
                     unsigned old_pos  = node_to_index[k];
                     unsigned new_pos  = interval_ptr[lambda_k] + interval_count[lambda_k] - 1;

                     node_to_index[index_to_node[old_pos]] = new_pos;
                     node_to_index[index_to_node[new_pos]] = old_pos;
                     swap(index_to_node[old_pos], index_to_node[new_pos]);

                     //update intervals
                     interval_count[lambda_k]   -= 1;
                     interval_count[lambda_k+1] += 1; //invalid write!
                     interval_ptr[lambda_k+1]    = new_pos;

                     //increment lambda_k
                     ++lambda[k];
                  }
               }
            }
         }

         //For each j in S_i /\ U
         jj_end = A_bar.row_ptr[i+1];
         for(unsigned jj = A_bar.row_ptr[i]; jj < jj_end; ++jj){
            unsigned j = A_bar.col_ind[jj];
            if ( T[j] == 0 ) {            //decrement lambda for node j
               if ( lambda[j] == 0 ) continue;

               //assert(lambda[j] > 0);//this would cause problems!

               //move j to the beginning of its current interval
               unsigned lambda_j = lambda[j];
               unsigned old_pos  = node_to_index[j];
               unsigned new_pos  = interval_ptr[lambda_j]; 

               node_to_index[index_to_node[old_pos]] = new_pos;
               node_to_index[index_to_node[new_pos]] = old_pos;
               std::swap(index_to_node[old_pos],index_to_node[new_pos]);

               //update intervals
               interval_count[lambda_j]   -= 1;
               interval_count[lambda_j-1] += 1;
               interval_ptr[lambda_j]     += 1;
               interval_ptr[lambda_j-1]    = interval_ptr[lambda_j] - interval_count[lambda_j-1];

               //decrement lambda_j
               --lambda[j];
            }
         }
      }
   }     


   // Assign C-points
   vector<unsigned> C;
   for (unsigned i = 0; i < T.size(); ++i) {
      if ( T[i] == 1 ) C.push_back(i); // XXX: valgrind shows this (maybe!) as uninitialized and gives an error message; not sure I belive that..
   }

   if ( C.size() == 0 ) {
      cerr << "error: image_seg.cpp:coarsen_AMG: Vector of C-points is empty!"
           << endl;
      exit(-1);
   }

   return C;
// }}}
}

vector<unsigned> strongly_influenced_by_j(const matrix_crs<double>& A_bar,
      const vector<unsigned>& T, const unsigned j) {
// {{{
// K is the set of nodes k such that T[k] == 0 and A_bar[k,j] > 0
// Such nodes are strongly influenced by node j
//
// TODO: This is slow in CRS.  I don't think A_bar is symmetric

   vector<unsigned> K;

   // Look down columns (which is slow in CRS!)
   // We check if T_k == 0 first, because that is much cheaper than 
   // checking if A_bar[k,j] > 0 (think "constraint propagation")
   for (unsigned k = 0; k < A_bar.m; ++k) {

      if ( T[k] != 0 ) continue;

      // we now know T[k] == 0, so check if A_bar[k,j] > 0
      
      for (unsigned cp = A_bar.row_ptr[k]; cp < A_bar.row_ptr[k+1]; ++cp) {
         unsigned c = A_bar.col_ind[cp];

         if ( c < j ) continue;
         else if ( c == j ) { // value exists at [k,j] in A_bar
            if ( A_bar.val[cp] > 0 ) {
               //cout << "[k,j] = [" << k << "," << j << "]" << endl;
               K.push_back(k);
               break;
            }
         }
         else { // c > j, so value doesn't exist in matrix
            break;
         }
      }

      // This is actually slower than the loop above
      //auto it_cp = lower_bound(A_bar.col_ind.begin() + A_bar.row_ptr[k], A_bar.col_ind.begin() + A_bar.row_ptr[k+1], j);
      //unsigned c = *it_cp;
      //if ( it_cp != A_bar.col_ind.begin() + A_bar.row_ptr[k+1] && c == j && A_bar.val[distance(A_bar.col_ind.begin(), it_cp)] > 0. ) { // A_bar[k,j] > 0
      ////if ( it_cp != A_bar.col_ind.begin() + A_bar.row_ptr[k+1] && 
      ////      c == j && A_bar.val[distance(A_bar.col_ind.begin(), it_cp)] != 0. ) { // XXX: Testing coarsen_AMG via test.cpp:fig8_4

      //   K.push_back(k);
      //}
   }

   //cout << "K = " << endl;
   //print_vector(K);
 
   return K;
// }}}
}

vector<unsigned> strongly_influenced_by_j_trans(const matrix_crs<double>& A_bar_trans,
      const vector<unsigned>& T, const unsigned j) {
// {{{
// K is the set of nodes k such that T[k] == 0 and A_bar_trans[j,k] > 0
// Such nodes are strongly influenced by node j
//
// This should give the same result as strongly_influenced_by_j, but should 
// be considerably more efficient

   vector<unsigned> K;

   // Look down rows of transpose(A_bar)
   for (unsigned kp = A_bar_trans.row_ptr[j]; kp < A_bar_trans.row_ptr[j+1]; ++kp) {
      unsigned k = A_bar_trans.col_ind[kp];

      if ( T[k] == 0 && A_bar_trans.val[kp] > 0 ) {
         K.push_back(k);
      }
   }
   
   return K;
// }}}
}

vector<unsigned> strongly_influence_k(const matrix_crs<double>& A_bar,
      const vector<unsigned>& T, const unsigned k) {
// {{{
// H is the set of nodes h such that T[h] == 0 and A_bar[k,h] > 0
// Such nodes strongly influence node k (k has become an F-point)

   vector<unsigned> H;

   // We look down row k of A_bar, which is cheap in CRS.

   for (unsigned hp = A_bar.row_ptr[k]; hp < A_bar.row_ptr[k+1]; ++hp) {
      unsigned h = A_bar.col_ind[hp];

      if ( T[h] == 0 && A_bar.val[hp] > 0. ) {
      //if ( T[h] == 0 && A_bar.val[hp] != 0. ) { // XXX: Testing coarsen_AMG via test.cpp:fig8_4
         //cout << "[k,h] = [" << k << "," << h << "]" << endl;
         H.push_back(h);
      }
   }

   return H;
// }}}
}

// }}}


///////////////////////
// Image Seg Drivers //
///////////////////////
// {{{

// main driver for image segmentation
matrix_crs<double> image_seg(const cv::Mat& img, seg_params& params) {
// {{{
  
   // check if grayscale (single channel)
   if ( img.channels() != 1 ) {
      cerr << "error: image_seg.cpp:image_seg: Image segmentation expects "
           << "an input matrix that has only one channel." << endl;
      exit(-1);
   }

   valarray<double> I = image_to_intensity(img, params);

   // First Level stuff
   list<image_level> levels(1);
   build_first_level(I, params, levels.begin());

   // Do V-Cycle
   unsigned l = 1;
   unsigned M = static_cast<unsigned>(I.size());
   matrix_crs<double> U = image_vcycle(l, M, levels, levels.begin(), params);

   // Assign pixels uniquely
   assign_uniquely(U);

   return U;
// }}}
}

// Overloaded drivers
/////////////////////
// this uses default segmentation parameters.  You probably shouldn't use this
matrix_crs<double> image_seg(const cv::Mat& img) {
// {{{
   // (not necessarily good...) default parameters 
   // eqn 21 of Inglis et al. 2010
   seg_params params(100.,100.,100.,0.1,0.1,0.15,5,1);

   return image_seg(img, params);
// }}}
}

// Helpers
//////////
void set_params(seg_params& params, const double alpha, 
      const double alpha_til, const double beta, const double theta,
      const double gamma, const double d1, const unsigned sigma,
      const unsigned rho) {
// {{{
   params.alpha     = alpha;
   params.alpha_til = alpha_til;
   params.beta      = beta;
   params.theta     = theta;
   params.gamma     = gamma;
   params.d1        = d1;
   params.sigma     = sigma;
   params.rho       = rho;
// }}}
}

void assign_uniquely(matrix_crs<double>& U) {
// {{{

   // Build up new CRS matrix (probably faster than erasing in place)
   vector<unsigned> row_ptr(U.m+1,0);
   vector<unsigned> col_ind; col_ind.reserve(U.m);
   vector<double> val; val.reserve(U.m);

   vector<double>::iterator Uval_it = U.val.begin();
   vector<double>::iterator Uval_it_end;
   for (unsigned i = 0; i < U.m; ++i) {
      row_ptr[i] = i;

      // Find max element in this row
      Uval_it = U.val.begin() + U.row_ptr[i];
      Uval_it_end = U.val.begin() + U.row_ptr[i+1];
      Uval_it = max_element(Uval_it, Uval_it_end);

      val.push_back(1.);
      col_ind.push_back(U.col_ind[Uval_it - U.val.begin()]);
   }
   row_ptr[U.m] = U.m;

   // assign new CRS components to U
   U.row_ptr = row_ptr;
   U.col_ind = col_ind;
   U.val = val;

// }}}
}

vector<valarray<double>> mat_to_vecs(const matrix_crs<double>& U) {
// {{{
// Convert a matrix U to a vector of valarrays representing the columns
// This is used to convert the segmentation Boolean matrix U to (hopefully
// not too many!) intensity vectors.

   vector<valarray<double>> intensities(U.n);

   // pre-allocate memory
   for (unsigned j = 0; j < U.n; ++j) {
      intensities[j].resize(U.m,0.); // init to zero
   }

   // copy over intensity data
   for (unsigned i = 0; i < U.m; ++i) {
      for (unsigned jp = U.row_ptr[i]; jp < U.row_ptr[i+1]; ++jp) {
         intensities[U.col_ind[jp]][i] = U.val[jp];
      }
   }

   return intensities;
// }}}
}

// }}}


///////////////////
// OpenCV things //
///////////////////
// {{{

cv::Mat load_image(const string& img_filename, const unsigned mode) {
// {{{
// Load image from file specified by filename into a CV Mat with uchar entries
// The image is blended to a grayscale image and sacled so the pixels
// have values between 0 and 1
// Matrices must be size n x n for the image segmentation method
//
// mode = 0 - default; grayscale
// mode = 1 - 3 color channels
//

   cv::Mat img;
   img = cv::imread(img_filename, mode); // uchar [0,255] vals

   // Check to see if we read the image properly
   if ( img.data == NULL ) {
      cerr << "error: image_seg.cpp:load_image: OpenCV imread loaded empty "
           << "matrix." << endl
           << "       filename = " << img_filename << endl;
      exit(-1);
   }

   // Check if the image is square
   if ( img.rows != img.cols && img.rows > 0 ) {
      cerr << "error: image_seg.cpp:load_image: the image \"" << img_filename 
           << "\" is not square." << endl;
      exit(-1);
   }

   return img;
// }}}
}


valarray<double> image_to_intensity(cv::Mat img, seg_params& params) {
// {{{
// Compute the intensity vector from a matrix read in via `load_image`
// The intensity vector has values scaled from 0 to 1 (0->0, 255->1)
// The entries of the vector are the entries of the image matrix sorted
// row-wise
   valarray<double> I(0., img.rows*img.cols);
  
   // Get the max pixel value of the image so we can rescale [0,1]
   double minval, maxval;
   cv::minMaxLoc(img, &minval, &maxval);

   // Set the size of the image in the parameters class
   // We'll use it later to build A on the first level
   params.n = static_cast<unsigned>(img.rows);

   for (int i=0; i < img.rows; ++i) {

      const uchar* img_i = img.ptr<uchar>(i); // pointer to row i
      for (int j=0; j < img.cols; ++j) {
         // Rachel pointed out that this is wrong
         // scale to [0,1] from range [0, 255] (uchar)
         //I[img.rows*i + j] = static_cast<double>( img_i[j] )/ UCHAR_MAX;

         // This should be the right scaling
         I[img.rows*i + j] = (static_cast<double>(img_i[j]) - minval)
            / (maxval - minval);
      }
   }

   //cout << "max = " << *max_element(begin(I), end(I)) << endl;
   //cout << "min = " << *min_element(begin(I), end(I)) << endl;
   return I;
// }}}
}


cv::Mat intensity_to_image1(const valarray<double>& I, const unsigned n) {
// {{{
// Convert an intensity vector (scaled 0 to 1) to a grayscale cv::Mat "image"
   cv::Mat img(n, n, CV_8UC1);

   if ( I.size() != n*n ) {
      cerr << "error: image_seg.cpp:intensity_to_image: intensity vector "
           << "and image size do not match" << endl;
      exit(-1);
   }

   uchar* img_i = img.ptr<uchar>(0); // pointer to row i, but since img
                                     // is continuous, the pointer extends
                                     // to the whole array
   for (unsigned ind = 0; ind < I.size(); ++ind) {
      img_i[ind] = cv::saturate_cast<uchar>(UCHAR_MAX*I[ind]);
   }

   return img;
// }}}
}


void write_seg_images(const cv::Mat& orig_img, const matrix_crs<double>& U,
      const string& filename_base, const unsigned mode) {
// {{{
// Write the segments specified by the matrix U to files with base filename
// given by filename_base.  The files will be something like
//    filename_base_seg_01.jpg
//    filename_base_seg_02.jpg
//    filename_base_seg_03.jpg
//    ...
// 
// Each image should be orig_img.rows x orig_img.cols ( square )
//
// mode = 0 - Just write images of the segments; don't blend with the
//            original image
//
// mode = 1 - Blend each segment with the original image
//
// TODO
// mode = 2 - Add (saturate) segment on top of original
//

   vector<valarray<double>> intensities = mat_to_vecs(U);

   if ( mode == 0 ) { // segment only
      
      for (unsigned i = 0; i < intensities.size(); ++i) {
         // Make filename
         string filename = filename_base;
         filename += "_seg_";
         // add leading zeros
         if ( i > 99 ) filename += "0";
         else if ( i > 9 ) filename += "00";
         else filename += "000";
         filename += to_string(i);
         filename += ".png";

         try {
            cout << "writing segment image file: " << filename << endl;
            cv::imwrite(filename, 
                  intensity_to_image1(intensities[i], orig_img.rows));
         }
         catch ( runtime_error& ex ) {
            cerr << "error: image_seg.cpp:write_seg_images: Something went "
                 << "wrong while writing images." << endl 
                 << "exception: " << ex.what();
            exit(-1);
         }
      }
   }

   else if ( mode == 1 ) { // original + segment blend

      // RGB images
      cv::Mat orig_rgb(orig_img.size(), CV_8UC3);
      cv::cvtColor(orig_img, orig_rgb, CV_GRAY2RGB);

      cv::Mat seg_rgb(orig_img.size(), CV_8UC3);
      cv::Mat blend_rgb(orig_img.size(), CV_8UC3);
      
      double alpha = 0.5; // proportion of segment vs. original
      
      uchar color[] = {128, 128, 0}; // BGR - teal
      
      for (unsigned i = 0; i < intensities.size(); ++i) {
         // Make filename
         string filename = filename_base;
         filename += "_seg_blend_";
         // add leading zeros
         if ( i > 99 ) filename += "0";
         else if ( i > 9 ) filename += "00";
         else filename += "000";
         filename += to_string(i);
         filename += ".png";

         // make the segment some color
         cv::cvtColor(intensity_to_image1(intensities[i], orig_img.rows),
               seg_rgb, CV_GRAY2RGB);

         cv::Mat channels[3];
         split(seg_rgb, channels);

         for (unsigned c = 0; c < 3; ++c) {

            unsigned ind_end = channels[c].rows*channels[c].cols;
            uchar* pix_ptr = channels[c].ptr<uchar>(0); // pointer to 
                                                        // beginning of img
            for ( unsigned ind = 0; ind < ind_end; ++ind ) {
               if ( pix_ptr[ind] > 0 ) {   // if the pixel is non-zero,
                  pix_ptr[ind] = color[c]; // set the color
               }
            }
         }

         merge(channels, 3, seg_rgb);

         // blend the segment and original
         // blend_rgb = alpha*seg_rgb + (1-alpha)*orig_rgb + 0.0
         cv::addWeighted(seg_rgb, alpha, orig_rgb, 1.-alpha, 0.0, blend_rgb);

         try {
            cout << "writing original+segment blend image file: "
                 << filename << endl;
            cv::imwrite(filename, blend_rgb);
         }
         catch ( runtime_error& ex ) {
            cerr << "error: image_seg.cpp:write_seg_images: Something went "
                 << "wrong while writing images." << endl 
                 << "exception: " << ex.what();
            exit(-1);
         }
      }
   }

   else {
      cerr << "error: image_seg.cpp:write_seg_images: Bad mode give: "
           << "mode = " << mode << endl;
      exit(-1);
   }

// }}}
}

// }}}


///////////////////////
// Algorithim Timing //
///////////////////////
// {{{

void runtime_scaling(const cv::Mat& img, seg_params& params, const string& out_filename) {

   cv::Mat img_resize;
   unsigned num_imgs = 20, num_samples = 5;
   unsigned max_size = img.rows, min_size = 10, pix_lin_step;
   vector<unsigned> img_sizes(num_imgs,0);

   // sample images sizes to pixel numbers are linear
   pix_lin_step = static_cast<double>(max_size*max_size-min_size*min_size)/static_cast<double>(num_imgs-1);
   img_sizes[0] = min_size;
   for (unsigned i = 1; i < num_imgs; ++i) {
      img_sizes[i] = static_cast<unsigned>(sqrt(img_sizes[i-1]*img_sizes[i-1] + pix_lin_step));
   }
   if ( img_sizes[num_imgs-1] != max_size ) {
      img_sizes.push_back(max_size);
   }
  
   // open data output file
   ofstream outfile;
   outfile.open(out_filename, ios::out | ios::trunc);
   if ( !outfile.is_open() ) {
      cerr << "error: image_seg.cpp:runtime_scaling: I had an issue opening the file \""
           << out_filename << "\"" << endl;
      exit(-1);
   }

   // rescale images and run num_samples runs of image segmentation
   matrix_crs<double> U;
   vector<unsigned> runtimes_ms(img_sizes.size(), 0);

   for (unsigned i = 0; i < img_sizes.size(); ++i) {
      cv::resize(img, img_resize, cv::Size(img_sizes[i], img_sizes[i]));
      cout << "i = " << i << "  rows = " << img_resize.rows << endl;

      // run num_samples image segmentations
      auto t_start = chrono::high_resolution_clock::now();
      for (unsigned samp = 0; samp < num_samples; ++samp) {
         U = image_seg(img_resize, params);
      }
      auto t_end = chrono::high_resolution_clock::now();

      // average time
      runtimes_ms[i] = static_cast<unsigned>(static_cast<double>(
               chrono::duration_cast<chrono::milliseconds>(t_end-t_start).count())
         / static_cast<double>(num_samples));
   
   }

   // write number of pixels and average runtimes to outfile for plotting with gnuplot
   for (unsigned i = 0; i < img_sizes.size(); ++i) {
      outfile << img_sizes[i]*img_sizes[i] << "  " << runtimes_ms[i] << endl;
   }

   outfile.close();
}

// }}}


int main(void) {

   cv::Mat img;
   seg_params params;
   matrix_crs<double> U;

   // Easy test cases
   //////////////////
   //img = load_image("../test_imgs/square.png");
   //set_params(params, 10., 5., 10., 0.1, 0.15, 0.15, 2, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/square", 1);

   //img = load_image("../test_imgs/square_inv.png");
   //set_params(params, 10., 5., 10., 0.1, 0.15, 0.15, 2, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/square_inv", 1);

   //img = load_image("../test_imgs/arrow_5.png");
   //set_params(params, 10., 10., 100., 0.1, 0.1, 0.15, 1, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/arrow_5", 1);
 
   //img = load_image("../test_imgs/squares.png");
   //set_params(params, 10., 5., 100., 0.1, 0.1, 0.15, 1, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/squares", 1);

   //img = load_image("../test_imgs/E_25.png");
   //set_params(params, 10., 5., 100., 0.1, 0.1, 0.15, 1, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/E", 1);

   //img = load_image("../test_imgs/arrow_25.png");
   //set_params(params, 10, 5., 100., 0.1, 0.1, 0.15, 3, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/arrow_25", 1);

   // Checker Disk
   ///////////////
   // Inglis et al.'s parameters
   //set_params(params, 10., 10., 10., 0.1, 0.1, 0.15, 5, 1);
   ////set_params(params, 20., 10., 10., 0.1, 0.1, 0.15, 5, 1); // a bit more contrast
   ////set_params(params, 10., 10., 10., 0.1, 0.1, 0.15, 5, 1); // change theta to use max-row or row-sum
   //img = load_image("../test_imgs/checker_disk_60.png");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/checker_disk_60", 1);
 
   //img = load_image("../test_imgs/checker_disk_120.png");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/checker_disk_120", 1);

   // Checker Disk Larger
   //////////////////////
   // Inglis et al.'s parameters
   //set_params(params, 10., 10., 10., 0.1, 0.1, 0.15, 5, 1);
   //set_params(params, 12., 10., 10., 0.1, 0.1, 0.15, 5, 1); // a bit more contrast
   //img = load_image("../test_imgs/checker_disk_larger_120.png");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/checker_disk_larger_120", 1);
 
   // Blob
   /////////
   //img = load_image("../test_imgs/blob_64.png");
   //set_params(params, 10., 1., 100.*0., 0.05, 0.10, 0.15, 5, 2);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/blob_64", 1);
 
   //img = load_image("../test_imgs/blob_128.png");
   //set_params(params, 10., 3., 10., 0.10, 0.10, 0.15, 5, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/blob_128", 1);
 
   // Spiral
   /////////
   //img = load_image("../test_imgs/spiral_64.png");
   //set_params(params, 10., 1., 10., 0.05, 0.10, 0.15, 5, 2);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/spiral_64", 1);
 
   //img = load_image("../test_imgs/spiral_256.png");
   ////set_params(params, 10., 1., 10., 0.05, 0.10, 0.15, 5, 2); // Inglis coarsening
   //set_params(params, 10., 1., 10., 0.025, 0.10, 0.15, 5, 2); // RS coarsening 
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/spiral_256", 1);
  
   // Peppers
   //////////
   //set_params(params, 50., 4., 10., 0.10, 0.15, 0.15, 5, 1);
   ////set_params(params, 10., 10., 10., 0.1, 0.1, 0.15, 5, 1);
   //img = load_image("../test_imgs/peppers_25.jpg");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/peppers_25", 1);

   //img = load_image("../test_imgs/peppers_50.jpg");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/peppers_50", 1);

   //img = load_image("../test_imgs/peppers_100.jpg");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/peppers_100", 1);

   //img = load_image("../test_imgs/peppers.jpg");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/peppers", 1);


   // Nazca Lines
   // This is probably a _very_ tough image to segment
   //img = load_image("../test_imgs/nazca3_200.jpg");
   //set_params(params, 30., 4., 1., 0.2, 0.15, 0.15, 5, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/nazca3_500", 1);
 
   //img = load_image("../test_imgs/nazca3_500.jpg");
   //set_params(params, 10., 1., 10., 0.05, 0.10, 0.15, 5, 2);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/nazca3_500", 1);

   // CK
   //img = load_image("../test_imgs/ck.png");
   //set_params(params, 90., 5., 5., 0.10, 0.15, 0.15, 5, 1);
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/ck", 1);

   // C. Elegans front
   // I want to see if we can segment the various parts of the worm.  I'd be surprised if this works.
   //set_params(params, 130., 0.5, 100., 0.10, 0.15, 0.15, 6, 2);
   //img = load_image("../test_imgs/c_elegans_front_256.png");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/c_elegans_front_256", 1);

   //img = load_image("../test_imgs/c_elegans_front_512.png");
   //U = image_seg(img, params);
   //write_seg_images(img, U, "gen_imgs/c_elegans_front_512", 1);
 
   // Two C. Elegans 
   // Can we tell the two worms apart?
   //set_params(params, 2., 23., 2., 0.05, 0.05, 0.10, 7, 1);
   set_params(params, 2., 25., 2., 0.05, 0.05, 0.10, 7, 1);
   img = load_image("../test_imgs/two_c_elegans.jpg");
   U = image_seg(img, params);
   write_seg_images(img, U, "gen_imgs/two_c_elegans", 1);


   
   // Algorithm timing stuff
   // We're currently sitting on a nice O(n^2) in the number of pixels.  This is almost entirely due to something slow thhat I'm doing in build_interp and coarsenAMG.  These two account for almost all of the cost of running the code!
   //set_params(params, 50., 4., 10., 0.10, 0.15, 0.15, 5, 1); // Inglis coarsening
   //set_params(params, 50., 4., 10., 0.10, 0.15, 0.15, 5, 1); // RS coarsening
   //img = load_image("../test_imgs/peppers.jpg");
   
   //set_params(params, 10., 1., 10., 0.05, 0.10, 0.15, 5, 2); // Inglis coarsening
   set_params(params, 10., 1., 10., 0.025, 0.10, 0.15, 5, 2); // RS coarsening 
   img = load_image("../test_imgs/spiral_512.png");
   
   //runtime_scaling(img, params, "gen_imgs/runtime_scaling.out");    


   return 0;
}
