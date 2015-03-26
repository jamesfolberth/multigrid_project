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
   it->S = zeros<double>(I.size(), I.size());
 
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
                    << " top row, left side" << endl;
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
                    << " top row, right side" << endl;
            }
 
            ++row_ptr[ind+1];
            col_ind.push_back(ind-1);
            val.push_back(exp(-a*abs(I[ind]-I[ind-1])));

            ++row_ptr[ind+1];
            col_ind.push_back(ind+n);
            val.push_back(exp(-a*abs(I[ind]-I[ind+n])));
         }

         else { // in the middle
            ++row_ptr[ind+1];
            if ( _DEBUG_ >= 2) {
               cout << "debug 2: ind = " << ind 
                    << " top row, middle side" << endl;
            }
 
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
                    << " middle row, left side" << endl;
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
                    << " middle row, right side" << endl;
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
                    << " middle row, middle side" << endl;
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
                    << " bottom row, left side" << endl;
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
                    << " bottom row, right side" << endl;
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
                    << " bottom row, middle side" << endl;
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
// in julia.
//
//TODO: this is a bottleneck in the code
//      look into std::lower_bound
//      http://www.cplusplus.com/reference/algorithm/lower_bound/

   it->L = -1.*it->A;

   for (unsigned i = 0; i < it->L.m; ++i) {
      
      // compute full row sum
      double sum = 0.;
      for (unsigned jp = it->L.row_ptr[i]; jp < it->L.row_ptr[i+1]; ++jp) {
         sum += it->L.val[jp];
      }

      // Insert into matrix
      for (unsigned jp = it->L.row_ptr[i]; jp < it->L.row_ptr[i+1]; ++jp) {
         unsigned j = it->L.col_ind[jp];

         if ( j < i && i != it->L.m-1 ) {
            continue;
         }

         else if ( j == i ) {
            // this spot already exists in the matrix,
            // so we don't need to insert

            it->L.val[jp] -= sum;
            break;
         }

         else { // j > i or i == m-1 (last row)
            // the spot j == i doesn't already exist in the matrix, 
            // so we need to insert a new element and adjust row_ptr
            
            if ( i != it->L.m-1 ) {
               it->L.col_ind.insert(it->L.col_ind.begin() + jp, i);
               it->L.val.insert(it->L.val.begin() + jp, -sum);
            }
            else {
               it->L.col_ind.insert(it->L.col_ind.end(), i);
               it->L.val.insert(it->L.val.end(), -sum);
            }

            // adjust row pointers
            for ( unsigned row = i+1; row <= it->L.m; ++row) {
               it->L.row_ptr[row] += 1;
            }

            break;
         }
      }
   }
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

   for (unsigned i = 0; i < it->G.m; ++i) {
      
      // compute full row sum
      double sum = 0.;
      for (unsigned jp = it->G.row_ptr[i]; jp < it->G.row_ptr[i+1]; ++jp) {
         sum += it->G.val[jp];
      }

      // Insert into matrix
      for (unsigned jp = it->G.row_ptr[i]; jp < it->G.row_ptr[i+1]; ++jp) {
         unsigned j = it->G.col_ind[jp];

         if ( j < i && i != it->G.m-1 ) {
            continue;
         }

         else if ( j == i ) {
            // this spot already exists in the matrix,
            // so we don't need to insert

            it->G.val[jp] -= sum;
            break;
         }

         else { // j > i or i == m-1 (last row)
            // the spot j == i doesn't already exist in the matrix, 
            // so we need to insert a new element and adjust row_ptr
            
            if ( i != it->G.m-1 ) {
               it->G.col_ind.insert(it->G.col_ind.begin() + jp, i);
               it->G.val.insert(it->G.val.begin() + jp, -sum);
            }
            else {
               it->G.col_ind.insert(it->G.col_ind.end(), i);
               it->G.val.insert(it->G.val.end(), -sum);
            }

            // adjust row pointers
            for ( unsigned row = i+1; row <= it->G.m; ++row) {
               it->G.row_ptr[row] += 1;
            }

            break;
         }
      }
   }
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

   vector<unsigned> row_ptr(M+1,0), col_ind;
   vector<double> val;

   unsigned C_ind = 0;
   for (unsigned i = 0; i < M; ++i) {
      row_ptr[i+1] = row_ptr[i];

      if ( i == C[C_ind] ) { // inject C-point up
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
         vector<unsigned> A_inds; // vector of indexes for A.col_ind/A.val
         A_inds.reserve( MIN((A.row_ptr[i+1]-A.row_ptr[i]), C.size()) );
         unsigned jp = A.row_ptr[i];

         vector<unsigned> C_inds; // vector of indexes for C
         C_inds.reserve( MIN((A.row_ptr[i+1]-A.row_ptr[i]), C.size()) );
         unsigned j = 0;
 
         double row_sum = 0.;
         while ( jp < A.row_ptr[i+1] && j < C.size() ) {
            if ( A.col_ind[jp] < C[j] ) ++jp;
            else if ( C[j] < A.col_ind[jp] ) ++j;
            else {
               //cout << A.col_ind[jp] << endl;
               row_sum += A.val[jp];

               A_inds.push_back(jp);
               C_inds.push_back(j);

               ++jp; ++j;
            }
         }

         //cout << "i = " << i << endl;
         //cout << "A row = " << endl;
         //for (auto it = A.col_ind.begin() + A.row_ptr[i];
         //      it != A.col_ind.begin() + A.row_ptr[i+1]; ++it) {
         //   cout << "        " << *it << endl;
         //}

         //cout << "C = " << endl;
         //print_vector(C);

         //cout << "A_inds = " << endl;
         //print_vector(A_inds);
 
         //cout << "C_inds = " << endl;
         //print_vector(C_inds);
 
         // now that we know the row sum and the proper indexes,
         // assign values to row of P
         for (j = 0; j < C_inds.size(); ++j) {
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

   return Ptil;
// }}}
}

matrix_crs<double> coarse_variance(const matrix_crs<double>& Sf, 
      const valarray<double>& Sc) {
// {{{
// Build the coarse-level variance matrix S by concatenating the restricted
// fine-level variance matrix Sf and coarse-level variance vector Sc
  
   // copy fine-level variance matrix
   matrix_crs<double> S = Sf;

   // Append the contents of Sc to the right-side of S
   unsigned num_added = 0;
   for (unsigned i = 0; i < S.m; ++i) {
      // We don't want to append a bunch of zeros (Sc usually has some, 
      // at least on the first level)
      if ( Sc[i] != 0. ) {
         S.col_ind.insert(S.col_ind.begin() + S.row_ptr[i+1], S.n);
         S.val.insert(S.val.begin() + S.row_ptr[i+1], Sc[i]);
         ++num_added;
         
         S.row_ptr[i+1] += num_added;
      }
      else continue;
   }

   // Adjust size of S
   ++S.n;

   return S;
// }}}
}

void rescale_coarse_coupling(matrix_crs<double>& A, 
      const list<image_level>::iterator it, const unsigned l_next,
      const seg_params& params) {
// {{{
// Rescale the input matrix A (= A^[r+1]) according to coarse-level variance 
// and also, if l_next >= rho, the multi-level variance
   
   double a_til = params.alpha_til;

   if ( l_next >= params.rho ) { // rescale using intensity and variance 
   
      double b = params.beta;
   
      // Store the address of S, which we'll use repeatedly later
      // XXX: Is this an aliasing issue?
      const matrix_crs<double> *S = &( next(it)->S );

      for (unsigned i = 0; i < A.m; ++i) {
         
         double Ii = next(it)->I[i];
         unsigned sip = S->row_ptr[i];
         unsigned sip_end = S->row_ptr[i+1];

         for (unsigned jp = A.row_ptr[i]; jp < A.row_ptr[i+1]; ++jp) {
            unsigned j = A.col_ind[jp];

            // rescale using coarse-level intensity
            A.val[jp] *= exp(-a_til*abs(Ii - next(it)->I[j]));

            // rescale using multi-level variance
            // Computing the 2-norms is probably slow
            // We're also not utilizing the symmetry of A; that's a factor of 2
            unsigned sjp = S->row_ptr[j];
            unsigned sjp_end = S->row_ptr[j+1];

            // XXX: this norm computation is not scaled like dnrm2
            double nrm = 0.;
            while ( sip < sip_end && sjp < sjp_end ) {
               if ( S->col_ind[sip] < S->col_ind[sjp] ) {
                  nrm += (S->val[sip]) * (S->val[sip]);
                  ++sip;
               }
               else if ( S->col_ind[sjp] < S->col_ind[sip] ) {
                  nrm += (S->val[sjp]) * (S->val[sjp]);
                  ++sjp;
               }
               else {
                  nrm += (S->val[sip]-S->val[sjp]) * (S->val[sip]-S->val[sjp]);
                  ++sip; ++sjp;
               }
            }
            nrm = sqrt(nrm);
            
            A.val[jp] *= exp(-b*nrm);
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
   next(it)->I = Ptil_trans*it->I;

   // coarse-level variance
   valarray<double> I2 = (it->I)*(it->I);
   valarray<double> nextI2 = (next(it)->I)*(next(it)->I);
   next(it)->S = coarse_variance(Ptil_trans*(it->S), Ptil_trans*I2 - nextI2);

   // coarse-level coupling matrix
   next(it)->A = (it->A)*P;
   next(it)->A = P_trans*(next(it)->A);
   rescale_coarse_coupling(next(it)->A, it, l_next, params);

   // coarse-level area matrices
   // W = next(it)->A;
   next(it)->V = (it->V)*P;
   next(it)->V = P_trans*(next(it)->V);

   // coarse-level boundary length matrices
   // TODO: I think we only need the diagonal of L,G
   build_L(next(it), params);
   build_G(next(it), params);

   // Build new saliency vector
   next(it)->Gamma.resize(M_next,0.);
   valarray<double> diagL = diag(next(it)->L);
   valarray<double> diagV = diag(next(it)->V);
   valarray<double> diagG = diag(next(it)->G);
   valarray<double> diagW = diag(next(it)->A);

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

   // We maybe introduced some zeros into U via P*U and also the sharpening.
   U.clean();

   return U;
// }}}
}


// Graph Coarsening
///////////////////
// XXX: I don't know if this is right!
vector<unsigned> coarsen_AMG(const list<image_level>::iterator it, 
      const seg_params& params) {
// {{{
// Coarsen the graph using something similar to the "standard" AMG coarsener
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

      // XXX untested, but we're not using it
      //double max_elem = *max_element(it->A.row_ptr.begin() + it->A.row_ptr[i],
      //      it->A.row_ptr.begin() + it->A.row_ptr[i+1]);

      // now add elements to A_bar
      for (unsigned jp = it->A.row_ptr[i]; jp < it->A.row_ptr[i+1]; ++jp) {
         unsigned j = it->A.col_ind[jp];
         double Aij = it->A.val[jp];

         if ( j != i && Aij >= params.theta*row_sum ) {
            ++row_ptr[i+1];
            col_ind.push_back(j);
            val.push_back(Aij);

            // count the number of nonzeros in each column
            lambda[j] += 1;
         }
      }
   }

   matrix_crs<double> A_bar(row_ptr, col_ind, val, M, M, 1);
   row_ptr.resize(0);
   col_ind.resize(0);
   val.resize(0);

   
   // start assigning C/F nodes
   // T[i] = 0  <- unassigned
   // T[i] = 1  <- C-point
   // T[i] = 2  <- F-point
   vector<unsigned> T(M,0);

   // if salient, designate as C node
   for (unsigned i = 0; i < M; ++i) {
      if ( it->Gamma[i] < params.gamma ) {
         T[i] = 1;
         lambda[i] = 0;
      }
   }

   // while any of the T[i] are zero
   // equivalently, while not all T[i] > 0
   // TODO: this can probably be optimized, but we should make sure it works
   // before we mess around with that

   unsigned j, max_lambda;
   vector<unsigned> K;
   while ( any_of(T.begin(), T.end(), 
            [] (unsigned Ti) -> bool {return Ti == 0;}) ) {

      max_lambda = *max_element(lambda.begin(), lambda.end());
      for ( j = 0; j < M; ++j ) {
         if ( T[j] == 0 && lambda[j] == max_lambda) break;
      }

      if ( j == M ) {
         cout << "T = " << endl;
         print_vector(T);

         cout << "lambda = " << endl;
         print_vector(lambda);
 
         exit(-1);
      }

      T[j] = 1;
      lambda[j] = 0;

      vector<unsigned> K = strongly_influenced_by_j(A_bar, T, j);
      //print_vector(K);

      vector<unsigned> H;
      for (auto k_it = K.begin(); k_it != K.end(); ++k_it) {
         
         T[*k_it] = 2; // nodes in K become F-points
         lambda[*k_it] = 0;

         H = strongly_influence_k(A_bar, T, *k_it);
         for (auto h_it = H.begin(); h_it != H.end(); ++h_it) {
            ++lambda[*h_it];
         }
      }
   }

   // Assign C-points
   vector<unsigned> C;
   for (unsigned i = 0; i < T.size(); ++i) {
      if ( T[i] == 1 ) C.push_back(i);
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
   
   vector<unsigned> K;

   //// Look down columns (which is slow in CRS!)
   //// We check if T_k == 0 first, because that is much cheaper than 
   //// checking if A_bar[k,j] > 0 (think "constraint propagation")
   // XXX: if done this way, this produces a bottleneck in the code
   //for (unsigned k = 0; k < A_bar.m; ++k) {

   //   if ( T[k] != 0 ) continue;

   //   // we now know T[k] == 0, so check if A_bar[k,j] > 0
   //   for (unsigned cp = A_bar.row_ptr[k]; cp < A_bar.row_ptr[k+1]; ++cp) {
   //      unsigned c = A_bar.col_ind[cp];

   //      if ( c < j ) continue;
   //      else if ( c == j ) { // value exists at [k,j] in A_bar
   //         if ( A_bar.val[cp] > 0 ) {
   //            //cout << "[k,j] = [" << k << "," << j << "]" << endl;
   //            K.push_back(k);
   //            break;
   //         }
   //      }
   //      else { // c > j, so value doesn't exist in matrix
   //         break;
   //      }
   //   }
   //}
   
   // A_bar is symmetric, so we could look down rows instead of columns
   // This should be faster in CRS, but only works for symmetric matrices
   //matrix_crs<double> Z = A_bar - transpose(A_bar);
   //Z.clean();
   //cout << Z << endl;
   //exit(-1);
   
   // Look down the jth row of A_bar, which is equivalent to checking down
   // the kth column, since A_bar is symmetric
   // We check if T_k == 0 first, because that is much cheaper than 
   // checking if A_bar[k,j] > 0 (think "constraint propagation")
   for (unsigned kp = A_bar.row_ptr[j]; kp < A_bar.row_ptr[j+1]; ++kp) {

      unsigned k = A_bar.col_ind[kp];
      if ( A_bar.val[kp] > 0 && T[k] == 0 ) {
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

      if ( T[h] == 0 && A_bar.val[hp] > 0 ) {
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
matrix_crs<double> image_seg(const string& img_filename, seg_params& params) {
// {{{
   cv::Mat image;

   image = load_image(img_filename);
   valarray<double> I = image_to_intensity(image, params);
  
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
// this uses default segmentation parameters
matrix_crs<double> image_seg(const string& img_filename) {
// {{{
   // (not necessarily good...) default parameters 
   // eqn 21 of Inglis et al. 2010
   seg_params params(100.,100.,100.,0.1,0.1,0.15,5,1);

   return image_seg(img_filename, params);
// }}}
}

// Helpers
//////////
void assign_uniquely(matrix_crs<double>& U) {
// {{{
   for (unsigned i = 0; i < U.m; ++i) {

      // if there is only one element in this row, the pixel assignment is
      // unique
      if ( U.row_ptr[i+1] - U.row_ptr[i] == 1 ) continue;
      else {
         // max_element returns (an iterator pointing to) the first
         // maximal element
         auto max_el_it = max_element(U.val.begin() + U.row_ptr[i],
               U.val.begin() + U.row_ptr[i+1]);
         unsigned max_jp = distance(U.val.begin(), max_el_it);
         U.val[max_jp] = 1.;

         // This is similar to matrix_crs::clean
         // erase any elements other than max_jp
         unsigned num_removed = 0;
         unsigned jp = U.row_ptr[i];
         unsigned jp_end = U.row_ptr[i+1];
         while ( jp < jp_end ) {
            if ( jp != max_jp ) {
               U.col_ind.erase(U.col_ind.begin() + jp);
               U.val.erase(U.val.begin() + jp);
               ++num_removed;
               --jp_end;
               --max_jp;
            }
            else ++jp;
         }
         
         // fix row pointers
         for (unsigned r = i+1; r <= U.m; ++r) {
            U.row_ptr[r] -= num_removed;
         }
      }
   }
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

cv::Mat load_image(const string& img_filename) {
// {{{
// Load image from file specified by filename into a CV Mat with uchar entries
// The image is blended to a grayscale image and sacled so the pixels
// have values between 0 and 1
// Matrices must be size n x n for the image segmentation method
   cv::Mat img;
   img = cv::imread(img_filename, 0); // load image to grayscale; uchar vals

   // Check to see if we read the image properly
   if ( img.data == NULL ) {
      cerr << "error: image_seg.cpp:load_image: OpenCV imread loaded empty "
           << "matrix.  filename = " << img_filename << endl;
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
   
   // Set the size of the image in the parameters class
   // We'll use it later to build A on the first level
   params.n = static_cast<unsigned>(img.rows);

   for (int i=0; i < img.rows; ++i) {

      const uchar* img_i = img.ptr<uchar>(i); // pointer to row i
      for (int j=0; j < img.cols; ++j) {
         // scale to [0,1] from range [0, 255] (uchar)
         I[img.rows*j + i] = static_cast<double>( img_i[j] )/ 255.;
      }
   }
  
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


void write_seg_images(const matrix_crs<double>& U, const unsigned n, 
      const string& filename_base, const unsigned mode) {
// {{{
// Write the segments specified by the matrix U to files with base filename
// given by filename_base.  The files will be something like
//    filename_base_seg_01.jpg
//    filename_base_seg_02.jpg
//    filename_base_seg_03.jpg
//    ...
// 
// Each image should be size n x n
//
// mode == 0 - write images with just the segments, no original image
// mode == 1 - blend each segment with the original image 
//             TODO

   vector<valarray<double>> intensities = mat_to_vecs(U);

   if ( mode == 0 ) {

      for (unsigned i = 0; i < intensities.size(); ++i) {
         string filename = filename_base;
         filename += "_seg_";
         // add leading zeros
         if ( i > 10 ) filename += "0";
         else filename += "00";
         filename += to_string(i);
         filename += ".png";

         try {
            cout << "writing image file " << filename << endl;
            cv::imwrite(filename, 
                  intensity_to_image1(intensities[i], n));
         }
         catch ( runtime_error& ex ) {
            cerr << "error: image_seg.cpp:write_seg_images: Something went "
                 << "wrong while writing images." << endl 
                 << "mode = " << mode << endl
                 << "exception: " << ex.what();
            exit(-1);
         }
      }
   }

   else if ( mode == 1 ) {
      cerr << "error: not implemented" << endl;
      exit(-1);
   }

   else {
      cerr << "error: image_seg.cpp:write_segment_images: Bad mode "
           << "specified: mode = " << mode << endl;
      exit(-1);
   }
// }}}
}

// }}}


int main(void) {

   //seg_params params(10.,100.,100.,0.1,0.1,0.15,2,1);
   //matrix_crs<double> U = image_seg("square.png", params);
   //write_seg_images(U, params.n, "gen_imgs/square", 0);

   seg_params params(10.,100.,100.,0.1,0.1,0.15,3,1);
   matrix_crs<double> U = image_seg("squares.png", params);
   write_seg_images(U, params.n, "gen_imgs/squares", 0);
   
   //seg_params params(10.,10.,10.,0.01,0.1,0.15,5,1);
   //matrix_crs<double> U = image_seg("peppers_25.jpg", params);
   //write_seg_images(U, params.n, "gen_imgs/peppers_25", 0);
   
   return 0;
}
