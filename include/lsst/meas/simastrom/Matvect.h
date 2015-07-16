#ifndef MATVECT__H
#define MATVECT__H

#include <stdlib.h> // for abort
#include <string.h> // for memset
#include <iostream>
#include <string>

namespace lsst {
namespace meas {
namespace simastrom {

class Vect;
class Mat;

#ifndef MATVECT_CHECK_BOUNDS
#define MATVECT_CHECK_BOUNDS
#endif

// Routines for solving linear systems + inversion 
//==================================================================



// solving linear system A.X = B 
// Uses lapack dposv_
// Matrix A is assumed to be symmetric (you'll get a core dump if it is not
// (actually this is not compulsory with dposv (see 'man dposv') but we do it
// for consistency).
// You just need to fill half (n*(n+1)/2 parameters) of the matrix
// if you have filled matrix M parameters for which x>=y (with M(x,y)=...), use UorL="L"
// else use UorL="U"
// Matrix A is modified in this routine (also B which is at the end the solution X)
int cholesky_solve(Mat &A, Vect &B, const char* UorL ); // "L" <=> A(i,j)!=0 for i>=j 

//! same as above but for a set of right hand sides.
int cholesky_solve(Mat &A, Mat &B, const char* UorL);

/* Inverts matrix A using the factorization done in cholesky_solve
 Uses lapack dpotri_
 Matrix A is assumed to be symmetric (you'll get a core dump if it is not
 (actually this is not compulsory (see 'man dptri') but we do it
 for consistency).
 This routine must be called after cholesky_solve, make sure the value of UorL
 is the same as that used with dposv.
 The MATRIX IS SYMMETRIZED on exit 
*/

int cholesky_invert(Mat &A, const char* UorL = "L"); // when cholesky_solve is called first



/* Routine to "solve" a **posdef** system with constraints.
   The found solution minimizes:
   X**t A X  -  B**t X   under the constraints : C**t X = Y

   - dimensions : A(n,n), B(n), C(n,m), Y(m)
   - return value : 0 == OK
   - On exit the solution is in B and the Lagrange multiplier values in Y

   The technique used assumes that the unconstrained problem is posdef
   (because Cholesky factorization is used).  In particular, if the
   constraints are mandatory for the problem to have a solution at
   all, this is the wrong routine.  You then have e.g. to build your
   constraints into the matrix and call symetric_solve.
*/


//! calls cholesky_solve with a dummy RHS and then cholesky_invert.
int posdef_invert(Mat &A, const char* UorL);


int cholesky_solve_with_constraints(Mat &A, Vect &B, 
				    Mat &Constraints, Vect& Y,
				    const char* UorL);


/*! When the system is posdef upto NCut. typically this happens with Lagrange
   Mulitpliers. No inversion routine (yet?) but it is possible */
//! solving linear systems posdef upto NCut.
int cholesky_solve_quasi_posdef(Mat &A, Vect &B, 
				const unsigned &NCut, const char* UorL);

//! same as above for several RHS
int cholesky_solve_quasi_posdef(Mat &A, Mat &B, 
				const unsigned &NCut, const char* UorL);

//! solve a symetric system (if posdef, use cholesky_solve)
int general_solve(Mat& A, Vect& B, bool invert_A, const char* UorL);


//! solve a general system
int really_general_solve(Mat &A, Vect &B);


//! invert a general matrix
int really_general_invert(Mat &A);



//! Symetric general system. return 0 when  OK. Undocumented invertion routine
int symetric_solve(Mat& A, Vect& B, const char* UorL);

//! same when there are several RHS (i.e. we have several B's)
int symetric_solve(Mat& A, Mat& B, const char* UorL);

//! eigenvalues using dsyevd. A is destroyed
int symetric_eigenvalues(Mat& A, Vect& EigenVals, const char* UorL);

//! eigenvectors and eigenvalues using dsyevd. A contains eigenvectors on output.
/* the eigenvector corresponding to the lowest eigenvalues reads v(i) = A(i,0) */
int symetric_diagonalize(Mat& A, Vect& EigenVals, const char* UorL);




// Diagonalization cernlib routine
//====================================================================
int DiagonalizeRealSymmetricMatrix(const Mat &matrix, Mat &eigenvectors , Vect &eigenvalues);

// Mat and Vect classes
//====================================================================


class Mat {
  
 private:
  double *data;
  unsigned int nx,ny;
  
  
 public:
  
  Mat() : data(NULL), nx(0), ny(0) {};  
  Mat(const unsigned int NX, const unsigned int NY);
  Mat(const Mat& other);
  ~Mat() { delete [] data;}
  
  void allocate(const unsigned int NX, const unsigned int NY);
  
  double operator () (const unsigned int i, const unsigned int j) const
  {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=nx || j>=ny) { 
    std::cout << "Mat::operator () overflow i,j nx,ny " 
	 << i << "," << j << " "
	 << nx << "," << ny << " "
	      << std::endl;
    abort();
  }
#endif
  return data[i+j*nx];
}

  double& operator () (const unsigned int i, const unsigned int j)
  {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=nx || j>=ny) { 
    std::cout << "Mat::operator () overflow i,j nx,ny " 
	 << i << "," << j << " "
	 << nx << "," << ny << " "
	      << std::endl;
    abort();
  }
#endif
  return data[i+j*nx];
}

  
  unsigned int SizeX() const { return nx;}
  unsigned int SizeY() const { return ny;}
  
  void Zero() {memset(data, 0, nx*ny*sizeof(double));};
  void Identity();

  const double* Data() const {return data;};
  double* NonConstData() {return data;};

  
  friend std::ostream& operator << (std::ostream &stream, const Mat &m)
    { m.writeASCII(stream); return stream;}
  
  // get a block of this matrix as a new matrix
  // size (x_max-x_min+1)*(y_max-y_min+1) (both min and max included)
  // remember  0 <= x < nx ,  0 <= y < ny
  Mat SubBlock
    (unsigned int x_min,unsigned int x_max,unsigned int y_min,unsigned int y_max) const;
  
  Mat WithoutRows(unsigned int y_min,unsigned int y_max) const;
  Mat WithoutColumns(unsigned int x_min,unsigned int x_max) const;


  // i/o in fits for matrices
  int readFits(const std::string &FitsName);
  int writeFits(const std::string &FitsName) const;
  int readASCII(std::istream& Stream);
  int readASCII(const std::string &FileName);
  int writeASCII(std::ostream& Stream) const;
  int writeASCII(const std::string &FileName) const;
  

  //! call it with the same argument as cholesky_sove or cholesky_invert
  void Symmetrize(const char* UorL);
  
  // inverts a posdef matrix using Cholesky factorization from lapack.
  int CholeskyInvert(const char *UorL);

  // operators
  Mat operator +(const Mat& Right) const;
  Mat operator -(const Mat& Right) const;
  Mat operator *(const Mat& Right) const;
  Vect operator *(const Vect& Right) const;
  Mat & operator =(const Mat& Right);
  
  void operator +=(const Mat& Right);
  void operator -=(const Mat& Right);
  void operator *=(const Mat& Right);
  
  Mat operator *(const double Right) const;
  friend Mat operator *(const double Left, const Mat &Right);
  void operator *=(const double Right);

  operator double() const;
  operator Vect() const;
  Mat transposed() const;
  void ExtractSubMat(const unsigned IndexMapping[], const unsigned NewSize,
		     Mat &SubMat) const;


};

class Vect {

 private:
  
  double *data;
  unsigned int n;
  
  
 public:
  
  Vect() : data(NULL), n(0) {};
  Vect(const unsigned int N);
  Vect(const Vect&);
  ~Vect() { delete [] data;}

  void allocate(const unsigned int N);
  
  double operator () (const unsigned int i) const
  {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=n) {
    std::cout << "Vec::operator () overflow i,n " 
	      << i << "," << n << std::endl;
    abort();
  }
#endif
  return data[i];
  }


  double& operator () (const unsigned int i)
  {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=n) {
    std::cout << "Vec::operator () overflow i,n " 
	      << i << "," << n << std::endl;
    abort();
  }
#endif
  return data[i];
}  
  unsigned int Size() const { return n;}
  unsigned int size() const { return n;}
  
  void Zero() {memset(data, 0, n*sizeof(double));};

  const double* Data() const {return data;};
  double* NonConstData() {return data;};

  int writeASCII(const std::string &FileName) const;
  int writeASCII(std::ostream& Stream) const;
  int readASCII(std::istream& Stream);


  void dump(std::ostream& Stream) const;
  friend std::ostream& operator << (std::ostream &stream, const Vect &v)
    { v.dump(stream); return stream;}

  // operators
  Vect operator +(const Vect& Right) const;
  Vect operator -(const Vect& Right) const;
  Vect & operator =(const Vect& Right);
  
  double operator *(const Vect& Right) const; // scalar product
  
  void operator +=(const Vect& Right);
  void operator -=(const Vect& Right);
  
  

  Vect operator *(const double Right) const;
  friend Vect operator *(const double Left, const Vect &Right);
  void operator *=(const double Right);

  Mat transposed() const;
  operator Mat() const;
  //  operator double() const;
};

}}}

#endif /*MATVECT__H */
