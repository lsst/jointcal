#include <iostream>
#include <iomanip> // for setprecision()
#include <memory> // for auto_ptr 
#include <fstream>
#include <fitsio.h>
#include <assert.h>
#include <string.h> // memcpy
#include "lsst/meas/simastrom/Matvect.h"

using namespace std;

#define dfact dfact_
#define dfinv dfinv_
#define dfeqn dfeqn_


// #define eisrs1d eisrs1d_

// // using cernstuff (from cernlib)
// extern "C" 
// {
//   void dfact(int *N, double *A, int *idim, double *r, int *ifail, double *Det, int *jfail);
//   void dfinv(int *n, double *A, int *idim, double *r);
//   void dfeqn(int *n, double *a, int *idim, double *r, int *k, double *b);
//   void eisrs1d(int *NM,int *N,double *AR,double *WR,double *ZR,int *IERR,double *WORK);
// }

// using lapack 
extern "C" {
  void dposv_(char *, int *, int *, double *, int *, double *, int *, int *);
  void dpotri_(char *, int *, double *, int *, int *);
  void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
  void dpotrs_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, 
	       double *B, int *LDB, int *INFO);

  void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
  void dsysv_(const char*, int* n, int* nrhs, 
	      double* a, int* lda, int* ipiv, 
	      double* b, int* ldb, 
	      double* work, int* lwork, int* info);
  void dsyevd_(char* JOBZ, char* UPLO, int* N,
	       double* A, int* LDA, 
	       double* W, 
	       double* WORK, int* LWORK, int* IWORK, int* LIWORK,
	       int* INFO);
  void dsytri_(const char*, int* n, double* a , int* lda, int* ipiv, double* work, int* info);

  void dgetri_(int *n, double *A, int * lda, int *ipiv, double *work, int *lwork, int *info);
  void dgetrf_(int *n, int* m, double *A, int *lda, int *ipiv, int *info);


#ifdef STORAGE
  void dgees_(char *jobvs, char *, void *select , 
	      int *n, double *a, int *lda, int * sdim, 
	      double *wr, double *wi, double *vs, int *ldvs, 
	      double * work, int *lwork, bool *bwork, int *info);

#endif
	      

};

// //==========================================================================
// int DiagonalizeRealSymmetricMatrix(const Mat &matrix, Mat &eigenvectors , Vect &eigenvalues) {
//   if ( matrix.SizeX() != matrix.SizeY() ) {
//     cerr << "in  DiagonalizeRealSymmetricMatrix, matrix is not square !!" << endl;
//     return 12;
//   }
//   int n = matrix.SizeX();
//   if ( int(eigenvectors.SizeX()) != n || int(eigenvectors.SizeY()) )
//     eigenvectors.allocate(n,n);
//   if ( int(eigenvalues.Size()) != n )
//     eigenvalues.allocate(n);
  
//   double * work = new double[n];
//   int ierr = 0;
  
//   Mat matrix_copy = matrix;
//   eisrs1d(&n,&n,matrix_copy.NonConstData(),eigenvalues.NonConstData(),eigenvectors.NonConstData(),&ierr,work);
  
//   delete [] work;
//   return ierr;
// }


//==========================================================================
namespace lsst {
namespace meas {
namespace simastrom {

int cholesky_solve(Mat &A, Vect &B, const char* UorL)
{
#ifdef FNAME
  cout << " >  cholesky_solve" << endl;
#endif  

  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_solve Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  if(B.Size() != A.SizeX()) {
    cout << "error in matvect.cc, cholesky_solve Vector B must have a dimension B.Size()=A.SizeY() and you have B.Size(),A.SizeY() = "
	 << B.Size() << "," << A.SizeY() << endl;
    abort();
  }

  double *a = A.NonConstData();
  double *b = B.NonConstData();
  
  int n = B.Size();
  int nhrs = 1, info = 0;
  char uorl[6]; strncpy(uorl,UorL,6);


  dposv_(uorl, &n, &nhrs, a, &n, b, &n, &info);

  if (info != 0) 
    cerr << " cholesky_solve(" << a << "," << b << "," << n
	 << ") : Warning: Cholesky factorization failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}

int cholesky_solve(Mat &A, Mat &B, const char* UorL)
{
#ifdef FNAME
  cout << " >  cholesky_solve" << endl;
#endif  

  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_solve Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  if(B.SizeX() != A.SizeX()) {
    cout << "error in matvect.cc, cholesky_solve Vector B must have a dimension B.SizeX()=A.SizeY() and you have B.SizeX(),A.SizeY() = "
	 << B.SizeX() << "," << A.SizeY() << endl;
    abort();
  }

  double *a = A.NonConstData();
  double *b = B.NonConstData();  
  int n = B.SizeX();
  int nhrs = B.SizeY(), info = 0;
  char uorl[6]; strncpy(uorl,UorL,6);

  dposv_(uorl, &n, &nhrs, a, &n, b, &n, &info);

  if (info != 0) 
    cerr << " cholesky_solve(" << a << "," << b << "," << n
	 << ") : Warning: Cholesky factorization failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}


//assumes that the matrix has been factorized by a call to cholesky_solve.
int cholesky_invert(Mat &A, const char* UorL)
{  
  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_invert Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  
  int info = 0;
  double *a = A.NonConstData();
  int n = A.SizeX();

  //  Now invert using the factorization done in dposv_

  char uorl[6]; strncpy(uorl,UorL,6);

  dpotri_(uorl, &n, a, &n, &info);

  if (info != 0) 
    cerr << " cholesky_invert(" << a << "," << n
	 << ") : Warning: inversion failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  else A.Symmetrize(UorL);
  return info;
}


int symetric_eigenvalues(Mat& A, Vect& EigenVals, const char* UorL)
{
  assert(A.SizeX() == A.SizeY());
  if (EigenVals.Size() < A.SizeX()) EigenVals.allocate(A.SizeX());
  
  int  info = 0;
  char* jobz = "N";
  int  n = A.SizeX();
  char uorl[6]; strncpy(uorl,UorL,6);
  double* a = A.NonConstData();
  double* eval  = EigenVals.NonConstData();
  int     lwork = 1 + 6*n + 2*n*n;
  double* work = new double[lwork];
  int     liwork = 3 + 5*n;
  int*    iwork = new int[liwork];
  
  dsyevd_(jobz, uorl, &n,
	  a, &n, 
	  eval, 
	  work, &lwork, iwork, &liwork, 
	  &info);

  delete[] work;
  delete[] iwork;

  
  if(info != 0)
    {
      cout << "[symetric_eigenvalues] ERROR in dsysevd, info=" << info
	   << endl;
    }
  
  return info;
}

int symetric_diagonalize(Mat& A, Vect& EigenVals, const char* UorL)
{
  assert(A.SizeX() == A.SizeY());

  int  info = 0;
  char* jobz = "V";
  int  n = A.SizeX();
  char uorl[6]; strncpy(uorl,UorL,6);
  double* a = A.NonConstData();
  double* eval  = EigenVals.NonConstData();
  int     lwork = 1 + 6*n + 2*n*n;
  double* work = new double[lwork];
  int     liwork = 3 + 5*n;
  int*    iwork = new int[liwork];
  
  dsyevd_(jobz, uorl, &n,
	  a, &n, 
	  eval, 
	  work, &lwork, iwork, &liwork, 
	  &info);
  
  if(info != 0)
    {
      cout << "[symetric_diagonalize] ERROR in dsysevd, info=" << info
	   << endl;
    }
  delete[] work;
  delete[] iwork;
  return info;
  
}


#ifdef STORAGE 
int posdef_invert(Mat &A, const char* UorL)
{
  Vect B(A.SizeX());
  int info = cholesky_solve(A, B, UorL);
  if (info == 0)
    info = cholesky_invert(A, UorL);
  return info;
}
#endif

#define DO10(A) A;A;A;A;A;A;A;A;A;A;

static double fast_scal_prod(double *x, double *y, const int size)
{
int nblock = size/10;
int remainder = size - nblock*10;
double sum = 0;
for (int i=0; i<nblock; ++i) 
  {
  DO10( sum+= (*x)*(*y); ++x; ++y;)
  }
 for (int i=0; i<remainder; ++i) {sum+= (*x)*(*y); ++x; ++y;}
return sum;
}



/*  minimum of 1/2 xT A X - Bx under constraint "Constraint x = ConstraintValue "
   minimum of 1/2 xT A X - Bx under constraint   CT x = Y

   lambda is solution of Ax - B + C lambda = 0;


   compute x' = A**(-1)*B
   compute lambda :   (CT * A**(-1)C) Lambda  = CT *A**(-1)*B - Y
   compute x = A**(-1)*(B - C*lambda)

you may have a look at (which uses different notations):
http://www-fp.mcs.anl.gov/Otc/Guide/OptWeb/continuous/constrained/qprog/


*/



int cholesky_solve_with_constraints(Mat &A, Vect &B, 
				    Mat &Constraints, Vect& ConstraintsValues,
				    const char* UorL)
{
#ifdef FNAME
  cout << " >  cholesky_solve_with_constraints" << endl;
#endif  

  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_solve_with_constraints Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  if(B.Size() != A.SizeX()) {
    cout << "error in matvect.cc, cholesky_solve_with_constraints " << endl 
	 << " Vector B must have a dimension B.Size()=A.SizeY() and you have B.Size(),A.SizeY() = "
	 << B.Size() << "," << A.SizeY() << endl;
    abort();
  }

  if (Constraints.SizeX() != A.SizeX())
    {
      cout << "error in matvect.cc, cholesky_solve_with_constraints " << endl
	 << " Constraints.SizeY() != A.Size  : "
	 << Constraints.SizeY() << ' ' << A.SizeX() << endl;
    abort();
    }

  if (Constraints.SizeY() != ConstraintsValues.Size())
    {
      cout << "error in matvect.cc, cholesky_solve_with_constraints " << endl
	   << " Constraints.SizeY() != ConstraintsValues.Size() : "
	   << Constraints.SizeY() << ' ' << ConstraintsValues.Size() << endl;
      abort();
    }

  double *a = A.NonConstData();
  double *b = B.NonConstData();
  
  int n = B.Size();
  int nrhs = 1, info = 0;
  char uorl[6]; strncpy(uorl,UorL,6);


  dposv_(uorl, &n, &nrhs, a, &n, b, &n, &info);

  if (info != 0) 
    {
      cerr << " cholesky_solve_with_constraints(" << a << "," << b << "," << n
	   << ") : Warning: Cholesky factorization failure . info =" 
	   << info <<  " (>0 is not pos.def)" << endl;
      
      return info;
    }
  

  unsigned int nconst = Constraints.SizeY();
  //need to keep a copy of the constraints
  Mat am1c(Constraints);
  double *c = am1c.NonConstData();
  nrhs = nconst;
  dpotrs_(uorl, &n, &nrhs, a, &n, c, &n, &info);

  Mat smalla(nconst,nconst);
  Vect smallb(nconst);
  for (unsigned int j=0; j <nconst; ++j)
    {
      for (unsigned int i = 0; i<=j ; ++i)
	{
	  smalla(i,j) = fast_scal_prod(&Constraints(0,i), &am1c(0,j), n);
	}
      smallb(j) = fast_scal_prod(&Constraints(0,j),&B(0), n) - ConstraintsValues(j);
    }

  // solve the small system
  if (cholesky_solve(smalla, smallb,"U") != 0)
    {
      cout << " cholesky_solve_with_constraints small sub-system could not be solved " << endl;
      return info;
    }

  // small b now contains the values of Lagrange multiplier

  B -=  am1c.transposed()*smallb;

  cout << " cholesky solve with constraints : residuals " << Constraints*B << endl;
  return 0;

}


// routines to solve problems where a subpart of A (0::NCut-1) is posdef,
// using cholesky for this part and other means for the remainder.

int cholesky_solve_quasi_posdef(Mat &A, Vect &B, 
				const unsigned &NCut, const char* UorL)
{
  unsigned s = B.Size();
  Mat VB(s,1);
  memcpy(VB.NonConstData(), B.NonConstData(), s*sizeof(double));
  int info = cholesky_solve_quasi_posdef(A, VB, NCut, UorL);
  memcpy(B.NonConstData(), VB.NonConstData(), s*sizeof(double));
  return info;
}



int cholesky_solve_quasi_posdef(Mat &A, Mat &B, 
				const unsigned &NCut, const char* UorL)
{
#ifdef FNAME
  cout << " >  cholesky_solve_quasi_posdef" << endl;
#endif

  if ( (*UorL != 'U') && ( *UorL != 'L'))
    {
      cout << "error in matvect.cc, cholesky_solve_quasi_posdef : " << endl
	   << " UorL should be either \"U\" or \"L\". got : " << UorL << endl;
      return -1;
    }

  if (A.SizeX() != A.SizeY()) {
    cout << "error in matvect.cc, cholesky_solve_quasi_posdef " << endl 
	 << " A must be square "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  if(B.SizeX() != A.SizeX()) {
    cout << "error in matvect.cc, cholesky_solve_quasi_posdef " << endl 
	 << " Vector B must have a dimension B.Size()=A.SizeY() and you have B.Size(),A.SizeY() = "
	 << B.SizeX() << "," << A.SizeY() << endl;
    abort();
  }

  if (NCut >= A.SizeX())
    {
      cout << "error in matvect.cc, ccholesky_solve_quasi_posdef " << endl
	 << " NCut >= A.Size() : NCut, A.Size  : "
	 << NCut << ' ' << A.SizeX() << endl;
    abort();
    }

  double *a = A.NonConstData();
  double *b = B.NonConstData();
  
  int n = B.SizeX();
  int posDefSize = NCut;
  unsigned int nrhs = B.SizeY();
  int inrhs = nrhs;
  int info = 0;
  int ldb = B.SizeX();
  char uorl[6]; strncpy(uorl,UorL,6);


  dposv_(uorl, &posDefSize, &inrhs, a, &n, b, &ldb, &info);

  if (info != 0) 
    {
      cerr << " cholesky_solve_quasi_posdef(" << a << "," << b << "," << n
	   << ") : Warning: Cholesky factorization failure . info =" 
	   << info <<  " (>0 is not pos.def)" << endl
	   << " NCut argument may be too large " << endl;      
      return info;
    }
  
  /*
    write the problem:

   (   D   E ) (X1)  = (B1)
   ( tr(E) F ) (X2)    (B2)

   where D is posdef (dim NCut)
   we have 
      (tr(E)D**(-1) E - F) X2 = tr(E)D**(-1)B1-B2
    - once D is factorized, we compute D**(-1) B1 (all done by dposv)
    - we save E.
    - we compute D**(-1) E
    - we compute G = tr(E)D**(-1) E - F
    - we compute tr(E)D**(-1)B1-B2
    - we solve for X2 (using dsysv)
    - we compute X1 =  D**(-1) B1 - D**(-1) E X2
    - done.


    For storing the factorization in view of being able to provide an inverse:
    I think that storing 
      (  Chol(D)          E            )
      (  ...      Fact(Ed**(-1)E -F)   ) 

    would enable subsequent inversion.

  */



  unsigned int nadd = n-NCut;
  //need to keep a copy of E. Extracting it from A is not mandatory...
  Mat E(NCut, nadd);
  if (*UorL == 'U')
    {
      for (unsigned j=0; j < nadd; ++j)
	for (unsigned i=0; i<NCut; ++i)
	  E(i,j) = A(i,j+NCut);
    }
  else
    {
	for (unsigned i=0; i<NCut; ++i)
	  for (unsigned j=0; j < nadd; ++j)
	    E(i,j) = A(j+NCut,i);
    }      

  Mat Dm1E(E);

  double *pE = Dm1E.NonConstData();
  int nrhs2 = nadd;
  dpotrs_(uorl, &posDefSize, &nrhs2, a, &n, pE, &posDefSize, &info);

  Mat G(nadd,nadd);
  Mat smallb(nadd,nrhs);

  //We decide arbitrarily that G is going to be of the "U" type
  for (unsigned int j=0; j < nadd; ++j)
    {
      for (unsigned int i = 0; i<=j ; ++i)
	G(i,j) = fast_scal_prod(&E(0,i), &Dm1E(0,j), NCut);
      for (unsigned k=0; k<nrhs; ++k)
	smallb(j,k) = fast_scal_prod(&E(0,j),&B(0,k), NCut) - B(j+NCut,k);
    }

      // we have to read the right part of of the subpart of A to put into F
  if (*UorL == 'U')
    {
      for (unsigned int j=0; j < nadd; ++j)
	for (unsigned int i = 0; i<=j ; ++i)
	  G(i,j) -= A(i+NCut,j+NCut);
    }
  else
    {
      for (unsigned int i = 0; i<nadd ; ++i)
	for (unsigned int j=i; j <nadd; ++j)
	  G(i,j) -= A(j+NCut,i+NCut);
    }
    

  // solve the small system
  
  if ((info = symetric_solve(G, smallb,"U")) != 0)
    {
      cout << " cholesky_solve_quasi_posdef: small sub-system could not be solved " << endl;
      return info;
    }

  // small b now contains X2
  Mat bprime = smallb*Dm1E;
      
  for (unsigned j=0; j<nrhs; ++j)
    {
      for (unsigned k=0; k<NCut; ++k)
	B(k,j) -= bprime(k,j);
      for (unsigned k=0; k<nadd; ++k)
	B(k+NCut,j) = smallb(k,j);
    }
  
  // store the partially inverted matrices in A
  // 
  //     ( D^-1    D^-1 E )
  // A = ( ...     (Et D^-1 E - F)^-1
  // 
  // will do that later.
  
  return 0;
}


// call this *after* a cholesky_solve_quasi_posdef (!)
int cholesky_invert_quasi_posdef(Mat& A, const unsigned& NCut, const char* UorL)
{
  // from the matrix A given by cholesky_solve_quasi_posdef, we can compute A^-1
  // 
  // 
  // A^-1 = ( D^-1 - D^-1 E (Et D^-1 E - F)^-1 Et D^-1     D^-1 E ()^-1 ) 
  //        (  ...                                         -(Et D^-1 E - F)^-1
  
  // or, if A_pre_invert =  (U  V)
  //                        (Vt W)
  // we have:
  // 
  //        A = (U     U-VWVt)
  //            ( ...    -W  )
  // 
  // should not be too complicated 
  // 
  // 
  return 0;
}


// calls dsytri after the use of dsysv
/* this is indeed almost the same routine as symetric_solve !!!! */
int general_solve(Mat& A, Vect& B, bool invert_A, const char* UorL)
{
  int nrhs=1;
  int n = B.Size();
  int lwork=n;
  int info;
  
  double* a = A.NonConstData();
  double* b = B.NonConstData();
  int*    ipiv = new int[n];
  double* work = new double[n];
  
  dsysv_(UorL, &n, &nrhs, 
	 a, &n, ipiv, b, &n, work, &lwork, &info);
  if(info!=0) {
    std::cout << " Pb solving sys: info="
	      << info << std::endl;
    delete[] ipiv;
    delete[] work;
    return info;
  }
  
  if(! invert_A) {
    delete[] ipiv;
    delete[] work;
    return info;
  }
  
  dsytri_(UorL, &n, a , &n, ipiv, work, &info);
  
  if(info!=0) {
    std::cout << " Pb inverting the matrix: info="
	      << info << std::endl;
  }
  
  delete[] ipiv;
  delete[] work;
  return info; 
}

int symetric_solve(Mat& A, Vect& B, const char* UorL)
{
  unsigned s = B.Size();
  Mat tmpB(s,1);
  memcpy(tmpB.NonConstData(), B.NonConstData(), s*sizeof(double));
  int info = symetric_solve(A, tmpB, UorL);
  memcpy(B.NonConstData(), tmpB.NonConstData(), s*sizeof(double));
  return info;
}


// symetric general system. return 0 when  OK
int symetric_solve(Mat& A, Mat& B, const char* UorL)
{
  int nrhs=B.SizeY();
  int n = B.SizeX();
  int lwork=n;
  int info;
    
    double* a = A.NonConstData();
    double* b = B.NonConstData();
    int*    ipiv = new int[n];
    double* work = new double[n];
    char uorl[6]; strncpy(uorl,UorL,6);

    dsysv_(uorl, &n, &nrhs, 
	   a, &n, ipiv, b, &n, work, &lwork, &info);
    
    delete[] ipiv;
    delete[] work;
    
    if(info!=0) {
      std::cout << " symetric_solve (dsysv) Problem factorizing the matrix: info="
		  << info << std::endl;
    }

    return info;
}


int really_general_solve(Mat &A, Vect &B)
{
  int n = B.Size();
  int nrhs = 1;
  int lda = A.SizeX();
  int *ipiv = new int[n];
  int ldb = n;
  int info;

  dgesv_(&n, &nrhs, A.NonConstData(), &lda, ipiv, B.NonConstData(), &ldb, &info);
  delete [] ipiv;

  if (info != 0)
    std::cout << " really_general_solve (dgesv) Problem factorizing the matrix: info="
	      << info << std::endl;

  return info;
}

int really_general_invert(Mat &A)
{
  int m = A.SizeX();
  int n = A.SizeY();
  int lda = A.SizeX();
  int *ipiv= new int[n];
  auto_ptr<int> pipiv(ipiv);

  int info;

  dgetrf_(&m, &n, A.NonConstData(), &lda, ipiv, &info);
  if (info != 0)
    {
      std::cout << " really_general_invert (dgetrf) Problem factorizing the matrix: info="
		<< info << std::endl;
      
      return info;
    }
  double tmpwork; // for workspace query
  int lwork = -1;
  dgetri_( &n, A.NonConstData() , &lda, ipiv, &tmpwork, &lwork, &info);  
  int nb= int(tmpwork);
  lwork = n*nb;
  double* work = new double[lwork];
  auto_ptr<double> pwork(work);
  dgetri_(&n,  A.NonConstData(), &lda, ipiv, work, &lwork, &info);  
  if (info != 0)
    {
      std::cout << " problem in dgetri : info = " << info << endl;
    }
  return info;
}

#ifdef STORAGE
/* eigenvalues of a general (i.e. non-symetric) matrix */
int general_eigenvalues(Mat &A, Vect &ER, Vect &EI)
{
  char *jobvs="N";
  char *sort="N";
  void * select = NULL;
  int n = A.SizeX();
  assert(A.SizeY() == n );
  int lda = n;
  int sdim; // 0 on output
  ER.allocate(n);
  EI.allocate(n);
  double *vs = NULL;
  int ldvs = 0;
  double work[4*n];
  int lwork = 4*n;
  bool *bwork = NULL;
  int info;
  dgees_(jobvs, sort, select, &n, A.NonConstData(),
	 &lda, &sdim, ER.NonConstData(), EI.NonConstData(),
	 vs, &ldvs, work, &lwork, bwork, &info);
  if (info != 0)
    std:: cout <<  "general_eigenvalues : something went wrong " << endl;
  return info;
}
#endif
  
  

  



//==================================================================




Mat::Mat(const unsigned int NX, const unsigned int NY) { 
 data=NULL;
 nx=ny=0;
 if(NX<0 || NY<0) {
   cout << "Mat::Mat ERROR NX,NY =" << NX << "," << NY <<  endl;
 }else{
   allocate(NX,NY);
 }
}

Mat::Mat(const Mat& other) {
  data=NULL;
  nx=ny=0;
  allocate(other.SizeX(),other.SizeY());
  memcpy(data,other.Data(),sizeof(double)*nx*ny); 
}

void Mat::allocate(const unsigned int NX, const unsigned int NY) {
  if(NX<0 || NY<0) {
    cout << "Mat::allocate ERROR NX,NY =" << NX << "," << NY <<  endl;
  }
  if(nx!=NX || ny!=NY) {
    nx=NX; 
    ny=NY; 
    if (data) 
      delete [] data;
    if(nx*ny>0)
      data = new double[nx*ny];
    else
      data = 0;
  } 
  Zero();
}


int Mat::writeASCII(const std::string &FileName) const
{
  std::ofstream S(FileName.c_str());
  if (!S)
    {
      cout << " Mat::writeASCII() : cannot open " << FileName << endl;
      return 0;
    }
  int status = writeASCII(S);
  S.close();
  return status;
}


int Mat::writeASCII(ostream& Stream) const {
  Stream << nx << " " << ny << endl;
  int oldprec = Stream.precision();
  Stream << setprecision(10);
  for(unsigned int j=0;j<ny;j++) {
    //Stream << "0.." << nx-1 << "," << j << ": ";
    for(unsigned int i=0;i<nx;i++)
      Stream << " " << (*this)(i,j);
    Stream << endl;
  }
  Stream << setprecision(oldprec);
  return 0;
}

int Mat::readASCII(const std::string &FileName)
{
  std::ifstream S(FileName.c_str());
  if (!S)
    {
      cout << " Mat::readASCII() : cannot open " << FileName << endl;
      return 0;
    }
  int status = readASCII(S);
  S.close();
  return status;
}


int Mat::readASCII(istream& Stream) {
  unsigned int  fnx,fny;
  
  Stream >> fnx >> fny;
  allocate(fnx,fny);
  
  double val;
  for(unsigned int j=0;j<ny;j++) {
    for(unsigned int i=0;i<nx;i++) {
      Stream >> val;
      (*this)(i,j)=val;
    }
  }
  return 0;
}


void Mat::Identity() {
  if(nx!=ny) {
    cout << "Mat::Identity ERROR nx!=ny" <<endl;
    abort();
  }
  Zero();
  for(unsigned int i=0;i<nx;++i)
    (*this)(i,i)=1.;
}

static bool same_size(const Mat& m1, const Mat& m2)
{
  if ((m1.SizeX() == m2.SizeX()) && (m1.SizeX() == m2.SizeY())) return true;
  cout << " matrices have different sizes" << endl;
  abort(); // we should in fact throw an exception.
  return false;
}

static bool same_size(const Vect& v1, const Vect& v2)
{
  if (v1.Size() == v2.Size()) return true;
  cout << " vectors have different sizes" << endl;
  abort(); // we should in fact throw an exception.
  return false;
}

Mat Mat::operator +(const Mat& Right) const
{
  same_size(*this,Right);
  Mat res = (*this);
  res += Right;
  return res;
}

Mat Mat::operator -(const Mat& Right) const
{
  same_size(*this,Right);
  Mat res = (*this);
  res -= Right;
  return res;
}

void Mat::operator +=(const Mat& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  unsigned int size = nx*ny;
  for(unsigned int i=0;i<size;++i, ++a, ++b)
    *a += *b;
}

void Mat::operator -=(const Mat& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  unsigned int size = nx*ny;
  for(unsigned int i=0;i<size;++i, ++a, ++b)
    *a -= *b;
}

Mat Mat::operator *(const double Right) const 
{
  Mat res = (*this);
  res *= Right;
  return res;
}

Mat operator *(const double Left, const Mat &Right)
{
  Mat res = Right;
  res *= Left;
  return res;
}
 
void Mat::operator *=(const double Right)
{
  unsigned int size = nx*ny;
  double *a = data;
  for(unsigned int i=0;i<size;++i, ++a)
    *a *= Right;
}

/* OLD VERY SLOW
Mat Mat::operator *(const Mat& Right) const
{
  if(nx != Right.SizeY()) {
    cout << "Mat::operator * ERROR nx != Right.SizeY() (" 
	 << nx << ',' << Right.SizeY() << ')' << endl;
    abort();
  }
  Mat res(Right.SizeX(),ny);
  
  for(unsigned int j=0;j<res.SizeY();++j) {
    for(unsigned int i=0;i<res.SizeX();++i) {  
      for(unsigned int k=0;k<nx;++k) {
	res(i,j) += (*this)(k,j)*Right(i,k);
      }
    }
  }
  return res;
}
*/


extern "C" {
void dgemm_(char &transa,
	    char &transb,
	    int &m, // number of rows of op(A)
	    int &n, // number of cols of op(B)
	    int &k, // number of columns of op(a) and numb of rows of op(B)
	    double &alpha,
	    const double *A, // const invented from the doc...
	    int &lda, // leading dim of A
	    const double *B, // const invented form the doc...
	    int &ldb,
	    double &beta,
	    double *C,
	    int &ldc);
	    
  // returns C = alpha A*B + beta C
  
  // note : if A(x,y)  in blas parlance :
  // number of rows = nx, number of columns = ny
  // which means that the "printout" reads (!= ours)
  // for (int i=0; i<nx; ++i)
  //   {
  //   for (int j=0;j<ny;++j)
  //       cout << A(i,j) << ' ';
  //   cout << endl;
  //   }

}


Mat Mat::operator *(const Mat& Right) const
{
  if(nx != Right.ny) {
    cout << "Mat::operator *= ERROR nx != Right.ny" << endl;
    abort();
  }
  Mat ResultMat(Right.nx,ny);

#define USE_DGEMM
#ifdef USE_DGEMM
  /* The definition chosen originally conflicts with usual definitions
      for(unsigned int k=0;k<nx;++k) {
	res(i,j) += (*this)(k,j)*Right(i,k);

	we will get the expected behaviour from DGEMM by swapping A and B
  */

  char transa='N';
  char transb='N';
  int m = Right.nx; // 
  int n = ny; //
  int k = nx; // = Right.ny 
  double alpha = 1;
  int lda = Right.nx;
  int ldb = nx;
  double beta = 0;
  int ldc = ResultMat.nx;
  
  if (Right.nx == 0) // we should throw an exception
    {
      std::cerr << " in matrix multiplication : Right.nx == 0 cannot work " << std::endl;
      abort ();
    }
	

  dgemm_(transa,transb,m,n,k,
	 alpha, Right.data,lda, 
	 data, ldb, 
	 beta, ResultMat.data, ldc);

#else
  double *res = ResultMat.data;
  for(size_t j=0;j<ResultMat.ny;++j) {
    for(size_t i=0;i<ResultMat.nx;++i,++res) {  
      double *left  = & data[j*nx];
      double *right = & Right.data[i];
      for(size_t k=0;k<nx;++k,++left,right+=Right.nx) {
	*res += (*left)*(*right);
      }
    }
  }
#endif /*USE_DGEMM  */
  return ResultMat;
}

Vect Mat::operator *(const Vect& Right) const
{
  if(nx != Right.Size()) {
    cout << "Mat::operator *= ERROR nx != Right.Size()" << endl;
    abort();
  }
  Vect res(ny);
  for(unsigned int j=0;j<ny;++j)
    {
      double acc = 0;
      for(unsigned int i=0;i<nx;++i)
	acc += (*this)(i,j)*Right(i);
    res(j) = acc;
    }
  return res;
}

void Mat::operator *=(const Mat& Right)
{
  Mat res = (*this)*Right;
  (*this) = res;
}

Mat & Mat::operator =(const Mat& Right){
  allocate(Right.SizeX(),Right.SizeY());
  memcpy(data,Right.Data(),nx*ny*sizeof(double));
  return (*this);
}

Mat::operator double() const
{
  if(nx!=1 || ny !=1) {
    cout << "Mat::operator double() error, nx=ny=1 needed, you have nx=" 
	 << nx <<" ny=" << ny << endl;
    abort();
  }
  return (*this)(0,0);
}

Mat::operator Vect() const
{
  if(nx!=1) {
    cout << "Mat::operator Vect() error, nx=1 needed, you have nx=" 
	 << nx << endl;
    abort();
  }
  Vect res(ny);
  for(unsigned int i=0;i<ny;i++) {
    res(i) = (*this)(0,i);
  }
  return res;
}

Mat Mat::transposed() const {
  Mat res(ny,nx);
  
  for(unsigned int i=0;i<nx;i++) {
    for(unsigned int j=0;j<ny;j++) {
      res(j,i)=(*this)(i,j);
    }
  }
  return res;
}


void Mat::ExtractSubMat(const unsigned IndexMapping[], const unsigned NewSize,
			Mat &SubMat) const
{
  SubMat.allocate(NewSize,NewSize);
  for (unsigned j=0; j< NewSize; ++j)
    {
      unsigned oldj = IndexMapping[j];
      for (unsigned i=0; i< NewSize; ++i)
	SubMat(i,j) = (*this)(IndexMapping[i],oldj);
    }
}

Mat Mat::SubBlock
(unsigned int x_min,unsigned int x_max,unsigned int y_min,unsigned int y_max) const {
  if( x_min<0 || x_max >= nx || y_min <0 || y_max>= ny ) {
    cout << "Mat::SubBlockFromIndexes ERROR, trying to get a sub-matrix with indices" << endl;
    cout << "x_min,x_max,y_min,y_max = "
	 << x_min << ","
	 << x_max << ","
    	 << y_min << ","
	 << y_max << endl;
    cout << "nx,ny = "<< nx << "," << ny << endl;
    abort();
  }
  unsigned int nx_new = (x_max-x_min+1);
  unsigned int ny_new = (y_max-y_min+1);
  Mat res(nx_new,ny_new);
  for(unsigned int j=0;j<ny_new;++j)
    for(unsigned int i=0;i<nx_new;++i)
      res(i,j) = (*this)(i+x_min,j+y_min);
  return res;
}

Mat Mat::WithoutRows(unsigned int y_min,unsigned int y_max) const {
  if( y_min <0 || y_max < y_min || y_max >= ny ) {
    cout << "Mat::WithoutRows ERROR " << endl;
    cout <<  y_min << " " <<  y_max << " " <<  ny << endl;
    abort();
  } 
  unsigned int nrows = y_max-y_min+1;
  Mat res(nx,ny-nrows);
  for(unsigned int j = 0 ;j<y_min ;j++)
    for(unsigned int i=0;i<nx;i++)
      res(i,j)=(*this)(i,j);
  for(unsigned int j = y_max+1 ;j<ny ;j++)
    for(unsigned int i=0;i<nx;i++)
      res(i,j-nrows)=(*this)(i,j);
  return res;
}

Mat Mat::WithoutColumns(unsigned int x_min,unsigned int x_max) const {
  if( x_min <0 || x_max < x_min || x_max >= nx ) {
    cout << "Mat::WithoutColumns ERROR " << endl;
    abort();
  } 
  unsigned int ncols = x_max-x_min+1;
  Mat res(nx-ncols,ny);
  for(unsigned int i=0;i<x_min;i++)
    for(unsigned int j = 0 ;j<ny ;j++)
      res(i,j)=(*this)(i,j);
  for(unsigned int i=x_max+1;i<nx;i++)
    for(unsigned int j = 0 ;j<ny ;j++)
      res(i-ncols,j)=(*this)(i,j);
  return res;
}





int Mat::readFits(const string &FitsName) {
  int status = 0;
  fitsfile *fptr = 0;
  fits_open_file(&fptr,FitsName.c_str(),0,&status);
  if (status)
   {
     cerr << " when opening file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  // first get the size of the image/matrix
  status=0;
  char value[256];
  fits_read_key(fptr, TSTRING, "NAXIS1", &value, NULL, &status);
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
      fits_report_error(stderr, status);
      return status;
    }
  int n1 = atoi(value);
  fits_read_key(fptr, TSTRING, "NAXIS2", &value, NULL, &status);
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
      fits_report_error(stderr, status);
      return status;
    }
  int n2 = atoi(value);
  if(n1<=0 || n2<=0) {
    cout << "Mat::readFits error NAXIS1,NAXIS2 = " << n1 << "," << n2 << endl;
    return -1;
  }
  allocate(n1,n2);

  status = 0;
  float nullval = 0;
  int anynull;
  fits_read_img(fptr, TDOUBLE, 1, nx*ny, &nullval,  data, &anynull, &status);
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
    }

  status = 0;
  fits_close_file(fptr, &status);
  if (status)
    {
     cerr << " when closing file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
   }
  return status;
}


int Mat::writeFits(const string &FitsName) const {
  // we cannot use fitsimage cause we want to write it in double precision
  
  int status = 0;
  fitsfile *fptr = 0;

  remove(FitsName.c_str());

  // enum FitsFileMode {RO = 0, RW = 1};
  fits_create_file(&fptr,FitsName.c_str(),  &status);
  if (status)
   {
     cerr << " when creating file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  // set a minimal header
  status = 0;
  long naxes[2];
  naxes[0]=nx;
  naxes[1]=ny;
  fits_write_imghdr(fptr,-64,2,naxes,&status);
  if (status)
   {
     cerr << " when writing minial header  ";
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }

  // say to cfitsio to take into account this new BITPIX
  status = 0;
  fits_flush_file(fptr,&status);
  if (status)
   {
     cerr << " when flushing  ";
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  status = 0;
  fits_write_img(fptr, TDOUBLE, 1, nx*ny, data, &status);
  if (status)
   {
     cerr << " when writing data ";
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  status = 0;
  fits_close_file(fptr, &status);
  if (status)
    {
     cerr << " when closing file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
   }
  return status;
}

void Mat::Symmetrize(const char* UorL) {
  
  if(nx!=ny) {
    cout << "Mat::Symmetrize ERROR nx!=ny nx,ny = " << nx << "," << ny << endl;
    abort();
  }
  
  
  for(unsigned int j=0;j<ny;j++)
    for(unsigned int i=j+1;i<nx;i++)
      if(UorL[0]=='L') { // x >= y
	(*this)(j,i)=(*this)(i,j);
      }else{
	(*this)(i,j)=(*this)(j,i);
      }
} 


int Mat::CholeskyInvert(const char *UorL) {
  if(nx!=ny) {
    cout << "Mat::CholeskyInvert ERROR nx!=ny nx,ny = " << nx << "," << ny << endl;
    abort();
  }

  char uorl[6]; strncpy(uorl,UorL,6);
  double *a = NonConstData();
  int n = SizeX();
  int info = 0;

  dpotrf_(uorl, &n, a,  &n, &info);
  if (info != 0) 
    {
      cout << "Mat::CholeskyInvert : could not factorize, info = " << info << endl;
      return info;
    }
  
  dpotri_(uorl, &n, a, &n, &info);
  if (info != 0) 
    {
      cout << "Mat::CholeskyInvert : could not invert, info = " << info << endl;
      return info;
    }
  return info;
}




//=================================================================

Vect::Vect(const unsigned int N) {
  data = NULL; 
  n=0;
  if(N<0) {
    cout << "Vect::Vect ERROR N = " << N <<  endl;
  }else{
    allocate(N);
  }
}

Vect::Vect(const Vect& other) {
  data = NULL; 
  n=0;
  allocate(other.Size());
  memcpy(data,other.Data(),sizeof(double)*n); 
}

void Vect::allocate(const unsigned int N) {
  if(N<0) {
    cout << "Vect::allocate ERROR N = " << N <<  endl;
  }
  if(n!=N) {
    n=N;
    if (data) 
      delete [] data;
    if(n>0)
      data = new double[n];
    else
      data = 0;
  }
  Zero();
};

int Vect::writeASCII(const std::string &FileName) const
{
  std::ofstream S(FileName.c_str());
  if (!S)
    {
      cout << " Vect::writeASCII() : cannot open " << FileName << endl;
      return 0;
    }
  int status = writeASCII(S);
  S.close();
  return status;
}



int Vect::writeASCII(ostream& Stream) const {
  Stream << n << endl;
  for(unsigned int i=0;i<n;i++)
    Stream << (*this)(i) << endl;
  return 0;
}

int Vect::readASCII(istream& Stream)  {
  unsigned int fn;
  Stream >> fn;
  allocate(fn);
  double val;
  for(unsigned int i=0;i<n;i++) {
    Stream >> val;
    (*this)(i) = val;
  }
  return 0;
}

void Vect::dump(ostream& Stream) const {
  for(unsigned int i=0;i<n;i++) {
    Stream << i << ":  " << float((*this)(i)) << endl;
  }
}

Vect Vect::operator *(const double Right) const 
{
  Vect res = (*this);
  res *= Right;
  return res;
}

Vect operator *(const double Left, const Vect &Right)
{
  Vect res = Right;
  res *= Left;
  return res;
}
 
void Vect::operator *=(const double Right)
{
  double *a = data;
  for(unsigned int i=0;i<n;++i, ++a)
    *a *= Right;
}


Vect Vect::operator +(const Vect& Right) const
{
  same_size(*this,Right);
  Vect res = (*this);
  res += Right;
  return res;
}

Vect Vect::operator -(const Vect& Right) const
{
  same_size(*this,Right);
  Vect res = (*this);
  res -= Right;
  return res;
}

void Vect::operator +=(const Vect& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  for(unsigned int i=0;i<n;++i, ++a, ++b)
    *a += *b;
}

void Vect::operator -=(const Vect& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  for(unsigned int i=0;i<n;++i, ++a, ++b)
    *a -= *b;
}

double Vect::operator *(const Vect& Right) const
{
  same_size(*this,Right);
  double res = 0;
  const double *a = data;
  const double *b = Right.Data();
  for(unsigned int i=0;i<n;++i, ++a, ++b)
    res += (*a)*(*b);
  return res;
}

Vect & Vect::operator =(const Vect& Right){
  allocate(Right.Size());
  memcpy(data,Right.Data(),n*sizeof(double));
  return (*this);
}

Mat Vect::transposed() const {
  Mat trans(n,1);
  for(unsigned int i=0;i<n;++i) {
    trans(i,0)=(*this)(i);
  }
  return trans;
}

#ifdef STORAGE /* causes troubles */
Vect::operator double() const
{
  if(n!=1) {
    cout << "Vect::operator double() error, n=1 needed, you have n=" 
	 << n << endl;
    abort();
  }
  return (*this)(0);
}
#endif


Vect::operator Mat() const {
//Mat Vect::asMat() const {
  Mat mat(1,n);
  for(unsigned int i=0;i<n;++i) {
    mat(0,i)=(*this)(i);
  }
  return mat;
}

}}}
