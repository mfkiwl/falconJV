
/** @file FiniteDefUtils.cpp
 *  @brief Implements finite deformation operations.
 *  
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2023], 
 */

#include <jem/base/array/tensor.h>

#include "BasicUtils.h"
#include "FiniteDefUtils.h"

using jem::System;
using jem::idx_t;
using jem::ALL;
using jem::END;
using jem::TensorIndex;
using jem::numeric::norm2;
using jem::io::endl;

//-----------------------------------------------------------------------
//   get1DShapeGradsTL
//-----------------------------------------------------------------------

void              get1DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f )

{
  JEM_PRECHECK ( b.size(0) == 1 &&
                 g.size(0) == 1 &&
                 f.size(0) == 1 &&
                 b.size(1) == g.size(1) );

  b = g * f;
}

//-----------------------------------------------------------------------
//   get2DShapeGradsTL
//-----------------------------------------------------------------------

void              get2DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f )

{
  JEM_PRECHECK ( b.size(0) == 4 &&
                 g.size(0) == 2 &&
                 f.size(0) == 2 &&
                 b.size(1) == 2 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    const int i0 = 2 * inode;
    const int i1 = i0 + 1;

    b(0,i0) = f(0,0) * g(0,inode);
    b(0,i1) = f(1,0) * g(0,inode);
    b(1,i0) = f(0,1) * g(1,inode);
    b(1,i1) = f(1,1) * g(1,inode);

    b(3,i0) = f(0,0) * g(1,inode) + f(0,1) * g(0,inode);
    b(3,i1) = f(1,1) * g(0,inode) + f(1,0) * g(1,inode);
  }
}

//-----------------------------------------------------------------------
//   get3DShapeGradsTL
//-----------------------------------------------------------------------

// The strain-displacement matrix B is for a strain vector stored
// as [epsilon_xx, epsilon_yy, epsilon_zz, epsilon_xy, epsilon_yz, epsilon_zx].

void              get3DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f )

{
  JEM_PRECHECK ( b.size(0) == 6 &&
                 g.size(0) == 3 &&
                 f.size(0) == 3 &&
                 b.size(1) == 3 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {

    const int i0 = 3 * inode;
    const int i1 = i0 + 1;
    const int i2 = i0 + 2;

    b(0,i0) = f(0,0) * g(0,inode);
    b(0,i1) = f(1,0) * g(0,inode);
    b(0,i2) = f(2,0) * g(0,inode);

    b(1,i0) = f(0,1) * g(1,inode);
    b(1,i1) = f(1,1) * g(1,inode);
    b(1,i2) = f(2,1) * g(1,inode);
    
    b(2,i0) = f(0,2) * g(2,inode);
    b(2,i1) = f(1,2) * g(2,inode);
    b(2,i2) = f(2,2) * g(2,inode);

    b(3,i0) = f(0,0) * g(1,inode) + f(0,1) * g(0,inode);
    b(3,i1) = f(1,0) * g(1,inode) + f(1,1) * g(0,inode);
    b(3,i2) = f(2,0) * g(1,inode) + f(2,1) * g(0,inode);

    b(4,i0) = f(0,1) * g(2,inode) + f(0,2) * g(1,inode);
    b(4,i1) = f(1,1) * g(2,inode) + f(1,2) * g(1,inode);
    b(4,i2) = f(2,1) * g(2,inode) + f(2,2) * g(1,inode);

    b(5,i0) = f(0,2) * g(0,inode) + f(0,0) * g(2,inode);
    b(5,i1) = f(1,2) * g(0,inode) + f(1,0) * g(2,inode);
    b(5,i2) = f(2,2) * g(0,inode) + f(2,0) * g(2,inode);

  }
}

//-----------------------------------------------------------------------
//   getShapeGradsTLFunc
//-----------------------------------------------------------------------

ShapeGradsTLFunc getShapeGradsTLFunc ( int rank )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  if      ( rank == 1 )
  {
    return & get1DShapeGradsTL;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeGradsTL;
  }
  else
  {
    return & get3DShapeGradsTL;
  }
}

//-----------------------------------------------------------------------
//    evalDeformationGradient (2D and 3D)
//-----------------------------------------------------------------------

void  evalDeformationGradient

  (       Matrix&  f,
    const Vector&  u,
    const Matrix&  g )
{
  int rank = g.size(0); f = 0.0;

  for ( int i = 0; i < rank; ++i )
  {
    for ( int j = 0; j < rank; ++j )
    {
      f(i,j) = dot ( g(j,ALL) , u[slice(i,END,rank)] );
    }
    f(i,i) += 1.;
  }
}
  
//-----------------------------------------------------------------------
//    getGreenLagrangeStrain (2D and 3D)
//-----------------------------------------------------------------------

void  getGreenLagrangeStrain

  (       Vector&  eps,
    const Matrix&  f )

{
  int rank = f.size(0); eps = 0.0;

  // Compute Green-Lagrange strain tensor: 1/2 * ( F^T*F - I )

  TensorIndex i,j,k;
  Matrix      tens( rank, rank );

  tens(i,j) = 0.5 * ( dot( f(k,i), f(k,j), k ) - where(i==j,1.,0.) );

  // Convert to vector (Voigt notation)
  
  eps = tensorUtils::tensor2voigtStrain ( tens, STRAIN_COUNTS[rank] );
}