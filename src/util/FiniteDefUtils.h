
/** @file FiniteDefUtils.h
 *  @brief Implements finite deformation operations.
 *  
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2023], 
 */

#ifndef FINITEDEF_UTILS_H
#define FINITEDEF_UTILS_H

#include <jem/base/System.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/utilities.h>

#include "Arrays.h"
#include "TensorUtils.h"

using jem::String;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;

//-----------------------------------------------------------------------
//   typedefs
//-----------------------------------------------------------------------

// A pointer to a function that computes the spatial derivatives of
// the interpolation matrix. This is the so-called B-matrix.
// It points to the corresponding function for 1D, 2D and 3D case for
// a Total Lagrangian formulation.

typedef void        (*ShapeGradsTLFunc)

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Functions to get modified gradient matrix B0 = B * F
//-----------------------------------------------------------------------

// These functions compute the B-matrix given the matrix of shape function 
// 1st gradients (in global coordinates). For 1D, 2D (plane stress/strain)
// and 3D problems.

void                  get1DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

void                  get2DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

void                  get3DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions for a Total 
// Lagrangian formulation.

ShapeGradsTLFunc        getShapeGradsTLFunc

  ( int                 rank );

//-----------------------------------------------------------------------
// Compute deformation gradient f 
// from nodal displacements u and shape function gradients g
//-----------------------------------------------------------------------

void   evalDeformationGradient

  (       Matrix&    f,
    const Vector&    u,
    const Matrix&    g );

//-----------------------------------------------------------------------
// Compute Green-LagrangeStrain vector from deformation gradient
//-----------------------------------------------------------------------

void   getGreenLagrangeStrain
  
  (       Vector&    eps,
    const Matrix&    f );

#endif