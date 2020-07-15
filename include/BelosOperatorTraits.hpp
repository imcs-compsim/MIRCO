//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef BELOS_OPERATOR_TRAITS_HPP
#define BELOS_OPERATOR_TRAITS_HPP

/// \file BelosOperatorTraits.hpp
/// \brief Class which defines basic traits for the operator type.
///
#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

namespace Belos {

  /// \class UndefinedOperatorTraits
  /// \brief Class used to require specialization of OperatorTraits.
  /// 
  /// This class is used by \c OperatorTraits to ensure that
  /// OperatorTraits cannot be used unless a specialization for the
  /// particular scalar, multivector, and operator types has been
  /// defined.
  template<class ScalarType, class MV, class OP>
  class UndefinedOperatorTraits
  {
  public:
    /// \brief Function that will not compile if instantiation is attempted.
    ///
    /// Any attempt to compile this function results in a compile-time
    /// error.  Such an error means that the template specialization
    /// of Belos::OperatorTraits class either does not exist for type
    /// <tt>OP</tt>, or is not complete.
    static inline void notDefined() { 
      OP::this_type_is_missing_a_specialization(); 
    };
  };
 
  /// \class OperatorTraits
  /// \brief Class which defines basic traits for the operator type.
  ///
  /// This class defines the interface that Belos' linear solvers use
  /// when interacting with operators (e.g., matrices and
  /// preconditioners).  OperatorTraits is a "traits class," which
  /// means that it uses compile-time polymorphism to interface
  /// between Belos' solvers and the particular implementation of
  /// matrices, preconditioners, and multivectors.
  ///
  /// A specialization of this traits class must exist for the
  /// <tt>ScalarType</tt>, <tt>MV</tt> and <tt>OP</tt> types.  If not,
  /// this class will produce a compile-time error.  (The \c
  /// UndefinedOperatorTraits class provides an intelligible error
  /// message for that case.)
  ///
  /// \tparam ScalarType The scalar type: that is, the data type for
  ///   which scalar multiplication by a (multi)vector is defined, and
  ///   the data type produced by a dot product of two vectors.
  ///
  /// \tparam MV The multivector type.  A "multivector" contains one
  ///   or more vectors.  All vectors in a multivector have the same
  ///   number of entries and belong to the same vector space.  (What
  ///   the latter phrase specifically means depends on the MV class.
  ///   Informally, it makes sense to add two vectors belonging to the
  ///   same vector space.)
  ///
  /// \tparam OP The operator type.  An operator behaves as a function
  ///   (not in the strict mathematical sense) that takes a
  ///   multivector X as input, and fills a multivector Y with the
  ///   result of applying the multivector to X.
  ///
  /// \ingroup belos_opvec_interfaces
  template <class ScalarType, class MV, class OP>
  class OperatorTraits 
  {
  public:

    /// \brief Apply Op to x, putting the result into y.
    ///
    /// If Op, x, and y are real-valued, then applying the conjugate
    /// transpose (trans = CONJTRANS) means the same thing as applying
    /// the transpose (trans = TRANS).  
    ///
    /// If Op does not support applying the transpose and you use
    /// trans != NOTRANS, or if there is some other error in applying
    /// the operator, this method throws a subclass of std::exception.
    static void 
    Apply (const OP& Op, 
	   const MV& x, 
	   MV& y, 
	   ETrans trans = NOTRANS)
    { 
      // This will result in a deliberate compile-time error, if a
      // specialization of OperatorTraits has not been defined for the
      // MV and OP types.
      UndefinedOperatorTraits<ScalarType, MV, OP>::notDefined(); 
    };

    /// \brief Whether this operator implements applying the transpose.
    ///
    /// The instance of OP which is the first argument of \c Apply()
    /// is not required to support applying the transpose (or
    /// Hermitian transpose, if applicable).  If it <i>does</i>
    /// support applying the transpose, this method should return
    /// true.  Otherwise, it should return false.
    ///
    /// If the operator is complex, "can apply its transpose" means
    /// that it can apply both its transpose and its Hermitian
    /// transpose.
    ///
    /// We provide a default implementation of this method that
    /// conservatively returns false.  If you want the specialization
    /// of OperatorTraits for OP to advertise that operators of type
    /// OP may implement applying the transpose, override the default
    /// implementation in the specialization.
    static bool HasApplyTranspose (const OP& Op) {
      return false;
    }
  };
  
} // end Belos namespace

#endif // BELOS_OPERATOR_TRAITS_HPP

// end of file BelosOperatorTraits.hpp
