// -*- C++ -*-
// $Id: unprec_one_flavor_ratio_rat_rat_monomial5d_w.h,v 3.1 2008-05-23 21:31:35 edwards Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 5D ferm monomials
 */

#ifndef __unprec_one_flavor_ratio_rat_rat_monomial5d_w_h__
#define __unprec_one_flavor_ratio_rat_rat_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_ratio_rat_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/one_flavor_ratio_rat_rat_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Wrapper class for 5D 2-flavor unprec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecOneFlavorWilsonTypeFermRatioRatRatMonomial5D :
    public  OneFlavorRatioRatRatExactUnprecWilsonTypeFermMonomial5D< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
  {
  public: 
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // Construct out of a parameter struct. Check against the desired FermAct name
    UnprecOneFlavorWilsonTypeFermRatioRatRatMonomial5D(const OneFlavorWilsonTypeFermRatioRatRatMonomialParams& param_);

  protected:

    //! Accessor for pseudofermion (read only)
    const multi1d< multi1d<T> >& getPhi(void) const {return phi;}

    //! mutator for pseudofermion
    multi1d< multi1d<T> >& getPhi(void) {return phi;}

    //! Get at fermion action
    const WilsonTypeFermAct5D<T,P,Q>& getNumerFermAct(void) const { 
      return *fermact_num;
    }

    //! Get at fermion action
    const WilsonTypeFermAct5D<T,P,Q>& getDenomFermAct(void) const { 
      return *fermact_den;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getNumerActionInvParams(void) const { 
      return actionInvParam_num;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getNumerForceInvParams(void) const { 
      return forceInvParam_num;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getDenomActionInvParams(void) const { 
      return actionInvParam_den;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getDenomForceInvParams(void) const { 
      return forceInvParam_den;
    }

    //! Return number of roots in used
    int getNPF() const {return num_pf;}

    //! Return the partial fraction expansion for the force calc
    const RemezCoeff_t& getNumerFPFE() const {return fpfe_num;}

    //! Return the partial fraction expansion for the action calc
    const RemezCoeff_t& getNumerSPFE() const {return spfe_num;}

    //! Return the partial fraction expansion for the heat-bath
    const RemezCoeff_t& getNumerSIPFE() const {return sipfe_num;}

    //! Return the partial fraction expansion for the force calc
    const RemezCoeff_t& getDenomFPFE() const {return fpfe_den;}

    //! Return the partial fraction expansion for the action calc
    const RemezCoeff_t& getDenomSPFE() const {return spfe_den;}

    //! Return the partial fraction expansion for the heat-bath
    const RemezCoeff_t& getDenomSIPFE() const {return sipfe_den;}

  private:
    // Hide empty constructor and =
    UnprecOneFlavorWilsonTypeFermRatioRatRatMonomial5D();
    void operator=(const UnprecOneFlavorWilsonTypeFermRatioRatRatMonomial5D&);

    // Pseudofermion field phi
    multi1d< multi1d<T> > phi;

    // A handle for the UnprecWilsonFermAct
    Handle<const UnprecWilsonTypeFermAct5D<T,P,Q> > fermact_num;

    // A handle for the UnprecWilsonFermAct
    Handle<const UnprecWilsonTypeFermAct5D<T,P,Q> > fermact_den;

    // The parameters for the inversion
    GroupXML_t actionInvParam_num;
    GroupXML_t forceInvParam_num;

    // The parameters for the inversion
    GroupXML_t actionInvParam_den;
    GroupXML_t forceInvParam_den;

    // Number of nth-roots
    int num_pf;

    // Coefficients and roots of partial fractions
    RemezCoeff_t  fpfe_num;
    RemezCoeff_t  spfe_num;
    RemezCoeff_t  sipfe_num;

    // Coefficients and roots of partial fractions
    RemezCoeff_t  fpfe_den;
    RemezCoeff_t  spfe_den;
    RemezCoeff_t  sipfe_den;
  };


} //end namespace chroma



#endif
