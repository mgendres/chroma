// -*- C++ -*-
// $Id: prec_dwf_fermact_array_sse_w.h,v 1.2 2004-09-08 02:48:25 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_fermact_array_sse_w_h__
#define __prec_dwf_fermact_array_sse_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace SSEEvenOddPrecDWFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Params for DWF
  struct SSEEvenOddPrecDWFermActArrayParams
  {
    SSEEvenOddPrecDWFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, SSEEvenOddPrecDWFermActArrayParams& param);
  void write(XMLReader& xml, const string& path, const SSEEvenOddPrecDWFermActArrayParams& param);


  //! 4D style even-odd preconditioned domain-wall fermion action
  /*! \ingroup fermact
   *
   * 4D style even-odd preconditioned domain-wall fermion action. 
   * Follows notes of Orginos (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  class SSEEvenOddPrecDWFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion>
  {
  public:
    //! General FermBC
    SSEEvenOddPrecDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
				 const Real& WilsonMass_, const Real& m_q_, int N5_) : 
      fbc(fbc_), WilsonMass(WilsonMass_), m_q(m_q_), N5(N5_) 
      {
	a5=1;
	QDPIO::cout << "Construct SSEEvenOddPrecDWFermActArray: WilsonMass = " << WilsonMass 
		    << "  m_q = " << m_q 
		    << "  N5 = " << N5 
		    << "  a5 = " << a5 
		    << endl;

	init();
      }

    //! Copy constructor
    SSEEvenOddPrecDWFermActArray(const SSEEvenOddPrecDWFermActArray& a) : 
      fbc(a.fbc), WilsonMass(a.WilsonMass), m_q(a.m_q), a5(a.a5), N5(a.N5) {}

    //! Assignment
    SSEEvenOddPrecDWFermActArray& operator=(const SSEEvenOddPrecDWFermActArray& a)
      {fbc=a.fbc; WilsonMass=a.WilsonMass; m_q=a.m_q; a5=a.a5; N5=a.N5; return *this;}

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real quark_mass() const {return m_q;}

    //! Produce a linear operator for this action
    const EvenOddPrecLinearOperator< multi1d<LatticeFermion> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

    //! Produce a linear operator for this action but with quark mass 1
    const LinearOperator< multi1d<LatticeFermion> >* linOpPV(Handle<const ConnectState> state) const;

    //! Optimized inverter - this is temporary
    /*! 
     * Solves  M.psi = chi
     *
     * Optimized inverter - will rename to qpropT
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     */
    void opt_qpropT(multi1d<LatticeFermion>& psi, 
		    Handle<const ConnectState> state, 
		    const multi1d<LatticeFermion>& chi, 
		    const InvertParam_t& invParam,
		    int& ncg_had) const;
  
    //! Destructor is automatic
    ~SSEEvenOddPrecDWFermActArray() {fini();}

  protected:
    //! Private internal initializer
    void init();

    //! Private internal destructor
    void fini();

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;
  };

}

using namespace Chroma;

#endif
