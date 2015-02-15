#include "chromabase.h"
#include "interpolate.h"
#include "util/gauge/reunit.h"
#include "util/gauge/su3proj.h"

namespace Chroma 
{

  void CoolInnerLinks( multi1d<LatticeColorMatrix> & u, int p, Double eps)
  {
    // p = 1 is plaquette inner links
    // p = 2 is cube inner links
    // p = 3 is hypercube inner links

    // These are the links to consider
    multi1d<LatticeInt> link(Nd);
    multi1d<LatticeBoolean> linkB(Nd);
    for(int mu=0; mu<Nd; ++mu) { link[mu] = 0; }

    // These are the plaquettes to consider
    multi2d<LatticeInt> plaquette(Nd,Nd);
    multi2d<LatticeBoolean> plaquetteB(Nd,Nd);
    for(int mu=0; mu<Nd; ++mu) {
      for(int nu=0; nu<Nd; ++nu) {
        plaquette[mu][nu] = 0;
      }
    }

    // Considered links are true
    for(int mu=0; mu<Nd; ++mu) {
      for(int sig=0; sig<Nd; ++sig) {
        if (sig!=mu) link[mu] += (Layout::latticeCoordinate(sig)%2);
      }
      linkB[mu] = (link[mu]==p);
    }

    // Considered plaquettes true
    for(int mu=0; mu<Nd; ++mu) {
      for(int nu=0; nu<Nd; ++nu) {
        for(int sig=0; sig<Nd; ++sig) {
          if ( (sig!=mu)&&(sig!=nu) ) plaquette[mu][nu] += (Layout::latticeCoordinate(sig)%2);
        }
        plaquetteB[mu][nu] = (plaquette[mu][nu]==p-1);
      }
    }

    // For each mu, sum staples in +-nu direction only
    LatticeColorMatrix tmp;
    multi2d<LatticeColorMatrix> staple(Nd,Nd);
    for(int mu=0; mu<Nd; ++mu ) {
      for(int nu=0; nu<Nd; ++nu ) staple[mu][nu]= zero;
      for(int nu = 0; nu<Nd; ++nu) {
        if(nu == mu)  continue;
        staple[mu][nu] += u[nu] * shift(u[mu], FORWARD, nu)  * adj(shift(u[nu], FORWARD, mu));
        tmp = adj(u[nu]) * u[mu] * shift(u[nu], FORWARD, mu);
        staple[mu][nu] += shift(tmp, BACKWARD, nu );
      }
    }

    // Then add up nu component of staples; only the unmasked ones contribute
    multi1d<LatticeColorMatrix> staple_sum(Nd);
    for(int mu=0; mu<Nd; ++mu ) {
      staple_sum[mu] = zero;
      for(int nu = 0; nu<Nd; ++nu) {
        if(nu == mu)  continue;
         staple_sum[mu] += where(plaquetteB[mu][nu], staple[mu][nu], LatticeColorMatrix(zero)); 
      }
    }

    // Then add masked links to current lattice
    for(int mu=0; mu<Nd; ++mu ) {
      u[mu] += eps * where(linkB[mu], staple_sum[mu], LatticeColorMatrix(zero));
    }

    // Then SU-project and Reunitarize
    LatticeColorMatrix u_unproj;
    for(int mu=0; mu<Nd; ++mu ) {
      u_unproj = adj(u[mu]);
      Double old_tr = sum(real(trace(u[mu] * u_unproj))) / toDouble(Layout::vol()*Nc);
      Double new_tr;
 
      int n_smr = 0;
      bool wrswitch = false;                      /* Write out iterations? */
      Double conver = 1;

      // THESE SHOULD NOT BE HARD CODED?
      // What are ideal values?
      const Real BlkAccu(0.001);
      int BlkMax=1000; 

      while ( toBool(conver > BlkAccu)  &&  n_smr < BlkMax )
      {
        ++n_smr;
    
        // Loop over SU(2) subgroup index
        for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
          su3proj(u[mu], u_unproj, su2_index);
    
        /* Reunitarize */
        reunit(u[mu]);
    
        /* Calculate the trace */
        new_tr = sum(real(trace(u[mu] * u_unproj))) / toDouble(Layout::vol()*Nc);
    
        if( wrswitch )
          QDPIO::cout << " BLOCK: " << n_smr << " old_tr= " << old_tr << " new_tr= " << new_tr;
    
        /* Normalized convergence criterion: */
        conver = fabs((new_tr - old_tr) / old_tr);
        old_tr = new_tr;
      }


    }


  }

  void DebugWrite(const std::string& file, const multi1d<LatticeColorMatrix>& u, multi1d<int>& nrow)
  {

    multi1d<int> coord(Nd);
    ColorMatrix mat;
    Complex cmpnt;
    Real re, im;

    QDPIO::cout << "Writing gauge field to : " << file << endl ;

    ofstream gaugeWrite(file.c_str());
    gaugeWrite << "x y z t mu row col re im" << endl ;

    for(coord(0)=0;coord(0)<nrow(0);++coord(0))
    for(coord(1)=0;coord(1)<nrow(1);++coord(1))
    for(coord(2)=0;coord(2)<nrow(2);++coord(2))
    for(coord(3)=0;coord(3)<nrow(3);++coord(3)) {
      for(int mu=0; mu<Nd; mu++) {
        mat = peekSite( u(mu), coord );
        for(int row=0; row<Nc; ++row)
        for(int col=0; col<Nc; ++col) {
          cmpnt=peekColor( mat, row, col );
          re = real(cmpnt);
          im = imag(cmpnt);
          gaugeWrite << coord(0) << " " \
                     << coord(1) << " " \
                     << coord(2) << " " \
                     << coord(3) << " " \
                     << mu << " " \
                     << row << " " \
                     << col << " " \
                     << re << " " \
                     << im << endl;
        }
      }
    }
    gaugeWrite.close();

  }



}  // namespace Chroma
