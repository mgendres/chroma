#include "chromabase.h"
#include "interpolate.h"
#include "util/gauge/reunit.h"
#include "util/gauge/su3proj.h"

namespace Chroma 
{


  void GetLinkMask(multi1d<LatticeBoolean>& linkB, int p)
  {

    linkB.resize(Nd);

    //  links takes the values:
    //
    //          -------------
    //        /      /      /|
    //       /      /      / |
    //       --------------  |
    //     /      /|     /| /|
    //    /      / 2    / |/ |
    //   /--0-------0--/  |  |
    //   |      | 2    | /| /
    //   0      1/     0/ |/
    //   |--1---|--1---|  /
    //   0      1      0 /
    //   |      |      |/
    //   ---0------0---
    //
    // where
    // 0 : plaquette boundary links
    // 1 : cube boundary links
    // 2 : hypercube boundary links
    // 3 : hypercube bulk links

    // This holds integer values used to identify links
    multi1d<LatticeInt> link(Nd);
    for(int mu=0; mu<Nd; ++mu) { link[mu] = zero; }

    // LinkB is true for links of value p
    for(int mu=0; mu<Nd; ++mu) {
      for(int sig=0; sig<Nd; ++sig) {
        if (sig==mu) continue;
        link[mu] += Layout::latticeCoordinate(sig)%2;
      }
      linkB[mu] = (link[mu]==p);
    }

  }

  void GetPlaquetteMask(multi2d<LatticeBoolean>& plaquetteB, int p)
  {

    plaquetteB.resize(Nd,Nd);

    //  plaquettes take the values:
    //
    //          -------------
    //        /      /      /|
    //       /      /      / |
    //       --------------  |
    //     /      /|     /| /|
    //    /      /1|    / |/ |
    //   /-------------/  |  |
    //   |      | /    | /| /
    //   |  0   |/ 0   |/ |/
    //   |------|------|  /
    //   |      |      | /
    //   |  0   |  0   |/
    //   --------------
    //

    // This holds integer values used to identify plaquettes
    multi2d<LatticeInt> plaquette(Nd,Nd);
    for(int mu=0; mu<Nd; ++mu) {
      for(int nu=0; nu<Nd; ++nu) {
        plaquette[mu][nu] = zero;
      }
    }

    // plaquetteB is true for plaquettes of value p;
    // staples on these plaquettes will be included in the cooling
    for(int mu=0; mu<Nd; ++mu) {
      for(int nu=0; nu<Nd; ++nu) {
        for(int sig=0; sig<Nd; ++sig) {
          if ( (sig==mu)||(sig==nu) ) continue;
          plaquette[mu][nu] += Layout::latticeCoordinate(sig)%2;
        }
        if (mu==nu) {
          plaquetteB[mu][nu] = false; // this probably doesn't matter
        } else { 
          plaquetteB[mu][nu] = (plaquette[mu][nu]==p);
        }
      }
    }

  }

  void MaskedMesPlaq(Double& m_plaq, multi1d<LatticeColorMatrix> & u, int p)
  {

    multi2d<LatticeBoolean> plaquetteB(Nd,Nd);
    GetPlaquetteMask(plaquetteB, p);

    m_plaq = 0.0;
    multi2d<LatticeDouble> plaq(Nd,Nd);
    for(int mu=1; mu < Nd; ++mu) {
      for(int nu=0; nu < mu; ++nu) {
        plaq[mu][nu] = where(plaquetteB[mu][nu], real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]))), Real(0.0));
        m_plaq += sum(plaq[mu][nu]);
      }
    }

    // Normalization; lazy man's approach
    Double n_plaq(0.0);
    for(int mu=1; mu < Nd; ++mu) {
      for(int nu=0; nu < mu; ++nu) {
        n_plaq += sum(where(plaquetteB[mu][nu], LatticeDouble(1.0), LatticeDouble(0.0)));
      }
    }


    m_plaq /= Nc*n_plaq;

  }

  void CoolCellInteriorLinks( multi1d<LatticeColorMatrix> & u, int p, Double eps, Real BlkAccu, int BlkMax)
  {

    if (p<1) {
      QDP_error_exit("CoolCellInteriorLinks : p must be greater than zero!");
    }

    multi1d<LatticeBoolean> linkB(Nd);
    GetLinkMask(linkB, p);

    multi2d<LatticeBoolean> plaquetteB(Nd,Nd);
    GetPlaquetteMask(plaquetteB, p-1);

    // For each mu, sum staples in +-nu direction only
    LatticeColorMatrix tmp;
    multi2d<LatticeColorMatrix> staple(Nd,Nd);
    for(int mu=0; mu<Nd; ++mu ) {
      for(int nu = 0; nu<Nd; ++nu ) staple[mu][nu]= zero;
      for(int nu = 0; nu<Nd; ++nu) {
        if(nu == mu)  continue;
        staple[mu][nu] += u[nu] * shift(u[mu], FORWARD, nu)  * adj(shift(u[nu], FORWARD, mu));
        tmp = adj(u[nu]) * u[mu] * shift(u[nu], FORWARD, mu);
        staple[mu][nu] += shift(tmp, BACKWARD, nu );
      }
    }

    // Then do a masked add of all the staples perp to mu
    multi1d<LatticeColorMatrix> u_smear = u;
    for(int mu=0; mu<Nd; ++mu ) {
      for(int nu = 0; nu<Nd; ++nu) {
        if(nu == mu)  continue;
        u_smear[mu] += eps * where(plaquetteB[mu][nu], staple[mu][nu], LatticeColorMatrix(zero)); 
      }
    }

    // Then SU-project and Reunitarize
    SU3proj( u_smear, BlkAccu, BlkMax);

    // Finally do a masked replace of the old links with the new ones
    // provided links originate from sites with siteB=true
    for(int mu=0; mu<Nd; ++mu ) {
      u[mu] = where(linkB[mu], u_smear[mu], u[mu]);
    }

  }

  void SU3proj( multi1d<LatticeColorMatrix> & u, Real BlkAccu, int BlkMax)
  {

    LatticeColorMatrix u_unproj;
    for(int mu=0; mu<Nd; ++mu ) { // Project links in mu direction

      u_unproj = adj(u[mu]);
      reunit(u[mu]);
      Double old_tr = sum(real(trace(u[mu] * u_unproj))) / toDouble(Layout::vol()*Nc);
      Double new_tr;

      int n_smr = 0;
      bool wrswitch = false;                      /* Write out iterations? */
      Double conver = 1;

      while ( toBool(conver > BlkAccu)  &&  n_smr < BlkMax )
      {
        ++n_smr;

        // Loop over SU(2) subgroup index
        for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
          su3proj(u[mu], u_unproj, su2_index);

        /* Reunitarize */
       // reunit(u[mu]);
        {
          LatticeBoolean bad;
          int numbad;
          reunit(u[mu], bad, numbad, REUNITARIZE_LABEL);
          if ( numbad > 0 )
          {
            QDPIO::cout << "BLOCK: WARNING unitarity violation\n";
            QDPIO::cout << "   n_smr = " << n_smr << "\n";
            QDPIO::cout << "   mu= " << mu << "  numbad= " << numbad << endl;
          }
        }

        /* Calculate the trace */
        new_tr = sum(real(trace(u[mu] * u_unproj))) / toDouble(Layout::vol()*Nc);

        if( wrswitch )
          QDPIO::cout << " BLOCK: " << n_smr << " old_tr= " << old_tr << " new_tr= " << new_tr << "\n";

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
