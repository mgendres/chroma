#include "chromabase.h"
#include "interpolate.h"

namespace Chroma 
{

  void CoolInnerLinks( multi1d<LatticeColorMatrix> & u, int p, Double eps)
  {
    // p = 1 is plaquette inner links
    // p = 2 is cube inner links
    // p = 3 is hypercube inner links

    // These are the links to consider
    LatticeInt link[Nd];
    LatticeBoolean linkB[Nd];
    for(int mu=0; mu<Nd; ++mu) { link[mu] = 0; }

    // These are the plaquettes to consider
    LatticeInt plaquette[Nd][Nd];
    LatticeBoolean plaquetteB[Nd][Nd];
    for(int mu=0; mu<Nd; ++mu) {
      for(int nu=0; nu<Nd; ++nu) {
        plaquette[mu][nu] = 0;
      }
    }

    // Considered links are true
    for(int mu=0; mu<Nd; ++mu) {
      for(int sig=0; sig<Nd; ++sig) {
        if (sig!=mu) link[mu] += Layout::latticeCoordinate(sig);
      }
      linkB[mu] = (link[mu]%2==p);
    }

    // Considered plaquettes true
    for(int mu=0; mu<Nd; ++mu) {
      for(int nu=0; nu<Nd; ++nu) {
        for(int sig=0; sig<Nd; ++sig) {
          if ( (sig!=mu)&&(sig!=nu) ) plaquette[mu][nu] += Layout::latticeCoordinate(sig);
        }
        plaquetteB[mu][nu] = (plaquette[mu][nu]%2==p);
      }
    }

    // Here, for each mu staples need to be computed given nu
    // Then plaquette mask is applied to the staples
    // Then add up nu component of staples; only the unmasked ones contribute
    // Then mask links
    // Then add masked links to current lattice
    // Then SU-project all
    // Finally reunitarize

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
