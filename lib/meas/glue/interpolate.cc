#include "chromabase.h"
#include "interpolate.h"

namespace Chroma 
{

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
