#ifndef __interpolate_h__
#define __interpolate_h__

namespace Chroma
{


  void CoolPlaquettes( multi1d<LatticeColorMatrix> & u);
  void CoolCubes( multi1d<LatticeColorMatrix> & u);
  void CoolHypercubes( multi1d<LatticeColorMatrix> & u);

  // file, lattice, nrow
  void DebugWrite(const std::string &,  const multi1d<LatticeColorMatrix> &, multi1d<int>&);

} // namespace Chroma

#endif
