#ifndef __interpolate_h__
#define __interpolate_h__

namespace Chroma
{

  void CoolInnerLinks( multi1d<LatticeColorMatrix> & u, int p, Double eps);

  // file, lattice, nrow
  void DebugWrite(const std::string &,  const multi1d<LatticeColorMatrix> &, multi1d<int>&);

} // namespace Chroma

#endif
