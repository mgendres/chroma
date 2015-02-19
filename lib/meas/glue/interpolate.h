#ifndef __interpolate_h__
#define __interpolate_h__

namespace Chroma
{

  // This cools the inner links of a p-cell
  // p = 2 : plaquette
  // p = 3 : cube
  // p = 4 : hypercube
  // Cooling is done by adding staple to link with coeef eps
  // then SU-projecting
  // BlkAccu and BlkMax are parameters cotrolling the SU-projection
  void CoolCellInteriorLinks( multi1d<LatticeColorMatrix> & u, int p, Double eps, Real BlkAccu, int BlkMax);

  // SU(3) projection
  void SU3proj( multi1d<LatticeColorMatrix> & u, Real BlkAccu, int BlkMax);

  // file, lattice, nrow
  void DebugWrite(const std::string &,  const multi1d<LatticeColorMatrix> &, multi1d<int>&);

} // namespace Chroma

#endif
