// $Id: readmilc_d.h,v 3.3 2009-05-14 20:11:31 bjoo Exp $

/*! \file
 *  \brief Read a MILC gauge configuration written in the 1997 format
 */

#ifndef __readmilc_d_h__
#define __readmilc_d_h__

#include "io/milc_io.h"

namespace Chroma {

//! Read a MILC gauge configuration written in the 1997 format
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC_d(XMLReader& xml, multi1d<LatticeColorMatrixF>& u, const string& cfg_file);

//! Read a MILC gauge configuration written in the 1997 format
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC_d(XMLReader& xml, multi1d<LatticeColorMatrixD>& u, const string& cfg_file);

//! Read a MILC gauge configuration written in the 1997 format
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC_d(MILCGauge_t& header, multi1d<LatticeColorMatrixF>& u, const string& cfg_file);


}  // end namespace Chroma

#endif
