/*! \file
 *  \brief Read a MILC gauge configuration written in the 1997 format and double the lattice
 */

#include "chromabase.h"
#include "io/milc_io.h"
#include "io/readmilc_d.h"
#include "qdp_util.h"    // from QDP

namespace Chroma {

//! Read a MILC configuration file
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC_d(MILCGauge_t& header, multi1d<LatticeColorMatrixF>& u, const string& cfg_file)
{
  START_CODE();

  u.resize(Nd);

  BinaryFileReader cfg_in(cfg_file); // for now, cfg_io_location not used

  int magic_number;
  read(cfg_in, magic_number);

  bool byterev = false ;
  if( magic_number != 20103){
    //need byte reversal
    byterev = true ;
    //QDP_error_exit("readMILC_d: unexpected byte order of file");
    QDPUtil::byte_swap((void *)&magic_number, sizeof(int), 1 );
  }
  if( magic_number != 20103){
    QDP_error_exit("readMILC_d: unexpected magic number");
  }
  

  // Check lattice size
  header.nrow.resize(Nd);
  read(cfg_in, header.nrow, Nd);
  if(byterev){
    QDPUtil::byte_swap((void *)&header.nrow[0], sizeof(int), Nd );
  }
  for(int j = 0; j < Nd; ++j)
    if ( 2*header.nrow[j] != Layout::lattSize()[j] )
      QDP_error_exit("readMILC_d: unexpected lattice size: header.nrow[%d]=%d",
                     j,header.nrow[j]);

  // Time stamp
  char date_tmp[65];
  cfg_in.readArray(date_tmp, 1, 64);
  date_tmp[64] = '\0';
  header.date = date_tmp;

  // Site order - only support non-sitelist format
  int order;
  read(cfg_in, order);
  if( order != 0)
    QDP_error_exit("readMILC_d: only support non-sitelist format");


  // Go ahead an read checksums, but will not use for now
  unsigned int sum29, sum31;
  read(cfg_in, sum29);
  read(cfg_in, sum31);
  if(byterev){
    QDPUtil::byte_swap((void *)&sum29,sizeof(int),1);
    QDPUtil::byte_swap((void *)&sum31,sizeof(int),1);
  }
  QDPIO::cout<<"Global sums (sum29, sum31): "<<sum29<<" "<<sum31<<endl; 

  /*
   * Read away...
   */

  for(int mu=0; mu < Nd; ++mu){
    u[mu]=1.0 ;
  }
  // MILC format has the directions inside the sites
  for(int site=0; site < Layout::vol(); ++site)
  {
    multi1d<int> coord = crtesn(site, Layout::lattSize()); // The coordinate
    // Read in Nd SU(3) matrices. 
    // NOTE: the su3_matrix layout should be the same as in QDP
    bool flag = true ;
    for(int mu=0; mu < Nd; ++mu){
      flag = flag&&(coord[mu]%2==0) ;
    }
    for(int mu=0; mu < Nd; ++mu){
      if(flag){
	read(cfg_in, u[mu], coord);    // read in a single site
	if(byterev)
	  QDPUtil::byte_swap((void *)&u[mu].elem(site).elem(),sizeof(Real),2*Nc*Nc);
      }//flag
    }//mu
  }//site

  cfg_in.close();
  
  END_CODE();
}



//! Read a MILC configuration file
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC_d(XMLReader& xml, multi1d<LatticeColorMatrixF>& u, const string& cfg_file)
{
  START_CODE();

  MILCGauge_t header;

  // Read the config and its binary header
  readMILC_d(header, u, cfg_file);

  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;
  write(xml_buf, "MILC", header);

  try 
  {
    xml.open(xml_buf);
  }
  catch(const string& e)
  { 
    QDP_error_exit("Error in readMILC_d: %s",e.c_str());
  }

  END_CODE();
}

//! Read a MILC configuration file
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC_d(XMLReader& xml, multi1d<LatticeColorMatrixD>& u, const string& cfg_file)
{
  START_CODE();

  // MILC configs only in single-prec, 
  multi1d<LatticeColorMatrixF> uu;
  readMILC_d(xml, uu, cfg_file);

  u.resize(uu.size());
  for(int mu=0; mu<uu.size(); ++mu)
    u[mu] = uu[mu];

  END_CODE();
}

}  // end namespace Chroma
