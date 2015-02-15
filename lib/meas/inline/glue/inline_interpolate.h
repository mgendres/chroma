// -*- C++ -*-
/*! \file
 * \brief Interpolating a configuration
 *
 *  ---More documentation here--- 
 */

#ifndef __inline_interpolate_h__
#define __inline_interpolate_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineInterpolateEnv 
  {

    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	int version ;
	Double tolP;
	Double tolC;
	Double tolH;
	Double epsP;
	Double epsC;
	Double epsH;
      } param;

      struct NamedObject_t
      {
	std::string     gauge_in;       /*!< Gauge fields */
	std::string     gauge_out;       /*!< Gauge fields */
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline task for running the interpolation
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // namespace InlineInterpolateEnv


} // namespace Chroma

#endif
