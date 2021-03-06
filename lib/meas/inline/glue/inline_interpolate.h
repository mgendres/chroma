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
        // P : plaquette
        // C : cube
        // H : hypercube

        // tol : stop cooling when relative change in average is below tol
	Double tolP;
	Double tolC;
	Double tolH;
        // eps : coefficent in front of staple at each stage of cool
	Double epsP;
	Double epsC;
	Double epsH;
        // These control the SU(2)/SU(3) projection
        Real BlkAccu;
        int BlkMax;
        // This writesout the field after each interpolation stage
        bool debug;
        std::string debug_file;

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
