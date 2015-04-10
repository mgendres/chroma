// -*- C++ -*-
/*! \file
 * \brief Determine if config satisfies admiss condition of Luscher
 *
 *  ---More documentation here--- 
 */

#ifndef __inline_admiss_h__
#define __inline_admiss_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineAdmissEnv 
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

        double s_min; // Perform a scan in s, starting at s_min
        double s_max; // Ending at s_max
        double n_steps; // number of steps to take between s_min and s_max

      } param;

      struct NamedObject_t
      {
	std::string     gauge_in;       /*!< Gauge fields */
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

  } // namespace InlineAdmissEnv


} // namespace Chroma

#endif
