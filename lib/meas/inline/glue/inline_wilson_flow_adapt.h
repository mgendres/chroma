// -*- C++ -*-
// inline_reweight_w.h
/*! \file
 * \brief Wilson Flow a configuration
 *
 *  ---More documentation here--- 
 */

#ifndef __inline_wilson_flow_adapt_h__
#define __inline_wilson_flow_adapt_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineWilsonFlowAdaptEnv 
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
	Real tol;
	Real wtime ;
	Real max_eps ;
	int t_dir; // the time direction of measurements 
      } param;

      struct NamedObject_t
      {
	std::string     eps_in;       /*!< step size */
	std::string     eps_out;       /*!< step size */

	std::string     gauge_in;       /*!< Gauge fields */
	std::string     gauge_out;       /*!< Gauge fields */
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline task for running the gauge Wilson Flow
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

  } // namespace InlineWilsonFlowAdaptEnv


} // namespace Chroma

#endif
