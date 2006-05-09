// -*- C++ -*-
// $Id: inline_stoch_baryon_w.h,v 3.1 2006-05-09 20:32:24 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic baryon operator
 *
 * Form-factors
 */

#ifndef __inline_stoch_baryon_h__
#define __inline_stoch_baryon_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStochBaryonEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineStochBaryonParams 
  {
    InlineStochBaryonParams();
    InlineStochBaryonParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
      int              mom2_max;           /*!< (mom)^2 <= mom2_max */
    } param;

    PropSourceSmear_t  source_smearing;
    PropSinkSmear_t    sink_smearing;

    struct Prop_t
    {
      //! Operators
      struct Operator_t
      {
	multi1d<std::string> soln_files;
      };

      std::string          op_file;
      multi1d<Operator_t>  op;
    };

    struct NamedObject_t
    {
      Prop_t                prop;
      std::string           gauge_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of stochastic baryon operators
  /*! \ingroup inlinehadron */
  class InlineStochBaryon : public AbsInlineMeasurement 
  {
  public:
    ~InlineStochBaryon() {}
    InlineStochBaryon(const InlineStochBaryonParams& p) : params(p) {}
    InlineStochBaryon(const InlineStochBaryon& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineStochBaryonParams params;
  };

}

#endif
