#ifndef __inline_refine_h__
#define __inline_refine_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  // A namespace for this particular   measurement
  namespace InlineRefineEnv 
  {
    extern const std::string name;
    bool registerAll();
  

    //! Parameter structure
    struct Params 
    {
      // Default constructor -- forward declaration
      Params();

      // Construct from XML -- forward declaration
      Params(XMLReader& xml_in, const std::string& path);

      // Write myself out
      void write(XMLWriter& xml_out, const std::string& path);

      // How often should I measure me in an HMC run
      unsigned long frequency;

      // I need the ID of the named object for the prop
      std::string gauge_id;
      std::string rf_gauge_id;

      std::string xml_file; // Support output to own private XML File
    }; // struct
  }; // namespace InlineMyMeasEnv

  class InlineRefine : public AbsInlineMeasurement 
  {
  public:
    // Constructor -- default -- do nothing
    ~InlineRefine() {}

    // Constructor -- from param struct -- copy param struct inside me
    InlineRefine(const InlineRefineEnv::Params& p) : params(p) {}

    // Constructor -- copy constructor -- copy param struct from argument
    InlineRefine(const InlineRefine& p) : params(p.params) {}

    // Boiler plate
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 


  private:
    void func(const unsigned long update_no,
	      XMLWriter& xml_out);

    InlineRefineEnv::Params params;
  };

};

#endif
