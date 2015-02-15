// 
/*! \file
 * \brief Interpolation
 *
 */

#include "inline_interpolate.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/glue/interpolate.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"


// The wilson flow code
#include "meas/glue/wilson_flow_w.h"


namespace Chroma 
{ 

  namespace InlineInterpolateEnv 
  {
    //! read input -- gauge fields
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_in", input.gauge_in);
      read(inputtop, "gauge_out", input.gauge_out);
    }

    //! write output -- gauge fields
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_in", input.gauge_in);
      write(xml, "gauge_out", input.gauge_out);

      pop(xml);
    }


    //! read input
    void read(XMLReader& xml, const string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "version", input.version);
      read(inputtop, "tolP", input.tolP);
      read(inputtop, "tolC", input.tolC);
      read(inputtop, "tolH", input.tolH);
      read(inputtop, "epsP", input.epsP);
      read(inputtop, "epsC", input.epsC);
      read(inputtop, "epsH", input.epsH);


    }

    //! write output
    void write(XMLWriter& xml, const string& path, const Params::Param_t& input)
    {
      push(xml, path);
    
      write(xml, "version", input.version);
      write(xml, "tolP", input.tolP);
      write(xml, "tolC", input.tolC);
      write(xml, "tolH", input.tolH);
      write(xml, "epsP", input.epsP);
      write(xml, "epsC", input.epsC);
      write(xml, "epsH", input.epsH);


      pop(xml);
    }


    //! read input
    void read(XMLReader& xml, const string& path, Params& input)
    {
      InlineInterpolateEnv::Params tmp(xml, path);
      input = tmp;
    }

    //! write output
    void write(XMLWriter& xml, const string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlineInterpolateEnv 


  namespace InlineInterpolateEnv 
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }
      
    const std::string name = "INTERPOLATE";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //--------------------------------------------------------------------------
    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for source construction
	read(paramtop, "Param", param);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "Interpolate");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    // Create the diluted source and apply Lanczos quarature 
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      push(xml_out, "Interpolate");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Interpolating a configuration" << endl;

      proginfo(xml_out);    // Print out basic program info

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Write out the input
      write(xml_out, "Input", params);

      // Test and grab a reference to the gauge field
      // -- we really need only two gauge fields --
      multi1d<LatticeColorMatrix> u ; 
     
      push(xml_out,"GaugeFieldInfo");
      XMLBufferWriter gauge_xml;
      try
	{
	  u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_in);
	  TheNamedObjMap::Instance().get(params.named_obj.gauge_in).getRecordXML(gauge_xml);
	}
      catch( std::bad_cast ) 
	{
	  QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	  QDP_abort(1);
	}
      catch (const string& e) 
	{
	  QDPIO::cerr << name << ": map call failed: " << e << endl;
	  QDP_abort(1);
	}
      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      pop(xml_out);
      

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "GaugeObservables", u);
      
///////// INTERPOLATION CODE ///////     

      multi1d<LatticeColorMatrix> interp_u = u ; 

      Double tol[Nd];
      tol[0] = 0; // this one isn't used
      tol[1] = params.param.tolP;
      tol[2] = params.param.tolC;
      tol[3] = params.param.tolH;

      Double eps[Nd];
      eps[0] = 0; // this one isn't used
      eps[1] = params.param.epsP;
      eps[2] = params.param.epsC;
      eps[3] = params.param.epsH;

      Double w_plaq, s_plaq, t_plaq, link; // only w_plaq is used
      multi2d<Double> plane_plaq;
      Double plaq, r;

      for (int p=1; p<Nd; ++p) {
        MesPlq(interp_u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
        do {
          plaq = w_plaq;
          CoolInnerLinks( interp_u, p, eps[p]);
//          wilson_flow(xml_out, interp_u, 1, 0.1 , -1); // this is just to test the stopping condition
          MesPlq(interp_u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
          r = w_plaq / plaq - 1.0;
          QDPIO::cout << "\tRelative plaq difference between sweeps for inner " \
                      <<  p << "-cell links : " <<  r << endl;
        } while ( toBool(r>tol[p]) );
      }

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "InterpolatedGaugeObservables", interp_u);


       // Now store the configuration to a memory object
      {
	XMLBufferWriter file_xml, record_xml;
	push(file_xml, "gauge");
	write(file_xml, "id", int(0));
	write(file_xml, "InterpolateParams", params.param);
	pop(file_xml);
	record_xml << gauge_xml;

	// Store the gauge field
	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_out);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_out) = interp_u;
	TheNamedObjMap::Instance().get(params.named_obj.gauge_out).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_out).setRecordXML(record_xml);
      }

      pop(xml_out);  // Interpolate
      
      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;
      
      QDPIO::cout << name << ": ran successfully" << endl;
      
      END_CODE();
    }

  }
  

} // namespace Chroma
