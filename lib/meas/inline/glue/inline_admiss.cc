// 
/*! \file
 * \brief Admiss
 *
 */

#include "inline_admiss.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"


namespace Chroma 
{ 

  namespace InlineAdmissEnv 
  {
    //! read input -- gauge fields
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_in", input.gauge_in);
    }

    //! write output -- gauge fields
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_in", input.gauge_in);

      pop(xml);
    }


    //! read input
    void read(XMLReader& xml, const string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "version", input.version);
      read(inputtop, "s_min", input.s_min);
      read(inputtop, "s_max", input.s_max);
      read(inputtop, "n_steps", input.n_steps);
    }

    //! write output
    void write(XMLWriter& xml, const string& path, const Params::Param_t& input)
    {
      push(xml, path);
    
      write(xml, "version", input.version);
      write(xml, "s_min", input.s_min);
      write(xml, "s_max", input.s_max);
      write(xml, "n_steps", input.n_steps);

      pop(xml);
    }


    //! read input
    void read(XMLReader& xml, const string& path, Params& input)
    {
      InlineAdmissEnv::Params tmp(xml, path);
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
  } // namespace InlineAdmissEnv 


  namespace InlineAdmissEnv 
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
      
    const std::string name = "ADMISS";

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

	push(xml_out, "Admiss");
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
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      push(xml_out, "Admiss");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Compute Admiss on a configuration" << endl;

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
      
///////// ADMISS CODE ///////     
      multi1d<int> nrow(4);
      nrow = Layout::lattSize();

      // Compute plaquettes
      multi2d<LatticeDouble> plaquette(Nd,Nd);
      plaquette = zero;
      for(int mu=1; mu < Nd; ++mu)
      for(int nu=0; nu < mu; ++nu) {
        plaquette[mu][nu] = Nc;
        plaquette[mu][nu] -= real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])));
        plaquette[mu][nu] /= Nc;
      }

      // Compute fraction of plaquettes over s, for a given s
      Double tot;
      int steps = params.param.n_steps;
      multi1d<Double> s_vec(steps); 
      multi1d<Double> frac_vec(steps); 
      double s = params.param.s_min;
      double s_inc = (params.param.s_max-params.param.s_min)/steps;
      for (int i = 0; i < steps; ++i) {
        s_vec[i] = s;
        tot = 0.0;
        for(int mu=1; mu < Nd; ++mu)
        for(int nu=0; nu < mu; ++nu) {
          tot += sum(where(plaquette[mu][nu]>s, LatticeDouble(1.0), LatticeDouble(0.0)));
        }
        frac_vec[i] = tot/Double( Layout::vol()*Nd*(Nd-1)/2 );
        s += s_inc;
      }

      // Write results
      push(xml_out, "admiss_results");
      write(xml_out,"s_vals", s_vec) ;
      write(xml_out,"frac_over", frac_vec) ;
      pop(xml_out);  // elem

      pop(xml_out);  // Admiss
      
      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;
      
      QDPIO::cout << name << ": ran successfully" << endl;
      
      END_CODE();
    }

  }
  

} // namespace Chroma
