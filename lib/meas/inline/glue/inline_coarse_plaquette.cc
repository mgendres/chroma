// $Id: inline_plaquette.cc,v 3.15 2007-09-27 04:33:49 edwards Exp $
/*! \file
 *  \brief Inline plaquette
 */

#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/glue/inline_coarse_plaquette.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/io/named_objmap.h"

#include "meas/inline/io/default_gauge_field.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma 
{  
  //! Plaquette input
  void read(XMLReader& xml, const string& path, InlineCoarsePlaquetteEnv::Params::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);
    param.cgs = CreateGaugeStateEnv::nullXMLGroup();

    switch (version) 
    {
    case 2:
      if (paramtop.count("GaugeState") != 0)
	param.cgs = readXMLGroup(paramtop, "GaugeState", "Name");
      break;

    default:
      QDPIO::cerr << "InlineCoarsePlaquetteEnv::Params::Param_t: " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }
  }

  //! Plaquette output
  void write(XMLWriter& xml, const string& path, const InlineCoarsePlaquetteEnv::Params::Param_t& param)
  {
    push(xml, path);

    int version = 2;
    write(xml, "version", version);
    xml << param.cgs.xml;

    pop(xml);
  }


  //! Plaquette input
  void read(XMLReader& xml, const string& path, InlineCoarsePlaquetteEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! Plaquette output
  void write(XMLWriter& xml, const string& path, const InlineCoarsePlaquetteEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  namespace InlineCoarsePlaquetteEnv 
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

    const std::string name = "COARSE_PLAQUETTE";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateGaugeStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    // Params
    Params::Params() 
    { 
      frequency = 0; 
      param.cgs          = CreateGaugeStateEnv::nullXMLGroup();
      named_obj.gauge_id = InlineDefaultGaugeField::getId();
      xml_file ="";
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Params
	read(paramtop, "Param", param);

	// Ids
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) {
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      if( params.xml_file != "" ) 
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);
	push( xml_out, "CoarsePlaquette");
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

    void 
    InlineMeas::func(const unsigned long update_no, 
		     XMLWriter& xml_out) 
    {
      START_CODE();
    
      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;

      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

	// Set the construct state and modify the fields
	{
	  std::istringstream  xml_s(params.param.cgs.xml);
	  XMLReader  gaugetop(xml_s);

	  Handle<CreateGaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	    cgs(TheCreateGaugeStateFactory::Instance().createObject(params.param.cgs.id,
								    gaugetop,
								    params.param.cgs.path));

	  Handle<GaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	    state((*cgs)(u));

	  // Pull the u fields back out from the state since they might have been
	  // munged with gaugeBC's
	  u = state->getLinks();
	}
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineCoarsePlaquetteEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineCoarsePlaquetteEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }

      push(xml_out, "CoarsePlaquette");
      write(xml_out, "update_no", update_no);

      multi1d<Double> w_plaq(Nd+1);
      multi1d<Double> s_plaq(Nd+1); 
      multi1d<Double> t_plaq(Nd+1); 
      multi1d<Double> link(Nd+1); 
      multi1d< multi2d<Double> > plane_plaq(Nd+1);

      for (int mu=0; mu<Nd+1; ++mu) {
        plane_plaq[mu].resize(Nd,Nd);
      }

      MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);

      for (int mu=0; mu<Nd+1; ++mu) {
        write(xml_out, "w_plaq_"+std::to_string(mu*1LL), w_plaq[mu]);
        write(xml_out, "s_plaq_"+std::to_string(mu*1LL), s_plaq[mu]);
        write(xml_out, "t_plaq_"+std::to_string(mu*1LL), t_plaq[mu]);
  
        if (Nd >= 2)
        {
          write(xml_out, "plane_01_plaq_"+std::to_string(mu*1LL), plane_plaq[mu][0][1]);
        }
  
        if (Nd >= 3)
        {
          write(xml_out, "plane_02_plaq_"+std::to_string(mu*1LL), plane_plaq[mu][0][2]);
          write(xml_out, "plane_12_plaq_"+std::to_string(mu*1LL), plane_plaq[mu][1][2]);
        }
  
        if (Nd >= 4)
        {
          write(xml_out, "plane_03_plaq_"+std::to_string(mu*1LL), plane_plaq[mu][0][3]);
          write(xml_out, "plane_13_plaq_"+std::to_string(mu*1LL), plane_plaq[mu][1][3]);
          write(xml_out, "plane_23_plaq_"+std::to_string(mu*1LL), plane_plaq[mu][2][3]);
        }
  
        write(xml_out, "link_"+std::to_string(mu*1LL), link[mu]);
      } 
    
      pop(xml_out); // pop("CoarsePlaquette");
    
      END_CODE();
    } 

  } // namespace InlineCoarsePlaquetteEnv

} // namespace Chroma
