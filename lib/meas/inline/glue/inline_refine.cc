
#include "inline_refine.h"

#include "chroma.h"
#include "debug_write.h"

using namespace QDP;

namespace Chroma 
{ 
  // Fill in the blanks
  namespace InlineRefineEnv 
  { 

    // Function to register with a factory
    // This is boilerplate stuff
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path) 
    {
      // Create an Params from the XML
      // Use it to create an InlineRefine
      return new InlineRefine( Params(xml_in, path) );
    }

    // The name of my measurement for the XML file
    // Change this for each measurement
    const std::string name = "REFINE";

    // Register the measurement in the measurement factory
    namespace { 
      bool registered = false;
    }

    bool registerAll()
    {
      bool success = true;
      if (! registered) { 
         success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
         registered = true;
      }
      return success;
    }



  // Param stuff
  Params::Params() { frequency = 0; gauge_id=""; rf_gauge_id=""; xml_file="";}

  Params::Params(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Get either the gauge_id in the XML or the ID of the default
      // field of no explicit ID exists.
      read(paramtop, "NamedObject/gauge_id", gauge_id);
      read(paramtop, "NamedObject/rf_gauge_id", rf_gauge_id);

      if( paramtop.count("xml_file") != 0 ) { 
	read(paramtop, "xml_file", xml_file);
      }
      else { 
	xml_file == "";
      }
	
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }

    // Read your own parameters from Param/ here

  }


  void
  Params::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Write out the output target/source configuration info
    QDP::write(xml_out, "gauge_id", gauge_id);
    QDP::write(xml_out, "rf_gauge_id", rf_gauge_id);

    if( xml_file != "" ){ 
      QDP::write(xml_out, "xml_file", xml_file);
    }


    pop(xml_out);
    }
  }

  void 
  InlineRefine::operator()(long unsigned int update_no, XMLWriter& xml_out) 
  {

    // This bit merely supports providing an external xml_file 
    // for this measurement
    if ( params.xml_file == "" ) { 
      
      func( update_no, xml_out );
    }
    else { 

      // Hash an XML file name from the user specified string
      // and the update number
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      // IN the global output, make a note that the output went
      // to this separate XML file
      push(xml_out, "refine");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
  }

  void 
  InlineRefine::func(unsigned long int update_no, XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch measure_time;
    measure_time.reset();
    measure_time.start();

    // Test that the gauge configuration exists in the map.
    XMLBufferWriter gauge_xml;
    try
    {
      // Try and get at the gauge field if it doesn't exist an exception will be thrown.
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);
      TheNamedObjMap::Instance().get(params.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      // If the object exists and is of the wrong type we will land in this catch.
      QDPIO::cerr << InlineRefineEnv::name << ": caught dynamic cast error" 
                  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      // If the object is not in the map we will land in this catch
      QDPIO::cerr << InlineRefineEnv::name << ": map call failed: " << e 
                  << endl;
      QDP_abort(1);
    }

    // If we got here, then the gauge field is in 
    // the map. The XML will have been captured.
    // Let us bind the references to a local name so 
    // we don't have to type the long lookup string every time.
    //
    // Note const here means we can't change the field
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.gauge_id);

    QDPIO::cout << InlineRefineEnv::name <<": Beginning " << endl;

    // Boilerplate stuff to the output XML
    push(xml_out, "refine");
    write(xml_out, "update_no", update_no);

    // Write info about the program
    proginfo(xml_out);

    // Write out the input
    params.write(xml_out, "Input");

    //----------------------------------------------------------------------
    // Compute something here....

    QDPIO::cout << "Gauge ID : " << params.gauge_id << endl;
    QDPIO::cout << "RF Gauge ID : " << params.rf_gauge_id << endl;

    multi1d<int> size(Nd);
    size = Layout::lattSize();

    QDPIO::cout << "Lattice shape : ";
    for (int i=0; i<Nd; ++i) QDPIO::cout << size[i] << " ";
    QDPIO::cout << endl;

    multi1d<int> rf_size(size);
    rf_size *= 2;

    QDPIO::cout << "RF Lattice shape : ";
    for (int i=0; i<Nd; ++i) QDPIO::cout << rf_size[i] << " ";
    QDPIO::cout << endl;

    QDPIO::cout << "Setting RF Lattice layout" << endl;
    Layout::setLattSize(rf_size);
    Layout::create();

    multi1d<LatticeColorMatrix> rf_u(Nd);

    DebugWrite("debug_gauge.dat", u, size);
    DebugWrite("debug_rf_gauge.dat", rf_u, rf_size);
    

    //---------------------------------------------------------------------

    // End of  your code.
    pop(xml_out);


    measure_time.stop();
    QDPIO::cout << InlineRefineEnv::name << ": total time = "
		<< measure_time.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineRefineEnv::name << ": ran successfully" << endl;
    END_CODE();
  } 

};
