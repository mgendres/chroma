// -*- C++ -*-
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_read_map_obj_disk_h__
#define __inline_read_map_obj_disk_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "util/ferm/map_obj/map_obj_disk.h"
namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineReadMapObjDiskEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineReadMapObjDiskParams 
  {
    unsigned int frequency;

    struct File {
      std::string   file_name;
    } file;

    struct NamedObject_t {
      std::string   object_id;
    } named_obj;

    InlineReadMapObjDiskParams(XMLReader& xml_in, const std::string& path);
    
  };

  void read(XMLReader& xml_in, const std::string& path, InlineReadMapObjDiskParams& p) ;

  void write(XMLWriter& xml_out, const std::string& path, const InlineReadMapObjDiskParams& p);


  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineReadMapObjDisk : public AbsInlineMeasurement 
  {
  public:
    ~InlineReadMapObjDisk() {}
    InlineReadMapObjDisk(const InlineReadMapObjDiskParams& p) : params(p) {}
    InlineReadMapObjDisk(const InlineReadMapObjDisk& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineReadMapObjDiskParams params;
  };

};

#endif
