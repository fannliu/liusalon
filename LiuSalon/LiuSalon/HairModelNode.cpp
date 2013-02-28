#include "HairModelNode.h"
#include <fstream>

#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg;				\
		return MS::kFailure;		\
	}

#define MNoVersionString
#define MNoPluginEntry
#include <maya/MFnPlugin.h>


#define MAKE_INPUT(attr) \
CHECK_MSTATUS(attr.setKeyable(true)); \
CHECK_MSTATUS(attr.setStorable(true)); \
CHECK_MSTATUS(attr.setReadable(true)); \
CHECK_MSTATUS(attr.setWritable(true));
#define MAKE_OUTPUT(attr) \
CHECK_MSTATUS(attr.setKeyable(false)); \
CHECK_MSTATUS(attr.setStorable(false)); \
CHECK_MSTATUS(attr.setReadable(true)); \
CHECK_MSTATUS(attr.setWritable(false));
#define MAKE_ADDR(attr) \
CHECK_MSTATUS(attr.setKeyable(false)); \
CHECK_MSTATUS(attr.setStorable(false)); \
CHECK_MSTATUS(attr.setReadable(true)); \
CHECK_MSTATUS(attr.setWritable(false)); \
CHECK_MSTATUS(attr.setHidden(true));

// TODO: add necessary MObjects
MTypeId HairModelNode::id( 0x80000 );

void* HairModelNode::creator()
{
	return new HairModelNode;
}

MStatus HairModelNode::initialize()
{
	// TODO: fill in attributes
	return MS::kSuccess;
}

MStatus HairModelNode::compute(const MPlug& plug, MDataBlock& data)

{
	// call this function whenever attribute changes
	MStatus returnStatus;
	// TODO: fill in
	return MS::kSuccess;
}
