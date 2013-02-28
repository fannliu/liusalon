#include <maya/MPxCommand.h>
#include <maya/MFnPlugin.h>
#include <maya/MIOStream.h>
#include <maya/MString.h>
#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MSimple.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MDGModifier.h>
#include <maya/MPlugArray.h>
#include <maya/MVector.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MStringArray.h>
#include <list>

#include "LiuSalonCmd.h"
//#include "LiuSalonNode.h"

MStatus initializePlugin( MObject obj )
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj, "LiuSalon", "1.0", "Any");

	// load external functions & ui
	MGlobal::executeCommand("source \"" + plugin.loadPath() + "/functions.mel\"");
    MGlobal::executeCommand("source \"" + plugin.loadPath() + "/ui.mel\"");
	status = plugin.registerUI("createLSystemUI", "deleteLSystemUI");

	// Register Command
	status = plugin.registerCommand( "LSystemCmd", LSystemCmd::creator, LSystemCmd::newSyntax );
    if (!status) {
        status.perror("registerCommand");
        return status;
    }

	// Register Node
	status = plugin.registerNode("LSystemNode", LSystemNode::id,
								  LSystemNode::creator, LSystemNode::initialize);
	if (!status) {
		status.perror("registerNode");
		return status;
	}

    return status;
    return status;
}

MStatus uninitializePlugin( MObject obj)
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj );

    status = plugin.deregisterCommand( "LSystemCmd" );
    if (!status) {
	    status.perror("deregisterCommand");
	    return status;
    }

	status = plugin.deregisterNode(LSystemNode::id);
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

    return status;
}


