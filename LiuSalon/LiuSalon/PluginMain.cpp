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
#include "HairModelNode.h"

MString filePath;
MStatus initializePlugin( MObject obj )
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj, "LiuSalon", "1.0", "Any");

	// load external ui
    MGlobal::executeCommand("source \"" + plugin.loadPath() + "/ui.mel\"");
	status = plugin.registerUI("createLiuSalonUI", "deleteLiuSalonUI");

	filePath = plugin.loadPath();

	// Register Command
	status = plugin.registerCommand( "LiuSalonCmd", LiuSalonCmd::creator, LiuSalonCmd::newSyntax );
    if (!status) {
        status.perror("registerCommand");
        return status;
    }

	// Register Node
	status = plugin.registerNode("HairModelNode", HairModelNode::id,
								  HairModelNode::creator, HairModelNode::initialize);
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

    status = plugin.deregisterCommand( "LiuSalonCmd" );
    if (!status) {
	    status.perror("deregisterCommand");
	    return status;
    }

	status = plugin.deregisterNode(HairModelNode::id);
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

    return status;
}


