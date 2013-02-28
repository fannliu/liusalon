#include "LiuSalonCmd.h"

#include <maya/MSimple.h>

//DeclareSimpleCommand( LiuSalon, "", "2012");

// command line args
// TODO: fill in

LiuSalonCmd::LiuSalonCmd() : MPxCommand()
{
}

LiuSalonCmd::~LiuSalonCmd() 
{
}

MSyntax LiuSalonCmd::newSyntax()
{
    MSyntax syntax;

    // TODO: fill in

	return syntax;
}

MStatus LiuSalonCmd::doIt( const MArgList& args )
{
	MStatus stat = MS::kSuccess;
	//setResult( "LiuSalon command executed!\n" );

	return stat;
}
