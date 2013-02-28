#include "LiuSalonCmd.h"
#include "cyHairFile.h"

#include <maya/MArgDatabase.h>
#include <list>
#include <sstream>

// command line args
const char *strandsFlag = "-s", *strandsLongFlag = "-strands";
const char *pointsFlag = "-p", *pointsLongFlag = "-points";
const char *lengthFlag = "-l", *lengthLongFlag = "-length";

LiuSalonCmd::LiuSalonCmd() : MPxCommand()
{
}

LiuSalonCmd::~LiuSalonCmd() 
{
}

MSyntax LiuSalonCmd::newSyntax()
{
    MSyntax syntax;

    syntax.addFlag(strandsFlag, strandsLongFlag, MSyntax::kLong);
	syntax.addFlag(pointsFlag, pointsLongFlag, MSyntax::kLong);
	syntax.addFlag(lengthFlag, lengthLongFlag, MSyntax::kDouble);

	return syntax;
}

MStatus LiuSalonCmd::doIt( const MArgList& args )
{
	MStatus stat = MS::kSuccess;
	setResult( "LiuSalon command executed!\n" );

	// create cyHairFile object for storing hair data
	cyHairFile* h = new cyHairFile();

	int numStrands;
	int numPoints;
	double hairLength;
	MArgDatabase argData(syntax(), args);

    if(argData.isFlagSet(strandsFlag)) {
		argData.getFlagArgument(strandsFlag, 0, numStrands);
		h->SetHairCount(numStrands);
	}
    if(argData.isFlagSet(pointsFlag)) {
		argData.getFlagArgument(pointsFlag, 0, numPoints);
		h->SetPointCount(numPoints);
	}
    if(argData.isFlagSet(lengthFlag))
		argData.getFlagArgument(lengthFlag, 0, hairLength);

	// now create the hair

    return MStatus::kSuccess;
	return stat;
}
