#include "LiuSalonCmd.h"
#include "cyHairFile.h"

#include <maya/MGlobal.h>
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

	int numStrands = 0;
	int numPoints = 0;
	double hairLength = 0;
	MArgDatabase argData(syntax(), args);
	h->Initialize();

    if(argData.isFlagSet(strandsFlag))
		argData.getFlagArgument(strandsFlag, 0, numStrands);
    if(argData.isFlagSet(pointsFlag))
		argData.getFlagArgument(pointsFlag, 0, numPoints);
    if(argData.isFlagSet(lengthFlag))
		argData.getFlagArgument(lengthFlag, 0, hairLength);

	h->CreateHair(numStrands);
	h->CreatePoints(numStrands*numPoints);
	unsigned short* segments = h->GetSegmentsArray();
	float* points = h->GetPointsArray();

	std::stringstream deg;
	deg << " -d " << numPoints-1;
	MString degreeString(deg.str().c_str());
	int index = 0;
	int x = 0;
	float ystart = -hairLength/2; // where hair strand starts drawing
	float segLength = hairLength/(numPoints-1); // hair segment length
	for (int i = 0; i < numStrands; i++)
	{
		std::stringstream pts;
		segments[i] = numPoints-1; // TODO: probably unnecessary
		for (int j = 0; j < numPoints; j++)
		{
			points[index] = x;
			points[index+1] = ystart + segLength*j;
			points[index+2] = 0;
			pts << " -p " << points[index] << " " << points[index+1] << " " << points[index+2];
			index += 3;
		}
		x++;
		MString pointString(pts.str().c_str());
		MGlobal::executeCommand("curve" + degreeString + pointString);
	}
	return MStatus::kSuccess;
}
