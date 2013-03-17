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
	// now create the hair
	cyHairFile* h = new cyHairFile();
	int hairCount = h->LoadFromFile("C:/Users/adair/Documents/GitHub/liusalon/LiuSalon/LiuSalon/hairFiles/straight.hair");
	MGlobal::executeCommand(MString("print ")+ hairCount);


	int currPtNum = 0;
	//unsigned short* temp = h->GetSegmentsArray();

	//unsigned short t2 = temp[0];
	//MGlobal::executeCommand(MString("print ")+ h->GetSegmentsArray()[0]);
	//-fpt true - the resulting surface at the path curve
	int hasSegments = h->GetHeader().arrays & CY_HAIR_FILE_SEGMENTS_BIT;
	int hasThickness = h->GetHeader().arrays & CY_HAIR_FILE_THICKNESS_BIT;

	int num_of_segments = h->GetHeader().d_segments;

	hairCount = 50;

	for(int i=0; i< hairCount; ++i){
		

		if(hasSegments)
			num_of_segments = h->GetSegmentsArray()[i];

		for(unsigned int j=0; j<num_of_segments;++j)
		{	
		
			MGlobal::executeCommand(MString("circle -radius ") + (hasThickness? h->GetThicknessArray()[currPtNum] : h->GetHeader().d_thickness) );
			
			int currPtIndex = currPtNum * 3;

			float scale = hasThickness? h->GetThicknessArray()[currPtNum+1] / h->GetThicknessArray()[currPtNum] : 1.0;

			MGlobal::executeCommand(MString("curve -d 1 -p ") 
				+ h->GetPointsArray()[currPtIndex] + " " + h->GetPointsArray()[currPtIndex+1] + " " + h->GetPointsArray()[currPtIndex+2]
				+ " -p " 
				+ h->GetPointsArray()[currPtIndex+3]  + " " + h->GetPointsArray()[currPtIndex+4]  + " " + h->GetPointsArray()[currPtIndex+5] 
				+" -k 0 -k 1 -name curve1;");
			MGlobal::executeCommand(MString("select -r nurbsCircle1 curve1;"));
			MGlobal::executeCommand(MString("extrude -po 1 -et 2 -ucp 1 -fpt true -upn true -sc ") + scale
										+ " -rsp 1 \"nurbsCircle1\" \"curve1\";");
			MGlobal::executeCommand(MString("select -r nurbsCircle1; doDelete;"));
			MGlobal::executeCommand(MString("select -r curve1; doDelete;"));
			
			currPtNum ++;
		}
		currPtNum ++;
	}

    return MStatus::kSuccess;
	return stat;
}
