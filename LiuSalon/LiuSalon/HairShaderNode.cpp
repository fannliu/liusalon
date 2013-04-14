#include "HairShaderNode.h"
#include <fstream>

#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg;				\
		return MS::kFailure;		\
	}

#define MNoVersionString
#define MNoPluginEntry
#include <maya/MFnPlugin.h>
#include <maya/MFnNurbsCurve.h>

#define MAKE_INPUT(attr) \
CHECK_MSTATUS(attr.setKeyable(true)); \
CHECK_MSTATUS(attr.setStorable(true)); \
CHECK_MSTATUS(attr.setReadable(true)); \
CHECK_MSTATUS(attr.setWritable(true));
#define MAKE_INPUT_COLOR(attr) \
CHECK_MSTATUS (attr.setDefault(1.0f, 1.0f, 1.0f)); \
CHECK_MSTATUS (attr.setUsedAsColor(true)); \
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

// primary highlight attributes (R)
MObject HairShaderNode::inColorR;
MObject HairShaderNode::outColorR;
MObject HairShaderNode::inIntensityR;
MObject HairShaderNode::outIntensityR;
MObject HairShaderNode::inLongShiftR;
MObject HairShaderNode::outLongShiftR;
MObject HairShaderNode::inLongWidthR;
MObject HairShaderNode::outLongWidthR;

// backlit light attributes (TT)
MObject HairShaderNode::inColorTT;
MObject HairShaderNode::outColorTT;
MObject HairShaderNode::inIntensityTT;
MObject HairShaderNode::outIntensityTT;
MObject HairShaderNode::inLongShiftTT;
MObject HairShaderNode::outLongShiftTT;
MObject HairShaderNode::inLongWidthTT;
MObject HairShaderNode::outLongWidthTT;
MObject HairShaderNode::inAzimuthWidthTT;
MObject HairShaderNode::outAzimuthWidthTT;

// secondary highlight attributes (TRT-G)
MObject HairShaderNode::inColorTRT; // this will affect glints as well
MObject HairShaderNode::outColorTRT;
MObject HairShaderNode::inIntensityTRT;
MObject HairShaderNode::outIntensityTRT;
MObject HairShaderNode::inLongShiftTRT;
MObject HairShaderNode::outLongShiftTRT;
MObject HairShaderNode::inLongWidthTRT;
MObject HairShaderNode::outLongWidthTRT;

// glint attributes (G)
MObject HairShaderNode::inIntensityG;
MObject HairShaderNode::outIntensityG;
MObject HairShaderNode::inAzimuthWidthG; // frequency of glints
MObject HairShaderNode::outAzimuthWidthG;

MTypeId HairShaderNode::id( 0x80001 );

void* HairShaderNode::creator()
{
	return new HairShaderNode;
}

MStatus HairShaderNode::initialize()
{
	MFnNumericAttribute nAttr;

	MStatus returnStatus;

	// user defined attributes for R
	HairShaderNode::inColorR = nAttr.createColor("delightcolorR", "cr");
	MAKE_INPUT_COLOR(nAttr);

	HairShaderNode::inIntensityR = nAttr.create("intensityR", "ir",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating intensity R attribute\n");
	MAKE_INPUT(nAttr);

	HairShaderNode::inLongShiftR = nAttr.create("longitudinalShiftR", "lsr", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating longitudinal shift R attribute\n");
	MAKE_INPUT(nAttr);

	HairShaderNode::inLongWidthR = nAttr.create("longitudinalWidthR", "lwr", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating longitudinal width R attribute\n");
	MAKE_INPUT(nAttr);

	// user defined attributes for TT
	HairShaderNode::inColorTT = nAttr.createColor("colorTT", "ctt");
	MAKE_INPUT_COLOR(nAttr);

	HairShaderNode::inIntensityTT = nAttr.create("intensityTT", "itt",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating intensity TT attribute\n");
	MAKE_INPUT(nAttr);

	HairShaderNode::inLongShiftTT = nAttr.create("longitudinalShiftTT", "lstt", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating longitudinal shift TT attribute\n");
	MAKE_INPUT(nAttr);

	HairShaderNode::inLongWidthTT = nAttr.create("longitudinalWidthTT", "lwtt", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating longitudinal width TT attribute\n");
	MAKE_INPUT(nAttr);	

	HairShaderNode::inAzimuthWidthTT = nAttr.create("azimuthalWidthTT", "awtt", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating azimuthal width TT attribute\n");
	MAKE_INPUT(nAttr);	

	// user defined attributes for TRT-G
	HairShaderNode::inColorTRT = nAttr.createColor("colorTRT", "ctrt");
	MAKE_INPUT_COLOR(nAttr);

	HairShaderNode::inIntensityTRT = nAttr.create("intensityTRT", "itrt",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating intensity TRT attribute\n");
	MAKE_INPUT(nAttr);

	HairShaderNode::inLongShiftTRT = nAttr.create("longitudinalShiftTRT", "lstrt",
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating longitudinal shift TRT attribute\n");
	MAKE_INPUT(nAttr);

	HairShaderNode::inLongWidthTRT = nAttr.create("longitudinalWidthTRT", "lwtrt", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating longitudinal width TRT attribute\n");
	MAKE_INPUT(nAttr);	

	// user defined attributes for G
	HairShaderNode::inIntensityG = nAttr.create("intensityG", "ig",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating intensity G attribute\n");
	MAKE_INPUT(nAttr);

	HairShaderNode::inAzimuthWidthG = nAttr.create("azimuthalWidthG", "lwg", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating azimuthal width G attribute\n");
	MAKE_INPUT(nAttr);

	// create the output attributes for R
	HairShaderNode::outColorR = nAttr.createColor("outColorR", "ocr");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outIntensityR = nAttr.create("outIntensityR", "oir",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating output intensity R attribute\n");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outLongShiftR = nAttr.create("outLongshiftR", "olsr", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output longitudinal shift R attribute\n");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outLongWidthR = nAttr.create("outLongwidthR", "olwr", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output longitudinal width R attribute\n");
	MAKE_OUTPUT(nAttr);
	
	// create the output attributes for T
	HairShaderNode::outColorTT = nAttr.createColor("outColorTT", "octt");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outIntensityTT = nAttr.create("outIntensityTT", "oitt",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating output intensity TT attribute\n");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outLongShiftTT = nAttr.create("outLongshiftTT", "olstt", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output longitudinal shift TT attribute\n");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outLongWidthTT = nAttr.create("outLongwidthTT", "olwtt", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output longitudinal width TT attribute\n");
	MAKE_OUTPUT(nAttr);	

	HairShaderNode::outAzimuthWidthTT = nAttr.create("outAzimuthwidthTT", "oawtt", 
									MFnNumericData::kDouble,
									0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output azimuthal width TT attribute\n");
	MAKE_OUTPUT(nAttr);	

	// create the output attributes for TRT-G
	HairShaderNode::outColorTRT = nAttr.createColor("outColorTRT", "octrt");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outIntensityTRT = nAttr.create("outIntensityTRT", "oitrt",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating output intensity TRT attribute\n");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outLongShiftTRT = nAttr.create("outLongshiftTRT", "olstrt",
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output longitudinal shift TRT attribute\n");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outLongWidthTRT = nAttr.create("outLongwidthTRT", "olwtrt", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output longitudinal width TRT attribute\n");
	MAKE_OUTPUT(nAttr);

	// create the output attributes for G
	HairShaderNode::outIntensityG = nAttr.create("outIntensityG", "oig",
								MFnNumericData::kDouble,
								0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating output intensity G attribute\n");
	MAKE_OUTPUT(nAttr);

	HairShaderNode::outAzimuthWidthG = nAttr.create("outAzimuthwidthG", "olwg", 
								MFnNumericData::kDouble,
								0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output azimuthal width G attribute\n");
	MAKE_OUTPUT(nAttr);

	// add R attributes
	returnStatus = addAttribute(HairShaderNode::inColorR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inColorR attribute\n");
	returnStatus = addAttribute(HairShaderNode::outColorR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outColorR attribute\n");
	returnStatus = addAttribute(HairShaderNode::inIntensityR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inIntensityR attribute\n");
	returnStatus = addAttribute(HairShaderNode::outIntensityR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outIntensityR attribute\n");
	returnStatus = addAttribute(HairShaderNode::inLongShiftR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inLongShiftR attribute\n");
	returnStatus = addAttribute(HairShaderNode::outLongShiftR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outLongShiftR attribute\n");
	returnStatus = addAttribute(HairShaderNode::inLongWidthR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inLongWidthR attribute\n");
	returnStatus = addAttribute(HairShaderNode::outLongWidthR);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outLongWidthR attribute\n");

	// add TT attributes
	returnStatus = addAttribute(HairShaderNode::inColorTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inColorTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outColorTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outColorTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::inIntensityTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inIntensityTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outIntensityTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outIntensityTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::inLongShiftTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inLongShiftTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outLongShiftTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outLongShiftTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::inLongWidthTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inLongWidthTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outLongWidthTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outLongWidthTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::inAzimuthWidthTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inAzimuthWidthTT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outAzimuthWidthTT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outAzimuthWidthTT attribute\n");

	// add TRT-G attributes
	returnStatus = addAttribute(HairShaderNode::inColorTRT); // this will affect glints as well
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inColorTRT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outColorTRT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outColorTRT attribute\n");
	returnStatus = addAttribute(HairShaderNode::inIntensityTRT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inIntensityTRT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outIntensityTRT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outIntensityTRT attribute\n");
	returnStatus = addAttribute(HairShaderNode::inLongShiftTRT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inLongShiftTRT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outLongShiftTRT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outLongShiftTRT attribute\n");
	returnStatus = addAttribute(HairShaderNode::inLongWidthTRT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inLongWidthTRT attribute\n");
	returnStatus = addAttribute(HairShaderNode::outLongWidthTRT);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outLongWidthTRT attribute\n");
	
	// add G attributes
	returnStatus = addAttribute(HairShaderNode::inIntensityG);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inIntensityG attribute\n");
	returnStatus = addAttribute(HairShaderNode::outIntensityG);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outIntensityG attribute\n");
	returnStatus = addAttribute(HairShaderNode::inAzimuthWidthG); // frequency of glints
	McheckErr(returnStatus, "ERROR adding HairShaderNode::inAzimuthWidthG attribute\n");
	returnStatus = addAttribute(HairShaderNode::outAzimuthWidthG);
	McheckErr(returnStatus, "ERROR adding HairShaderNode::outAzimuthWidthG attribute\n");

	// attribute affects R
	returnStatus = attributeAffects(HairShaderNode::inColorR,
								    HairShaderNode::outColorR);
	returnStatus = attributeAffects(HairShaderNode::inIntensityR,
								    HairShaderNode::outIntensityR);
	returnStatus = attributeAffects(HairShaderNode::inLongShiftR,
								    HairShaderNode::outLongShiftR);
	returnStatus = attributeAffects(HairShaderNode::inLongWidthR,
								    HairShaderNode::outLongWidthR);

	// attribute affects TT
	returnStatus = attributeAffects(HairShaderNode::inColorTT,
								    HairShaderNode::outColorTT);
	returnStatus = attributeAffects(HairShaderNode::inIntensityTT,
								    HairShaderNode::outIntensityTT);
	returnStatus = attributeAffects(HairShaderNode::inLongShiftTT,
								    HairShaderNode::outLongShiftTT);
	returnStatus = attributeAffects(HairShaderNode::inLongWidthTT,
								    HairShaderNode::outLongWidthTT);
	returnStatus = attributeAffects(HairShaderNode::inAzimuthWidthTT,
								    HairShaderNode::outAzimuthWidthTT);

	// attribute affects TRT
	returnStatus = attributeAffects(HairShaderNode::inColorTRT,
								    HairShaderNode::outColorTRT);
	returnStatus = attributeAffects(HairShaderNode::inIntensityTRT,
								    HairShaderNode::outIntensityTRT);
	returnStatus = attributeAffects(HairShaderNode::inLongShiftTRT,
								    HairShaderNode::outLongShiftTRT);
	returnStatus = attributeAffects(HairShaderNode::inLongWidthTRT,
								    HairShaderNode::outLongWidthTRT);

	// attribute affects glints
	returnStatus = attributeAffects(HairShaderNode::inIntensityG,
								    HairShaderNode::outIntensityG);
	returnStatus = attributeAffects(HairShaderNode::inAzimuthWidthG,
								    HairShaderNode::outAzimuthWidthG);

	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	return MS::kSuccess;
}

MStatus HairShaderNode::compute(const MPlug& plug, MDataBlock& data)
{
	// call this function whenever attribute changes
	// recompute the output mesh

	MStatus returnStatus;

	// giant if statement for all the input/outputs
	if (plug == outIntensityG) {

	}
	else
		return MS::kUnknownParameter;
	/*if (plug == outputMesh) {
		cyHairFile* h = new cyHairFile();

		// get file
		MDataHandle fileData = data.inputValue( file, &returnStatus ); 
		McheckErr(returnStatus, "Error getting file data handle\n");
		std::string fileName = fileData.asString().asChar();
		
		if (!fileName.empty())
		{
			MObject pluginObj = MFnPlugin::findPlugin("LiuSalon");
			MFnPlugin plugin(pluginObj);
			
			// find the file
			cout << fileName.c_str() << endl;
			int hairCount = h->LoadFromFile(fileName.c_str());

			// get hair strands
			MDataHandle strandsData = data.inputValue( numStrands, &returnStatus ); 
			McheckErr(returnStatus, "Error getting strands data handle\n");
			int st = strandsData.asInt();

			if ( st > 0 && st < hairCount ) {
				h->SetHairCount(st); //set number of hair strand
			}
		}
		else
		{
			// if no file then get other data for user created hair
			// get num points per strand
			MDataHandle pointsData = data.inputValue( numPoints, &returnStatus ); 
			McheckErr(returnStatus, "Error getting points data handle\n");
			int pts = pointsData.asInt();

			// get hair strands
			MDataHandle strandsData = data.inputValue( numStrands, &returnStatus ); 
			McheckErr(returnStatus, "Error getting strands data handle\n");
			int st = strandsData.asInt();

			// get hair length
			MDataHandle lengthData = data.inputValue( hairLength, &returnStatus ); 
			McheckErr(returnStatus, "Error getting length data handle\n");
			double len = lengthData.asDouble();

			// get curves array
			MArrayDataHandle inputArray = data.inputArrayValue( inputCurve,
																&returnStatus );
			McheckErr(returnStatus, "Error getting curve data handle\n");
			int numCurves = inputArray.elementCount();
			// now use the curves values
			if (numCurves >= 1)
			{
				createHairCurve(inputArray, *h);
			}
			else
			{
				h->SetDefaultSegmentCount(pts - 1);
				h->SetHairCount(st);
				h->CreatePoints(st*pts);
				// fill in random points
				int index = 0;
				int x = 0;
				float ystart = -len/2; // where hair strand starts drawing
				float segLength = len/(pts-1); // hair segment length

				float* points = h->GetPointsArray();

				for (int i = 0; i < st; i++)
				{
					for (int j = 0; j < pts; j++)
					{
						points[index] = x;
						points[index+1] = ystart + segLength*j;
						points[index+2] = 0;
						index += 3;
					}
					x++;
				}
			}
		}
		MDataHandle outputHandle = data.outputValue(outputMesh, &returnStatus);
		McheckErr(returnStatus, "ERROR getting polygon data handle\n");

		MFnMeshData dataCreator;
		MObject newOutputData = dataCreator.create(&returnStatus);
		McheckErr(returnStatus, "ERROR creating outputData");
			
		createMesh(newOutputData, returnStatus, *h);
		McheckErr(returnStatus, "ERROR creating new strand");
		
		outputHandle.set(newOutputData);
		data.setClean( plug );
	}
	else
		return MS::kUnknownParameter;
		*/
	return MS::kSuccess;

}
