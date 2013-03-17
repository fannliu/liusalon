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

MObject HairModelNode::outputMesh;
MObject HairModelNode::numStrands;
MObject HairModelNode::numPoints;
MObject HairModelNode::hairLength;
MObject HairModelNode::file;

MTypeId HairModelNode::id( 0x80000 );

void* HairModelNode::creator()
{
	return new HairModelNode;
}

MStatus HairModelNode::initialize()
{
	MFnNumericAttribute strandsAttr;
	MFnNumericAttribute pointsAttr;
	MFnNumericAttribute lengthAttr;
	MFnTypedAttribute fileAttr;
	MFnTypedAttribute outAttr;

	MStatus returnStatus;

	HairModelNode::numStrands = strandsAttr.create( "strands", "s",
											MFnNumericData::kInt,
											0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating HairModelNode strands attribute\n");
	MAKE_INPUT(strandsAttr);

	HairModelNode::numPoints = pointsAttr.create( "pointsPerStrand", "p",
										  MFnNumericData::kInt,
										  0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating HairModelNode points attribute\n");
	MAKE_INPUT(pointsAttr);

	HairModelNode::hairLength = lengthAttr.create( "lengthOfStrand", "l",
										  MFnNumericData::kDouble,
										  0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating HairModelNode length attribute\n");
	MAKE_INPUT(lengthAttr);

	HairModelNode::file = fileAttr.create( "HAIRfile", "f",
											MFnData::kString,
											&returnStatus ); 
	McheckErr(returnStatus, "ERROR creating HairModelNode output attribute\n");
	MAKE_INPUT(fileAttr);

	HairModelNode::outputMesh = outAttr.create( "outputMesh", "out",
												 MFnData::kMesh,
												 &returnStatus ); 
	McheckErr(returnStatus, "ERROR creating HairModelNode output attribute\n");
	outAttr.setStorable(false);

	returnStatus = addAttribute(HairModelNode::outputMesh);
	McheckErr(returnStatus, "ERROR adding outputMesh attribute\n");

	returnStatus = addAttribute(HairModelNode::numPoints);
	McheckErr(returnStatus, "ERROR adding points attribute\n");

	returnStatus = addAttribute(HairModelNode::numStrands);
	McheckErr(returnStatus, "ERROR adding strands attribute\n");

	returnStatus = addAttribute(HairModelNode::hairLength);
	McheckErr(returnStatus, "ERROR adding length attribute\n");

	returnStatus = addAttribute(HairModelNode::file);
	McheckErr(returnStatus, "ERROR adding file attribute\n");

	returnStatus = attributeAffects(HairModelNode::numPoints,
								    HairModelNode::outputMesh);

	returnStatus = attributeAffects(HairModelNode::numStrands,
									HairModelNode::outputMesh);

	returnStatus = attributeAffects(HairModelNode::hairLength,
								    HairModelNode::outputMesh);

	returnStatus = attributeAffects(HairModelNode::file,
								    HairModelNode::outputMesh);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	return MS::kSuccess;
}

MObject HairModelNode::createMesh(MObject& outData, MStatus& stat, cyHairFile& hair)
{
	MFnMesh mesh;
	MFnMeshData meshData;
	MObject meshObject = meshData.create();
	int numVertices;
	int numPolygons;

	MPointArray vertexArray;  // The vertices of the geometry.
	MIntArray polygonCounts;  // Array of # of vertices per face.
	MIntArray polygonConnects;  // Basically an Index Buffer Object.

	float* points = hair.GetPointsArray();
	unsigned short* segments = hair.GetSegmentsArray();
	cyHairFileHeader header = hair.GetHeader();
	int numStrands = header.hair_count;
	int index = 0;
	int numSegments;

	// The magic to fill vertexArray, polygonCounts, and polygonConnects goes here.
	for (int i = 0; i < numStrands; i++) 
	{
		// get number of segments on this strand
		numSegments = segments[i];
		for (int j = 0; j < numSegments; j++)
		{
			MPointArray cvertexArray;  // The vertices of the geometry.
			MIntArray cpolygonCounts;  // Array of # of vertices per face.
			MIntArray cpolygonConnects;  // Basically an Index Buffer Object.
			CylinderMesh newCylinder(MPoint(points[index], points[index+1], points[index+2]), MPoint(points[index+3], points[index+4], points[index+5]), 0.5, 0.5);
			newCylinder.appendToMesh(vertexArray, polygonCounts, polygonConnects);
		
			index += 3;
		}

		index += 3; // skip the last point
	}
	numVertices = vertexArray.length();
	numPolygons = polygonCounts.length();

	MObject newMesh = mesh.create(numVertices, numPolygons, vertexArray, polygonCounts,
								  polygonConnects, outData, &stat);
	return newMesh;
}

MStatus HairModelNode::compute(const MPlug& plug, MDataBlock& data)
{
	// call this function whenever attribute changes
	// recompute the output mesh
	MStatus returnStatus;
	if (plug == outputMesh) {
		cyHairFile* h = new cyHairFile();

		// TODO: handle HAIR file
		// get file
		MDataHandle fileData = data.inputValue( file, &returnStatus ); 
		McheckErr(returnStatus, "Error getting file data handle\n");
		std::string fileName = fileData.asString().asChar();

		if (!fileName.empty())
		{
			// find the actual file
			MObject pluginObj = MFnPlugin::findPlugin("LiuSalon");
			MFnPlugin plugin(pluginObj);
		}
		else
		{
			// if no file then get other data for user created hair
			// get num points
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

			// TODO: add a cube as an input attribute

			h->CreateHair(st);
			h->CreatePoints(st*pts);
			unsigned short* segments = h->GetSegmentsArray();
			float* points = h->GetPointsArray();

			// fill in random points
			int index = 0;
			int x = 0;
			float ystart = -len/2; // where hair strand starts drawing
			float segLength = len/(pts-1); // hair segment length

			for (int i = 0; i < st; i++)
			{
				segments[i] = pts-1; // TODO: probably unnecessary
				for (int j = 0; j < pts; j++)
				{
					points[index] = x;
					points[index+1] = ystart + segLength*j;
					points[index+2] = 0;
					index += 3;
				}
				x++;
			}

			/* Get output object */
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
	}
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}
