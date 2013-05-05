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
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>

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
#define MAKE_CURVE_OUTPUT(attr) \
CHECK_MSTATUS(attr.setArray(true)); \
CHECK_MSTATUS(attr.setReadable(true)); \
CHECK_MSTATUS(attr.setWritable(false)); \
CHECK_MSTATUS (attr.setUsesArrayDataBuilder(true));
#define MAKE_ADDR(attr) \
CHECK_MSTATUS(attr.setKeyable(false)); \
CHECK_MSTATUS(attr.setStorable(false)); \
CHECK_MSTATUS(attr.setReadable(true)); \
CHECK_MSTATUS(attr.setWritable(false)); \
CHECK_MSTATUS(attr.setHidden(true));

MObject HairModelNode::inputCurve;
MObject HairModelNode::outputCurves;
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
	MFnTypedAttribute curveAttr;

	MStatus returnStatus;

	// for selecting an input curve and translating to framework
	HairModelNode::inputCurve = curveAttr.create("inputCurve", "in", 
										MFnData::kNurbsCurve, 
										&returnStatus);
	McheckErr(returnStatus, "ERROR creating curve attribute\n");
	MAKE_INPUT(curveAttr);

	// set it as an array to get an array of curves
	returnStatus = curveAttr.setArray(true);
	McheckErr(returnStatus, "ERROR setting curve attribute array\n");

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
    
	HairModelNode::outputCurves = outAttr.create( "outputCurves", "oc",
													MFnNurbsCurveData::kNurbsCurve,
													&returnStatus );
	
	//HairModelNode::outputMesh = outAttr.create( "outputMesh", "out",
	//											 MFnData::kMesh,
	//											 &returnStatus ); 
												 
	McheckErr(returnStatus, "ERROR creating HairModelNode output attribute\n");
	//MAKE_OUTPUT(outAttr);
	MAKE_CURVE_OUTPUT(outAttr);

	returnStatus = addAttribute(HairModelNode::inputCurve);
	McheckErr(returnStatus, "ERROR adding inputCurve attribute\n");

	returnStatus = addAttribute(HairModelNode::outputCurves);
	McheckErr(returnStatus, "ERROR adding outputCurves attribute\n");

//	returnStatus = addAttribute(HairModelNode::outputMesh);
//	McheckErr(returnStatus, "ERROR adding outputMesh attribute\n");
	

//	returnStatus = addAttribute(HairModelNode::numPoints);
//	McheckErr(returnStatus, "ERROR adding points attribute\n");

	returnStatus = addAttribute(HairModelNode::numStrands);
	McheckErr(returnStatus, "ERROR adding strands attribute\n");

//	returnStatus = addAttribute(HairModelNode::hairLength);
//	McheckErr(returnStatus, "ERROR adding length attribute\n");

	returnStatus = addAttribute(HairModelNode::file);
	McheckErr(returnStatus, "ERROR adding file attribute\n");

	returnStatus = attributeAffects(HairModelNode::inputCurve,
								    HairModelNode::outputCurves);

//	returnStatus = attributeAffects(HairModelNode::numPoints,
//								    HairModelNode::outputCurves);

	returnStatus = attributeAffects(HairModelNode::numStrands,
									HairModelNode::outputCurves);

//	returnStatus = attributeAffects(HairModelNode::hairLength,
//								    HairModelNode::outputCurves);

	returnStatus = attributeAffects(HairModelNode::file,
								    HairModelNode::outputCurves);

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
	float* thickness = hair.GetThicknessArray();
	unsigned short* segments = hair.GetSegmentsArray();
	cyHairFileHeader header = hair.GetHeader();

	int numStrands = header.hair_count; 
	int index = 0;
	int numSegments = header.d_segments;

	int hasSegments = ((header.arrays & CY_HAIR_FILE_SEGMENTS_BIT) || !hair.use_default_seg);
	int hasThickness = header.arrays & CY_HAIR_FILE_THICKNESS_BIT;
	
	if (hasSegments)
		cout << "has segments!" << endl;
	cout << numStrands << endl;
	// The magic to fill vertexArray, polygonCounts, and polygonConnects goes here.
	for (int i = 0; i < numStrands; i++) 
	{
		// get number of segments on this strand
		if(hasSegments)
			numSegments = segments[i];
		for (int j = 0; j <= numSegments; j++)
		{
			MPointArray cvertexArray;  // The vertices of the geometry.
			MIntArray cpolygonCounts;  // Array of # of vertices per face.
			MIntArray cpolygonConnects;  // Basically an Index Buffer Object.
			float r1 = header.d_thickness;
			float r2 = r1;

			if(hasThickness){
				r1 = thickness[ index/3     ];
				r2 = thickness[ index/3 + 1 ];
			}
			CylinderMesh newCylinder(MPoint(points[index], points[index+1], points[index+2]), MPoint(points[index+3], points[index+4], points[index+5]), r1, r2);
			newCylinder.appendToMesh(vertexArray, polygonCounts, polygonConnects);

			index += 3;
		}
	}

	numVertices = vertexArray.length();
	numPolygons = polygonCounts.length();

	MObject newMesh = mesh.create(numVertices, numPolygons, vertexArray, polygonCounts,
								  polygonConnects, outData, &stat);
	return newMesh;
}

// this is translating a user defined curve to a HAIR curve
MStatus HairModelNode::createHairCurve(MArrayDataHandle &inputArray, cyHairFile& hair)
{
	MStatus returnStatus;
	int numCurves = inputArray.elementCount();
	// LOOK: Count the number of CVs, uncomment if necessary 
	int numCVs;
	int numKnots;

	inputArray.jumpToElement(0);
	MDataHandle elementHandle = inputArray.inputValue(&returnStatus);
	McheckErr(returnStatus, "createHairCurve:: input value\n");
	MObject countCurve (elementHandle.asNurbsCurve());
	MFnNurbsCurve countCurveFn (countCurve);
	numKnots = countCurveFn.numKnots (&returnStatus);
	McheckErr(returnStatus, "createHairCurve:: counting knots\n");

	if (numCurves >= 1)
	{
		// number of strands = number of curves
		hair.CreateHair(numCurves);

		// need to count total number of knots for allocating array
		int totalKnots = 0;
		int totalCVs = 0;
		int currKnots;
		int currCVs;
		
		for (int curveNum = 0; curveNum < numCurves; curveNum++)
		{
			MObject curve (inputArray.inputValue().asNurbsCurve());
			MFnNurbsCurve curveFn( curve, &returnStatus );
			McheckErr(returnStatus, "createHairCurve:: error creating curve function set\n");
			
			MObject countCurve (elementHandle.asNurbsCurve());
			MFnNurbsCurve countCurveFn (countCurve);
			currKnots = countCurveFn.numKnots (&returnStatus);
			McheckErr(returnStatus, "createHairCurve:: counting knots at curve loop\n");
			currCVs = countCurveFn.numCVs(&returnStatus);
			McheckErr(returnStatus, "createHairCurve:: counting CVs at curve loop\n");

			totalCVs += currCVs;
			totalKnots += currKnots;
			returnStatus = inputArray.next();
		}

		// allocate array
		hair.CreatePoints(totalKnots);

		// go back to the first element in array
		inputArray.jumpToElement(0);
		
		unsigned short* segments = hair.GetSegmentsArray();
		float* points = hair.GetPointsArray();
		int index = 0;

		// TODO: add CVs if necessary, use CVs and numSegments for smoothing?
		for (int curveNum = 0; curveNum < numCurves; curveNum++)
		{
			MObject curve (inputArray.inputValue().asNurbsCurve());
			MFnNurbsCurve curveFn( curve, &returnStatus );
			McheckErr(returnStatus, "createHairCurve:: error creating curve function set\n");
			MDoubleArray knots;
			
			returnStatus = curveFn.getKnots(knots);
			McheckErr(returnStatus, "createHairCurve:: error getting curve knots\n");

			// knots are parameters, use to get points
			// build the hair array
			// get the number of knots on this particular curve
			numKnots = curveFn.numKnots (&returnStatus);
			numCVs = curveFn.numCVs(&returnStatus);
			McheckErr(returnStatus, "createHairCurve:: counting knots in loop\n");

			// fill in number of segments using number of knots
			segments[curveNum] = numKnots - 1;

			for (int i = 0; i < numKnots; i++)
			{
				MPoint p;
				//curveFn.getCV(i, p, MSpace::kWorld);
				curveFn.getPointAtParam(knots[i], p, MSpace::kWorld);
				points[index] = p[0];
				points[index+1] = p[1];
				points[index+2] = p[2];
				index += 3;
			}

			if (knots.length() != (unsigned)numKnots)
				returnStatus = MS::kFailure;
			McheckErr(returnStatus, "createHairCurve:: inconsistent number of knots, rebuild curves\n");

			returnStatus = inputArray.next();
		}
	}
	return MS::kSuccess;
}

MStatus HairModelNode::compute(const MPlug& plug, MDataBlock& data)
{
	// call this function whenever attribute changes
	// recompute the output curves

	MStatus returnStatus;
	cout << "testtesttest" << endl;
	if (plug == outputCurves)
	{
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
			cerr << hairCount;

			// get hair strands
			MDataHandle strandsData = data.inputValue( numStrands, &returnStatus ); 
			McheckErr(returnStatus, "Error getting strands data handle\n");
			int st = strandsData.asInt();

			if ( st > 0 && st < hairCount )
			{
				h->SetHairCount(st); //set number of hair strand
			}
			MArrayDataHandle outputArray = data.outputArrayValue( outputCurves, &returnStatus );
			McheckErr(returnStatus, "Error getting output data handle\n");

			MArrayDataBuilder builder(outputCurves, st, &returnStatus);
			McheckErr(returnStatus, "Error creating builder\n");
 
			int hasSegments = h->GetHeader().arrays & CY_HAIR_FILE_SEGMENTS_BIT;
			int num_of_segments = h->GetHeader().d_segments;	
			
			int currPtIndex= 0;
			for (int curveNum = 0; curveNum < st; curveNum++)
			{
				MDataHandle outHandle = builder.addElement(curveNum);
				MFnNurbsCurveData dataCreator;
				MObject outCurveData = dataCreator.create();

				if(hasSegments)
					num_of_segments = h->GetSegmentsArray()[curveNum];

				MFnNurbsCurve newCurve;
				const unsigned deg = 3;
				const unsigned ncvs = num_of_segments + 1;
				const unsigned spans = ncvs - deg;   // Number of spans
				const unsigned nknots = spans+2*deg-1; // Number of knots

				MDoubleArray knotSequences;
				for (int i = 0; i < nknots; i++)
					knotSequences.append((double) i);

				MPointArray controlVertices;
				for (int j = 0; j < ncvs; j++)
				{
					float x = h->GetPointsArray()[currPtIndex  ];
					float y = h->GetPointsArray()[currPtIndex+2];
					float z = h->GetPointsArray()[currPtIndex+1];
					MPoint p(x, y, z);
					controlVertices.append(p);
					currPtIndex += 3;
				}
				newCurve.create(controlVertices, knotSequences, deg,
								MFnNurbsCurve::kOpen, false, false,
								outCurveData, &returnStatus);
				McheckErr(returnStatus, "Error creating curve\n");

				outHandle.set(outCurveData);
			}

			// Set the builder back into the output array.  This statement
			// is always required, no matter what constructor was used to
			// create the builder.
			returnStatus = outputArray.set(builder);
			McheckErr(returnStatus, "Error setting the builder\n");

			// Since we compute all the elements of the array, instead of
			// just marking the plug we were asked to compute as clean, mark
			// every element of the array as clean to prevent further calls
			// to this compute method during this DG evaluation cycle.
			returnStatus = outputArray.setAllClean();
			McheckErr(returnStatus, "Error cleaning outputCurves\n");
		}
	}
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
	/* THIS IS FOR CYLINDERS
	if (plug == outputMesh) {
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
			cout << "hair file is " << fileName.c_str() << endl;
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
		
	//	createCurves(newOutputData, returnStatus, *h);
		createMesh(newOutputData, returnStatus, *h);
		McheckErr(returnStatus, "ERROR creating new strand");
		
		outputHandle.set(newOutputData);
		data.setClean( plug );
	} 
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;*/
}
