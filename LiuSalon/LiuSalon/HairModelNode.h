#include <maya/MTime.h>
#include <maya/MFnMesh.h>
#include <maya/MPoint.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MPointArray.h>
#include <maya/MPxNode.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MFnMeshData.h>
#include "cyHairFile.h"
#include "cylinder.h"

class HairModelNode : public MPxNode
{
public:
					HairModelNode() {};
	virtual 		~HairModelNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static  void*	creator();
	static  MStatus initialize();

	static MObject inputCurve;
	static MObject outputMesh;
	static MObject numStrands;
	static MObject numPoints;
	static MObject hairLength;
	static MObject file;
	static MTypeId	id;

protected:
	MObject createMesh(MObject& outData, MStatus& stat, cyHairFile& hair);
	MStatus createHairCurve(MArrayDataHandle &inputArray, cyHairFile& hair);
};