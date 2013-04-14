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

class HairShaderNode : public MPxNode
{
public:
					HairShaderNode() {};
	virtual 		~HairShaderNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static  void*	creator();
	static  MStatus initialize();

	// primary highlight attributes (R)
	static MObject inColorR;
	static MObject outColorR;
	static MObject inIntensityR;
	static MObject outIntensityR;
	static MObject inLongShiftR;
	static MObject outLongShiftR;
	static MObject inLongWidthR;
	static MObject outLongWidthR;

	// backlit light attributes (TT)
	static MObject inColorTT;
	static MObject outColorTT;
	static MObject inIntensityTT;
	static MObject outIntensityTT;
	static MObject inLongShiftTT;
	static MObject outLongShiftTT;
	static MObject inLongWidthTT;
	static MObject outLongWidthTT;
	static MObject inAzimuthWidthTT;
	static MObject outAzimuthWidthTT;

	// secondary highlight attributes (TRT-G)
	static MObject inColorTRT; // this will affect glints as well
	static MObject outColorTRT;
	static MObject inIntensityTRT;
	static MObject outIntensityTRT;
	static MObject inLongShiftTRT;
	static MObject outLongShiftTRT;
	static MObject inLongWidthTRT;
	static MObject outLongWidthTRT;

	// glint attributes (G)
	static MObject inIntensityG;
	static MObject outIntensityG;
	static MObject inAzimuthWidthG; // frequency of glints
	static MObject outAzimuthWidthG;

	static MTypeId	id;

protected:
	MObject createMesh(MObject& outData, MStatus& stat, cyHairFile& hair);
	MStatus createHairCurve(MArrayDataHandle &inputArray, cyHairFile& hair);
};