#ifndef cylinder_H_
#define cylinder_H_

#include <maya/MPoint.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MVector.h>
#include <maya/MVectorArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>

class CylinderMesh
{
public:
    CylinderMesh(const MPoint& start, const MPoint& end, double r1, double r2);
    ~CylinderMesh();

    void getMesh(
        MPointArray& points, 
        MIntArray& faceCounts, 
        MIntArray& faceConnects);
    
    void appendToMesh(
        MPointArray& points, 
        MIntArray& faceCounts, 
        MIntArray& faceConnects);

protected:
    void transform(MPointArray& points, MVectorArray& normals);
    MPoint mStart;
    MPoint mEnd;
    double r1, r2;

    // Creates a unit cylinder from (0,0,0) with radius r1 to (1,0,0) with radius r2
	static void initCylinderMesh(double r1, double r2);
    static MPointArray gPoints;
    static MVectorArray gNormals;
    static MIntArray gFaceCounts;
    static MIntArray gFaceConnects;
};

#endif