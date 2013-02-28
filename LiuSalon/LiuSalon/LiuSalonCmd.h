#ifndef CreateLiuSalonCmd_H_
#define CreateLiuSalonCmd_H_

#include <maya/MPxCommand.h>
#include <string>
#include <maya/MSyntax.h>

class LiuSalonCmd : public MPxCommand
{
public:
    LiuSalonCmd();
    virtual ~LiuSalonCmd();
    static void* creator() { return new LiuSalonCmd(); }
    MStatus doIt( const MArgList& args );
	static MSyntax newSyntax();
};

#endif