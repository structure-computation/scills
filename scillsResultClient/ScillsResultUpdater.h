#ifndef SCILLSRASULTUPDATER_H
#define SCILLSRASULTUPDATER_H

#include "DisplayScene.h"

#include <Soca/Updater.h>
class QDataStream;

/**
*/
class ScillsResultUpdater : public Updater {
protected:
    bool run( MP mp );
    virtual QString type() const { return "ScillsResultUpdater"; }
    
};

#endif // SCILLSRASULTUPDATER_H
