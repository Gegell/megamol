#pragma once

#include "datatools/AbstractParticleManipulator.h"

namespace megamol::trialvolume {

class IdentityToICol : public datatools::AbstractParticleManipulator {
public:
    /** Return module class name */
    static const char* ClassName(void) {
        return "IdentityToICol";
    }

    /** Return module class description */
    static const char* Description(void) {
        return "Copies the Identity stream of MPDC into the ICol stream";
    }

    /** Module is always available */
    static bool IsAvailable(void) {
        return true;
    }

    /** Ctor */
    IdentityToICol(void);

    /** Dtor */
    virtual ~IdentityToICol(void);

protected:
    /**
     * Manipulates the particle data
     *
     * @remarks the default implementation does not changed the data
     *
     * @param outData The call receiving the manipulated data
     * @param inData The call holding the original data
     *
     * @return True on success
     */
    virtual bool manipulateData(geocalls::MultiParticleDataCall& outData, geocalls::MultiParticleDataCall& inData);

private:
    std::vector<float> intensities;
};

} /* end namespace megamol::trialvolume */
