/*
 * vtrFileReader.h
 *
 * Copyright (C) 2020-2021 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOL_TRIALVOLUME_VTRFILEREADER_H_INCLUDED
#define MEGAMOL_TRIALVOLUME_VTRFILEREADER_H_INCLUDED

#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/Module.h"

#include "geometry_calls/VolumetricDataCall.h"
#include <vector>

namespace megamol::trialvolume {


/**
 * Data source module for mmvtkm files
 */
class vtrFileReader : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) {
        return "vtrFileReader";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) {
        return "File loader module for vtr files.";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) {
        return true;
    }

    /** Ctor. */
    vtrFileReader(void);

    /** Dtor. */
    virtual ~vtrFileReader(void);

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    virtual bool create(void);

    /**
     * Implementation of 'Release'.
     */
    virtual void release(void);

    /**
     * Callback receiving the update of the file name parameter.
     *
     * @param slot The updated ParamSlot.
     *
     * @return Always 'true' to reset the dirty flag.
     */
    bool filenameChanged(core::param::ParamSlot& slot);

    /**
     * Gets the data from the source.
     *
     * @param caller The calling call.
     *
     * @return 'true' on success, 'false' on failure.
     */
    bool getDataCallback(core::Call& caller);

    /**
     * Gets the meta data from the source.
     *
     * @param caller The calling call.
     *
     * @return 'true' on success, 'false' on failure.
     */
    bool getMetaDataCallback(core::Call& caller);

    bool loadFile();

private:
    
    bool dummyCallback(core::Call& caller);

    /** The data update hash */
    std::size_t dataHash = 0;

    /** The file name  */
    core::param::ParamSlot filename_;

    /** The slot for requesting data */
    core::CalleeSlot getDataCalleeSlot_;

    /** The volume data */
    std::vector<float> volume;

    /** The minimum value of the volume */
    float minValue = 0.0f;

    /** The maximum value of the volume */
    float maxValue = 0.0f;

    /** The vtkm data holder metadata */
    geocalls::VolumetricDataCall::Metadata metadata;

    /** The vtr data file name */
    std::string vtrFilename;

    /** Used as flag if file has changed */
    bool fileChanged_;
    
    /** The bounding box */
    vislib::math::Cuboid<float> bbox;
};

} /* end namespace megamol */

#endif /* MEGAMOL_TRIALVOLUME_VTRFILEREADER_H_INCLUDED */
