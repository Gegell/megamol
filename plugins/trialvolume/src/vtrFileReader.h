/*
 * VtrFileReader.h
 *
 * Copyright (C) 2020-2021 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#pragma once

#include <vector>

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

namespace megamol::trialvolume {


/**
 * Data source module for mmvtkm files
 */
class VtrFileReader : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) {
        return "VtrFileReader";
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
    VtrFileReader(void);

    /** Dtor. */
    virtual ~VtrFileReader(void);

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
    std::size_t data_hash_ = 0;

    /** The file name  */
    core::param::ParamSlot filename_;

    /** The slot for requesting data */
    core::CalleeSlot get_data_callee_slot_;

    /** The volume data */
    std::vector<float> volume_;

    /** The minimum value of the volume */
    float min_value_ = 0.0f;

    /** The maximum value of the volume */
    float max_value_ = 0.0f;

    /** The vtkm data holder metadata */
    geocalls::VolumetricDataCall::Metadata metadata_;

    /** The vtr data file name */
    std::string vtr_filename_;

    /** Used as flag if file has changed */
    bool file_changed_;

    /** The bounding box */
    vislib::math::Cuboid<float> bbox_;
};

} // namespace megamol::trialvolume
