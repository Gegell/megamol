#pragma once

#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/Module.h"

#include <vector>

namespace megamol::trialvolume {

class MeshSegmentation : public core::Module {
public:
    struct Segment {
        std::shared_ptr<std::vector<float>> base_vertices;
        std::shared_ptr<std::vector<float>> base_normals;
        std::shared_ptr<std::vector<unsigned int>> base_indices;

        std::vector<size_t> vertices;
        std::vector<size_t> triangle_offsets;
    };

    /**
     * Answer the name of the objects of this description.
     */
    static const char* ClassName() {
        return "MeshSegmentation";
    }

    /**
     * Gets a human readable description of the module.
     */
    static const char* Description() {
        return "Segments a mesh by loose parts.";
    }

    /**
     * Answers whether this module is available on the current system.
     */
    static bool IsAvailable() {
        return true;
    }

    /** Constructor */
    MeshSegmentation();

    /** Destructor */
    ~MeshSegmentation() override;

private:
    /** Gets called every time the module is created. */
    bool create() override;

    /** Gets called every time the module is released. */
    void release() override;

    /**
     * Gets the segments found in the mesh.
     */
    bool getSegmentationCallback(core::Call& call);

    /**
     * Triggers the segmentation manually.
     */
    bool buttonPressedCallback(core::param::ParamSlot& slot);

    /** The slot for the mesh data. */
    core::CallerSlot mesh_data_slot_;

    /** The slot for the segmentation data. */
    core::CalleeSlot segmentation_data_slot_;

    /** Force update button. */
    core::param::ParamSlot force_update_slot_;

    /** Last found segments. */
    std::shared_ptr<std::vector<Segment>> segments_;

    /** Hash for the last calculated segments. */
    size_t hash_;
};

} // namespace megamol::trialvolume
