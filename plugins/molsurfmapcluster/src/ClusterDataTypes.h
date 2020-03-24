#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include "vislib/graphics/BitmapImage.h"

namespace megamol {
namespace molsurfmapcluster {

struct ClusterNode_2 {
    /** The ID of this node */
    int64_t id = -1;

    /** The ID of the left child. -1 if no left child exists */
    int64_t left = -1;

    /** The ID of the right child. -1 if no right child exists */
    int64_t right = -1;

    /** The ID of the parent node. -1 if no parent node exists */
    int64_t parent = -1;

    /** The ID of the representative node. -1 or id, if this node is the representative */
    int64_t representative = -1;

    /** The PDB id of the represented protein */
    std::string pdbID = "";

    /** Path to the represented picture file*/
    std::string picturePath = "";

    /** Pointer to the represented image */
    std::weak_ptr<vislib::graphics::BitmapImage> picturePtr;
};

/** Enum representing all available clustering methods */
enum class ClusteringMethod { IMAGE_MOMENTS = 0, COLOR_MOMENTS = 1 };

/** Enum representing all available distance measures */
enum class DistanceMeasure {
    EUCLIDEAN_DISTANCE = 0,
    L3_DISTANCE = 1,
    COSINUS_DISTANCE = 2,
    DICE_DISTANCE = 3,
    JACCARD_DISTANCE = 4
};

/** Enum representing all available linkage modes */
enum class LinkageMethod { CENTROID_LINKAGE = 0, SINGLE_LINKAGE = 1, AVERAGE_LINKAGE = 2 };

} // namespace molsurfmapcluster
} // namespace megamol
