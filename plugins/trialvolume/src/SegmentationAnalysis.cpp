#include "SegmentationAnalysis.h"
#include "MeshSegmentationCall.h"

#include "mmcore/param/ButtonParam.h"
#include "vislib/math/Cuboid.h"

#include <Eigen/Dense>

using namespace megamol::trialvolume;

SegmentationAnalysis::SegmentationAnalysis()
        : mesh_slot_("mesh", "The mesh data.")
        , button_slot_("button", "The button to force analysis.") {
    // Initialize the input mesh data slot
    mesh_slot_.SetCompatibleCall<MeshSegmentationCall::segmentation_description>();
    MakeSlotAvailable(&mesh_slot_);

    // Initialize the output segmentation data slot
    button_slot_ << new core::param::ButtonParam();
    button_slot_.SetUpdateCallback(&SegmentationAnalysis::buttonPressedCallback);
    MakeSlotAvailable(&button_slot_);
}

SegmentationAnalysis::~SegmentationAnalysis() {
    release();
}

bool SegmentationAnalysis::create(void) {
    return true;
}

void SegmentationAnalysis::release(void) {
    // TODO release any bound resources here
}

bool SegmentationAnalysis::buttonPressedCallback(core::param::ParamSlot& slot) {
    // Call the segmentation to get the segments
    auto* segmentation_call = mesh_slot_.CallAs<MeshSegmentationCall>();
    if (segmentation_call == nullptr) {
        return false;
    }
    (*segmentation_call)(0);

    // Get the segments
    auto segments = *segmentation_call->GetSegments();

    // Compute the metrics
    int id = 0;
    for (auto& segment : segments) {
        computeMetrics(segment, id++);
    }
    return true;
}
using Eigen::Matrix3Xf;
using Eigen::Vector3f;

inline float signedVolumeOfTriangle(Vector3f p1, Vector3f p2, Vector3f p3) {
    return p1.dot(p2.cross(p3)) / 6.0f;
}

inline float surfaceAreaOfTriangle(Vector3f p1, Vector3f p2, Vector3f p3) {
    return (p2 - p1).cross(p3 - p1).norm() / 2.0f;
}

SegmentationAnalysis::SegmentMetadata SegmentationAnalysis::computeMetrics(
    const MeshSegmentation::Segment& segment, const int id) {
    auto base_vertices = *segment.base_vertices;
    auto base_indices = *segment.base_indices;

    typedef Eigen::Map<Vector3f> Vector3fMap;

    Matrix3Xf vertices_matrix(3, segment.vertices.size());
    for (size_t i = 0; i < segment.vertices.size(); i++) {
        auto offset = segment.vertices[i] * 3;
        vertices_matrix.col(i) = Vector3fMap(&base_vertices[offset]);
    }

    Vector3f min = vertices_matrix.rowwise().minCoeff();
    Vector3f max = vertices_matrix.rowwise().maxCoeff();
    Vector3f center = (min + max) / 2.0f;
    Vector3f extents = max - min;

    vertices_matrix = (vertices_matrix.colwise() - center).eval();

    float volume = 0.0f;
    float surface_area = 0.0f;

    for (auto tri : segment.triangle_offsets) {
        auto v0 = Vector3fMap(&base_vertices[base_indices[tri + 0] * 3]) - center;
        auto v1 = Vector3fMap(&base_vertices[base_indices[tri + 1] * 3]) - center;
        auto v2 = Vector3fMap(&base_vertices[base_indices[tri + 2] * 3]) - center;

        volume += signedVolumeOfTriangle(v0, v1, v2);
        surface_area += surfaceAreaOfTriangle(v0, v1, v2);
    }

    // Compute the singular values
    Eigen::JacobiSVD<Eigen::Matrix3Xf> svd(vertices_matrix);

    // Generate the metadata
    struct SegmentMetadata metadata = {id, {center[0], center[1], center[2]},
        vislib::math::Cuboid(min[0], min[1], min[2], max[0], max[1], max[2]), volume, surface_area,
        {svd.singularValues()[0], svd.singularValues()[1], svd.singularValues()[2]}};

    // Log information
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis] Segment %d:", metadata.id);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Vertices: %d",
        segment.vertices.size());
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Triangles: %d",
        segment.triangle_offsets.size());
    core::utility::log::Log::DefaultLog.WriteInfo(
        "[SegmentationAnalysis]   Centroid: (%f, %f, %f)", center[0], center[1], center[2]);
    core::utility::log::Log::DefaultLog.WriteInfo(
        "[SegmentationAnalysis]   Extents: (%f, %f, %f)", extents[0], extents[1], extents[2]);
    // core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Volume: %f", volume);
    // core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Surface area: %f", surface_area);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Volume to surface area: %f",
        volume / surface_area);
    // core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Singular values: (%f, %f, %f)",
        // metadata.singular_vals[0], metadata.singular_vals[1], metadata.singular_vals[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Singular value ratio (1st to 3rd): %f",
        metadata.singular_vals[0] / metadata.singular_vals[2]);
    return metadata;
}