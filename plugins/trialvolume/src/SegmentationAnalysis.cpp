#include "SegmentationAnalysis.h"
#include "MeshSegmentationCall.h"

#include "mmcore/param/ButtonParam.h"
#include "vislib/math/Cuboid.h"

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
        auto metadata = SegmentMetadata();
        metadata.id = id;
        computeMetrics(segment, metadata);
    }
    return true;
}

using Eigen::Matrix3Xf;
using Eigen::Vector3f;

bool SegmentationAnalysis::computeMetrics(const MeshSegmentation::Segment& segment, SegmentMetadata& metadata) {
    auto base_vertices = *segment.base_vertices;
    auto base_indices = *segment.base_indices;

    typedef Eigen::Map<Vector3f> Vector3fMap;

    Matrix3Xf vertices_matrix(3, segment.vertices.size());
    for (size_t i = 0; i < segment.vertices.size(); i++) {
        auto offset = segment.vertices[i] * 3;
        // vertices_matrix.col(i) = Vector3fMap(base_vertices.data() + offset);
        vertices_matrix.col(i) = Vector3f(base_vertices[offset], base_vertices[offset + 1], base_vertices[offset + 2]);
    }

    Vector3f min = vertices_matrix.rowwise().minCoeff();
    Vector3f max = vertices_matrix.rowwise().maxCoeff();
    Vector3f center = (min + max) / 2.0f;
    Vector3f extents = max - min;

    vertices_matrix = (vertices_matrix.colwise() - center).eval();

    float volume = 0.0f;
    float surface_area = 0.0f;

    // for (auto tri : segment.triangle_offsets) {
    //     auto v0 = vertices_matrix.col(tri / 3);
    //     auto v1 = vertices_matrix.col(tri / 3 + 1);
    //     auto v2 = vertices_matrix.col(tri / 3 + 2);

    //     auto n0 = v0 - center;
    //     auto n1 = v1 - center;
    //     auto n2 = v2 - center;

    //     volume += signedVolumeOfTriangle(n0, n1, n2);
    //     surface_area += surfaceAreaOfTriangle(n0, n1, n2);
    // }

    // Generate the metadata
    metadata.centroid = center;
    metadata.bounds = vislib::math::Cuboid<float>(min[0], min[1], min[2], max[0], max[1], max[2]);

    metadata.volume = volume;
    metadata.surface_area = surface_area;

    // Compute the singular values
    Eigen::JacobiSVD<Eigen::Matrix3Xf> svd(vertices_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    metadata.singular_vals << svd.singularValues()[0], svd.singularValues()[1], svd.singularValues()[2];

    // Log information
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis] Segment %d:", metadata.id);
    core::utility::log::Log::DefaultLog.WriteInfo(
        "[SegmentationAnalysis]   Centroid: (%f, %f, %f)", center[0], center[1], center[2]);
    core::utility::log::Log::DefaultLog.WriteInfo(
        "[SegmentationAnalysis]   Extents: (%f, %f, %f)", extents[0], extents[1], extents[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Volume: %f", volume);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Surface area: %f", surface_area);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Singular values: (%f, %f, %f)",
        metadata.singular_vals[0], metadata.singular_vals[1], metadata.singular_vals[2]);
    return true;
}

inline float SegmentationAnalysis::signedVolumeOfTriangle(Vector3f p1, Vector3f p2, Vector3f p3) {
    return p1.dot(p2.cross(p3)) / 6.0f;
}

inline float SegmentationAnalysis::surfaceAreaOfTriangle(Vector3f p1, Vector3f p2, Vector3f p3) {
    return (p2 - p1).cross(p3 - p1).norm() / 2.0f;
}