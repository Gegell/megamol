#include "SegmentationAnalysis.h"
#include "MeshSegmentationCall.h"

#include "datatools/table/TableDataCall.h"
#include "mmcore/param/ButtonParam.h"
#include "vislib/math/Cuboid.h"

#include <Eigen/Dense>

using namespace megamol::trialvolume;

SegmentationAnalysis::SegmentationAnalysis()
        : mesh_slot_("mesh", "The mesh data.")
        , tabular_output_slot_("tabular_output", "The analysis in tabular form.")
        , button_slot_("button", "The button to force analysis.") {
    // Initialize the input mesh data slot
    mesh_slot_.SetCompatibleCall<MeshSegmentationCall::segmentation_description>();
    MakeSlotAvailable(&mesh_slot_);

    // Initialize the tabular output slot
    tabular_output_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(0), &SegmentationAnalysis::analyzeSegmentsCallback);
    tabular_output_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(1), &SegmentationAnalysis::getHashCallback);
    MakeSlotAvailable(&tabular_output_slot_);

    // Initialize the output segmentation data slot
    button_slot_ << new core::param::ButtonParam();
    button_slot_.SetUpdateCallback(&SegmentationAnalysis::buttonPressedCallback);
    MakeSlotAvailable(&button_slot_);

    // Generate the table columns
    std::vector<std::string> column_names = {"ID", "Vertices", "Triangles", "CentroidX", "CentroidY", "CentroidZ",
        "ExtendsX", "ExtendsY", "ExtendsZ", "Volume", "Surface", "Sphericity", "Singular1", "Singular2", "Singular3"};

    table_columns_ = {};
    for (auto& column_name : column_names) {
        table_columns_.emplace_back();
        table_columns_.back().SetName(column_name);
        table_columns_.back().SetType(datatools::table::TableDataCall::ColumnType::QUANTITATIVE);
        table_columns_.back().SetMinimumValue(std::numeric_limits<float>::lowest());
        table_columns_.back().SetMaximumValue(std::numeric_limits<float>::max());
    }
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

bool SegmentationAnalysis::getHashCallback(core::Call& call) {
    auto& c = dynamic_cast<datatools::table::TableDataCall&>(call);
    c.SetDataHash(0);
    return true;
}

bool SegmentationAnalysis::buttonPressedCallback(core::param::ParamSlot& slot) {
    return recalculateMetrics();
}

bool SegmentationAnalysis::analyzeSegmentsCallback(core::Call& call) {
    if (!recalculateMetrics()) {
        return false;
    }

    auto* table_call = dynamic_cast<datatools::table::TableDataCall*>(&call);
    table_content_.clear();
    table_content_.reserve(segment_metrics_.size() * table_columns_.size());
    for (auto segment : segment_metrics_) {
        table_content_.push_back(static_cast<float>(segment.id));
        table_content_.push_back(static_cast<float>(segment.num_vertices));
        table_content_.push_back(static_cast<float>(segment.num_triangles));
        table_content_.push_back(segment.centroid[0]);
        table_content_.push_back(segment.centroid[1]);
        table_content_.push_back(segment.centroid[2]);
        table_content_.push_back(segment.extents[0]);
        table_content_.push_back(segment.extents[1]);
        table_content_.push_back(segment.extents[2]);
        table_content_.push_back(segment.volume);
        table_content_.push_back(segment.surface_area);
        table_content_.push_back(segment.sphericity);
        table_content_.push_back(segment.singular_vals[0]);
        table_content_.push_back(segment.singular_vals[1]);
        table_content_.push_back(segment.singular_vals[2]);
    }

    // Set the actual data call
    table_call->SetDataHash(hash_);
    table_call->Set(table_columns_.size(), segment_metrics_.size(), table_columns_.data(), table_content_.data());
    table_call->SetFrameCount(1);

    return true;
}

bool SegmentationAnalysis::recalculateMetrics() {
    // Call the segmentation to get the segments
    auto* segmentation_call = mesh_slot_.CallAs<MeshSegmentationCall>();
    if (segmentation_call == nullptr) {
        return false;
    }

    // Check whether we have to recalculate the metrics
    if (hash_ == segmentation_call->DataHash()) {
        return true;
    }
    (*segmentation_call)(0);

    // Get the segments
    auto segments = *segmentation_call->GetSegments();

    // Compute the metrics
    segment_metrics_.clear();
    for (auto& segment : segments) {
        int segment_id = segment_metrics_.size();
        segment_metrics_.push_back(computeMetrics(segment, segment_id));
    }
    hash_ = segmentation_call->DataHash();
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

    float sphericity = std::pow(36.0f * 3.1415926f * volume * volume, 1.0f / 3.0f) / surface_area;

    // Compute the singular values
    Eigen::JacobiSVD<Eigen::Matrix3Xf> svd(vertices_matrix);

    // Generate the metadata
    struct SegmentMetadata metadata = {id, segment.vertices.size(), segment.triangle_offsets.size(),
        {center[0], center[1], center[2]}, {extents[0], extents[1], extents[2]}, volume, surface_area, sphericity,
        {svd.singularValues()[0], svd.singularValues()[1], svd.singularValues()[2]}};

    // Log information
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis] Segment %d:", metadata.id);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Vertices: %d", segment.vertices.size());
    core::utility::log::Log::DefaultLog.WriteInfo(
        "[SegmentationAnalysis]   Triangles: %d", segment.triangle_offsets.size());
    core::utility::log::Log::DefaultLog.WriteInfo(
        "[SegmentationAnalysis]   Centroid: (%f, %f, %f)", center[0], center[1], center[2]);
    core::utility::log::Log::DefaultLog.WriteInfo(
        "[SegmentationAnalysis]   Extents: (%f, %f, %f)", extents[0], extents[1], extents[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Volume: %f", volume);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Surface area: %f", surface_area);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Sphericity: %f", sphericity);
    // core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Singular values: (%f, %f, %f)",
    // metadata.singular_vals[0], metadata.singular_vals[1], metadata.singular_vals[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("[SegmentationAnalysis]   Singular value ratio (1st to 3rd): %f",
        metadata.singular_vals[0] / metadata.singular_vals[2]);
    return metadata;
}