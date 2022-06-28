#include "SegmentationAnalysis.h"
#include "MeshSegmentationCall.h"

#include "datatools/table/TableDataCall.h"
#include "mmcore/param/ButtonParam.h"
#include "mmcore/param/TransferFunctionParam.h"
#include "vislib/math/Cuboid.h"
#include "mesh/MeshDataCall.h"

#include <Eigen/Dense>

using namespace megamol::trialvolume;

SegmentationAnalysis::SegmentationAnalysis()
        : mesh_slot_("mesh", "The mesh data.")
        , mesh_data_slot_("mesh_data", "The analysis data associated with the mesh.")
        , tabular_output_slot_("tabular_output", "The analysis in tabular form.")
        , transfer_function_slot_("transfer_function", "The transfer function to use for the datasets.")
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

    // Initialize the transfer function slot
    transfer_function_slot_ << new core::param::TransferFunctionParam();
    transfer_function_slot_.SetUpdateCallback(&SegmentationAnalysis::transferFunctionCallback);
    MakeSlotAvailable(&transfer_function_slot_);

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

    // Initialize the mesh data slot
    mesh_data_slot_.SetCallback(mesh::MeshDataCall::ClassName(), mesh::MeshDataCall::FunctionName(0),
        &SegmentationAnalysis::meshDataCallCallback);
    mesh_data_slot_.SetCallback(mesh::MeshDataCall::ClassName(), mesh::MeshDataCall::FunctionName(1),
        &SegmentationAnalysis::meshMetadataCallCallback);
    MakeSlotAvailable(&mesh_data_slot_);

    // Initialize the output data sets
    for (int i = 0; i < output_data_sets_.size(); i++) {
        output_data_sets_[i] = std::make_shared<mesh::MeshDataCall::data_set>();
        output_data_sets_[i]->data = std::make_shared<std::vector<float>>();
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
    c.SetDataHash(hash_);
    return true;
}

bool SegmentationAnalysis::transferFunctionCallback(core::param::ParamSlot& param) {
    param.ResetDirty();
    for (auto i = 0; i < output_data_sets_.size(); i++) {
        output_data_sets_[i]->transfer_function = param.Param<core::param::TransferFunctionParam>()->Value();
        output_data_sets_[i]->transfer_function_dirty = true;
    }
    return true;
}

bool SegmentationAnalysis::buttonPressedCallback(core::param::ParamSlot& slot) {
    return recalculateMetrics();
}

bool SegmentationAnalysis::analyzeSegmentsCallback(core::Call& call) {
    auto* table_call = dynamic_cast<datatools::table::TableDataCall*>(&call);
    if (table_call == nullptr) {
        return false;
    }
    if (table_call->DataHash() == hash_) {
        return true;
    }

    if (!ensureFreshMetrics()) {
        return false;
    }

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

bool SegmentationAnalysis::meshMetadataCallCallback(core::Call& _call) {
    auto* call = dynamic_cast<mesh::MeshDataCall*>(&_call);
    if (call == nullptr) {
        return false;
    }

    // The following data is computed:
    // id                           - Unique ID of the segment
    call->set_data("id");
    // surface                      - Surface area of the segment
    call->set_data("surface");
    // volume                       - Volume of the segment
    call->set_data("volume");
    // sphericity                   - Sphericity of the segment
    call->set_data("sphericity");
    // singular_vals_largest_to_smallest - Ratio of largest to smallest singular value of the segment
    call->set_data("singular_vals_largest_to_smallest");

    return true;
}


bool SegmentationAnalysis::meshDataCallCallback(core::Call& call) {
    auto* mdc = dynamic_cast<mesh::MeshDataCall*>(&call);

    if (mdc == nullptr) {
        return false;
    }
    if (mdc->DataHash() == hash_) {
        return true;
    }
    if (!ensureFreshMetrics()) {
        return false;
    }

    mdc->set_data("id", output_data_sets_[0]);
    mdc->set_data("surface", output_data_sets_[1]);
    mdc->set_data("volume", output_data_sets_[2]);
    mdc->set_data("sphericity", output_data_sets_[3]);
    mdc->set_data("singular_vals_smallest_to_largest", output_data_sets_[4]);

    for (auto i = 0; i < output_data_sets_.size(); i++) {
        output_data_sets_[i]->data->resize(input_data_.vertices->size() / 3, -1);
        output_data_sets_[i]->min_value = 0;
        output_data_sets_[i]->max_value = 0;
    }

    for (auto metrics : segment_metrics_) {
        auto segment = input_data_.segments->at(metrics.id);
        for (auto vert : segment.vertices) {
            (*output_data_sets_[0]->data)[vert] = metrics.id;
            (*output_data_sets_[1]->data)[vert] = metrics.surface_area;
            (*output_data_sets_[2]->data)[vert] = metrics.volume;
            (*output_data_sets_[3]->data)[vert] = metrics.sphericity;
            (*output_data_sets_[4]->data)[vert] = metrics.singular_vals[2] / metrics.singular_vals[0];
        }
        output_data_sets_[1]->max_value = std::max(output_data_sets_[1]->max_value, metrics.surface_area);
        output_data_sets_[2]->max_value = std::max(output_data_sets_[2]->max_value, metrics.volume);
        output_data_sets_[3]->max_value = std::max(output_data_sets_[3]->max_value, metrics.sphericity);
        output_data_sets_[4]->max_value = std::max(output_data_sets_[4]->max_value, metrics.singular_vals[2] / metrics.singular_vals[0]);
    }
    output_data_sets_[0]->max_value = segment_metrics_.size() - 1;

    mdc->SetDataHash(hash_);

    return true;
}

bool SegmentationAnalysis::ensureFreshMetrics() {
    return refreshInput() && recalculateMetrics();
}

bool SegmentationAnalysis::refreshInput() {
    // Call the segmentation to get the segments
    auto* segmentation_call = mesh_slot_.CallAs<MeshSegmentationCall>();
    if (segmentation_call == nullptr) {
        return false;
    }

    // Check whether we have to recalculate the metrics
    if (segmentation_call->DataHash() == input_data_.hash) {
        return true;
    }
    (*segmentation_call)(0);

    // Get the segments
    input_data_.segments = segmentation_call->GetSegments();
    input_data_.vertices = segmentation_call->GetBaseVertices();
    input_data_.indices = segmentation_call->GetBaseIndices();
    input_data_.hash = segmentation_call->DataHash();
    return true;
}

bool SegmentationAnalysis::recalculateMetrics() {
    // Exit early if the hash is the same as the last time we called this function
    if (input_data_.hash == hash_) {
        return true;
    }

    // Compute the metrics
    segment_metrics_.clear();
    for (auto& segment : *input_data_.segments) {
        int segment_id = segment_metrics_.size();
        segment_metrics_.push_back(computeMetrics(segment, segment_id, *input_data_.vertices, *input_data_.indices));
    }

    hash_ = input_data_.hash;
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
    const MeshSegmentation::Segment& segment, const int id, const std::vector<float>& vertices,
    const std::vector<unsigned int>& indices) {

    typedef Eigen::Map<const Vector3f> Vector3fMap;

    Matrix3Xf vertices_matrix(3, segment.vertices.size());
    for (size_t i = 0; i < segment.vertices.size(); i++) {
        auto offset = segment.vertices[i] * 3;
        vertices_matrix.col(i) = Vector3fMap(&vertices[offset]);
    }

    Vector3f min = vertices_matrix.rowwise().minCoeff();
    Vector3f max = vertices_matrix.rowwise().maxCoeff();
    Vector3f center = (min + max) / 2.0f;
    Vector3f extents = max - min;

    vertices_matrix = (vertices_matrix.colwise() - center).eval();

    float volume = 0.0f;
    float surface_area = 0.0f;

    for (auto tri : segment.triangle_offsets) {
        auto v0 = Vector3fMap(&vertices[indices[tri + 0] * 3]) - center;
        auto v1 = Vector3fMap(&vertices[indices[tri + 1] * 3]) - center;
        auto v2 = Vector3fMap(&vertices[indices[tri + 2] * 3]) - center;

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
#ifdef TRIALVOLUME_VERBOSE
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
#endif
    return metadata;
}