#include "MeshSegmentation.h"
#include "MeshSegmentationCall.h"

#include "mesh/TriangleMeshCall.h"
#include "mmcore/param/ButtonParam.h"

#include <vector>
#include <map>

using namespace megamol::trialvolume;

MeshSegmentation::MeshSegmentation()
        : mesh_data_slot_("mesh_data", "The mesh data to segment")
        , segmentation_data_slot_("segmentation_data", "The segmentation data")
        , force_update_slot_("force_segmentation", "Force segmentation") {
    // Initialize the input mesh data slot
    mesh_data_slot_.SetCompatibleCall<megamol::mesh::TriangleMeshCall::triangle_mesh_description>();
    MakeSlotAvailable(&mesh_data_slot_);

    // Initialize the output segmentation data slot
    segmentation_data_slot_.SetCallback(MeshSegmentationCall::ClassName(),
        MeshSegmentationCall::FunctionName(0), &MeshSegmentation::getSegmentationCallback);
    MakeSlotAvailable(&segmentation_data_slot_);

    // Initialize the force update slot
    force_update_slot_ << new core::param::ButtonParam();
    force_update_slot_.SetUpdateCallback(&MeshSegmentation::buttonPressedCallback);
    MakeSlotAvailable(&force_update_slot_);
}

MeshSegmentation::~MeshSegmentation() {
    release();
}

bool MeshSegmentation::create() {
    return true;
}

void MeshSegmentation::release() {
    // TODO release any bound resources here
}

bool MeshSegmentation::buttonPressedCallback(core::param::ParamSlot& slot) {
    auto dummy_call = new MeshSegmentationCall();
    getSegmentationCallback(*dummy_call);
    core::utility::log::Log::DefaultLog.WriteInfo("Segmentation updated");
    return true;
}

int merge(std::map<size_t, size_t> &parent, size_t index) {
    if (parent.find(index) == parent.end()) {
        parent[index] = index;
        return index;
    }
    auto root = index;
    while (parent[root] != root) {
        root = parent[root];
    }
    while (parent[index] != root) {
        auto next = parent[index];
        parent[index] = root;
        index = next;
    }
    return root;
}

bool MeshSegmentation::getSegmentationCallback(core::Call& call) {
    // Get the mesh data
    auto* meshCall = mesh_data_slot_.CallAs<megamol::mesh::TriangleMeshCall>();
    if (meshCall == nullptr) {
        return false;
    }
    (*meshCall)(0);

    // Generate the segmentation data call
    auto segmentationCall = dynamic_cast<MeshSegmentationCall*>(&call);
    if (segmentationCall == nullptr) {
        return false;
    }

    if (meshCall->DataHash() == hash_) {
        return true;
    }
    hash_ = meshCall->DataHash();

    auto indices = *(meshCall->get_indices());

    // Perform union-find algorithm to find connected components

    // Store the parent for each vertex
    // TODO make it so that it could also use the edges as connections, not just vertices
    auto parent = std::map<size_t, size_t>();
    // Iterate over every triangle
    for (size_t i = 0; i < indices.size(); i += 3) {
        // Get the indices of the three vertices of the triangle
        size_t a = indices[i];
        size_t b = indices[i + 1];
        size_t c = indices[i + 2];
        // Find the parent of each vertex
        parent[merge(parent, b)] = merge(parent, a);
        parent[merge(parent, c)] = merge(parent, a);
    }
    // Extract the connected components
    segments_ = std::make_shared<std::vector<Segment>>();
    auto parent_to_segment_index = std::map<size_t, size_t>();
    for (size_t i = 0; i < parent.size(); i++) {
        parent[indices[i]] = merge(parent, indices[i]);
        if (parent_to_segment_index.find(parent[indices[i]]) == parent_to_segment_index.end()) {
            parent_to_segment_index[parent[indices[i]]] = segments_->size();
            segments_->push_back(Segment());
        }
        auto segment = &(*segments_)[parent_to_segment_index[parent[indices[i]]]];
        segment->vertices.push_back(indices[i]);
    }
    for (size_t i = 0; i < indices.size(); i += 3) {
        auto segment = &(*segments_)[parent_to_segment_index[merge(parent, indices[i])]];
        segment->triangle_offsets.push_back(i);
    }

    // Populate the segmentation data call
    segmentationCall->SetSegmentation(segments_);
    segmentationCall->SetDataHash(hash_);

    // Print some information about the segmentation
    core::utility::log::Log::DefaultLog.WriteInfo("Segmentation has %d segments", segments_->size());
    for (size_t i = 0; i < segments_->size(); i++) {
        core::utility::log::Log::DefaultLog.WriteInfo("Segment %d has %d vertices and %d triangles",
            i, (*segments_)[i].vertices.size(), (*segments_)[i].triangle_offsets.size());
    }

    return true;
}
