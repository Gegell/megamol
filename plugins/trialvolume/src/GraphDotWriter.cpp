#include "trialvolume/GraphDotWriter.h"
#include "trialvolume/GraphCall.h"

#include "mmcore/param/ButtonParam.h"
#include "mmcore/param/FilePathParam.h"

#include <fstream>

using namespace megamol::trialvolume;

GraphDotWriter::GraphDotWriter()
        : graph_slot_("in_graph", "The graph to write")
        , filename_slot_("filename", "The filename to write to")
        , write_button_slot_("write_button", "The button to write the graph") {
    // Initialize the input graph slot
    graph_slot_.SetCompatibleCall<GraphCall::graph_description>();
    MakeSlotAvailable(&graph_slot_);

    // Initialize the filename slot
    filename_slot_ << new core::param::FilePathParam("", core::param::FilePathParam::Flag_File_ToBeCreated);
    MakeSlotAvailable(&filename_slot_);

    // Initialize the write button slot
    write_button_slot_ << new core::param::ButtonParam();
    write_button_slot_.SetUpdateCallback(&GraphDotWriter::writeButtonCallback);
    MakeSlotAvailable(&write_button_slot_);
}

bool GraphDotWriter::create() {
    return true;
}

void GraphDotWriter::release() {}

bool GraphDotWriter::writeButtonCallback(core::param::ParamSlot& slot) {
    return writeDotFile();
}

bool GraphDotWriter::writeDotFile() {
    // Get the graph data
    auto* graph_call = graph_slot_.CallAs<GraphCall>();
    if (graph_call == nullptr) {
        core::utility::log::Log::DefaultLog.WriteError("GraphDotWriter: No graph call connected to in_graph slot.");
        return false;
    }

    if (!(*graph_call)(0)) {
        core::utility::log::Log::DefaultLog.WriteError("GraphDotWriter: Could not get graph data.");
        return false;
    }

    // Open the file
    auto filename = filename_slot_.Param<core::param::FilePathParam>()->Value();
    if (filename.empty()) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[ParticleClusterTracking] No filename specified for dot file.");
        return false;
    }
    std::ofstream dot_file(filename.generic_u8string().c_str());
    if (!dot_file.is_open()) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[ParticleClusterTracking] Could not open dot file \"%s\" for writing.",
            filename.generic_u8string().c_str());
        return false;
    }

    // Write the graph
    dot_file << "digraph G {" << std::endl;

    // Write the nodes
    for (auto& node : *graph_call->GetNodes()) {
        dot_file << node.ToDot() << std::endl;
    }

    // Write the edges
    for (auto& edge : *graph_call->GetEdges()) {
        dot_file << edge.ToDot() << std::endl;
    }

    dot_file << "}" << std::endl;
    return true;
}
