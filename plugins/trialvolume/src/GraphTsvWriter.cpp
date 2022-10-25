#include "trialvolume/GraphTsvWriter.h"
#include "trialvolume/GraphCall.h"

#include "mmcore/param/ButtonParam.h"
#include "mmcore/param/FilePathParam.h"

#include <fstream>

using namespace megamol::trialvolume;

GraphTsvWriter::GraphTsvWriter()
        : graph_query_slot_("in_graph", "The graph to write")
        , graph_receiving_slot_("receive_graph", "The graph to write")
        , filename_slot_("filename", "The filename to write to")
        , write_button_slot_("write_button", "The button to write the graph") {
    // Initialize the input graph slot
    graph_query_slot_.SetCompatibleCall<GraphCall::graph_description>();
    MakeSlotAvailable(&graph_query_slot_);

    // Initialize the receiving graph slot
    graph_receiving_slot_.SetCallback(
        GraphCall::ClassName(), GraphCall::FunctionName(0), &GraphTsvWriter::getDataCallback);
    MakeSlotAvailable(&graph_receiving_slot_);

    // Initialize the filename slot
    filename_slot_ << new core::param::FilePathParam("", core::param::FilePathParam::Flag_Directory_ToBeCreated);
    MakeSlotAvailable(&filename_slot_);

    // Initialize the write button slot
    write_button_slot_ << new core::param::ButtonParam();
    write_button_slot_.SetUpdateCallback(&GraphTsvWriter::writeButtonCallback);
    MakeSlotAvailable(&write_button_slot_);
}

bool GraphTsvWriter::create() {
    return true;
}

void GraphTsvWriter::release() {}

bool GraphTsvWriter::writeButtonCallback(core::param::ParamSlot& slot) {
    // Unmark the button as dirty
    slot.ResetDirty();

    // Get the graph data
    auto* graph_call = graph_query_slot_.CallAs<GraphCall>();
    if (graph_call == nullptr) {
        core::utility::log::Log::DefaultLog.WriteError("GraphTsvWriter: No graph call connected to in_graph slot.");
        return false;
    }

    if (!(*graph_call)(0)) {
        core::utility::log::Log::DefaultLog.WriteError("GraphTsvWriter: Could not get graph data.");
        return false;
    }

    return writeTsvFiles(graph_call);
}

bool GraphTsvWriter::getDataCallback(core::Call& call) {
    auto* graph_call = dynamic_cast<GraphCall*>(&call);
    if (graph_call == nullptr)
        return false;

    return writeTsvFiles(graph_call);
}

bool GraphTsvWriter::writeTsvFiles(GraphCall* graph_call) {
    // Open the files
    auto dirname = filename_slot_.Param<core::param::FilePathParam>()->Value();
    if (dirname.empty()) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[%s] No directory specified for tsv files.", ClassName());
        return false;
    }

    if (!std::filesystem::exists(dirname)) {
        std::filesystem::create_directories(dirname);
    }

    std::ofstream clusters_file((dirname / "nodes.tsv").generic_u8string().c_str());
    std::ofstream connections_file((dirname / "edges.tsv").generic_u8string().c_str());

    // Get the data
    auto const& graph = graph_call->GetGraph();
    auto const& nodes = graph->nodes_;
    auto const& edges = graph->edges_;

    if (nodes.empty()) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("[%s] No nodes in graph.", ClassName());
        return true;
    } else {
        nodes.at(0)->writeTsvHeaders(clusters_file);
        clusters_file << std::endl;
        for (auto const& node : nodes) {
            node->writeTsv(clusters_file);
            clusters_file << std::endl;
        }
    }

    if (edges.empty()) {
        megamol::core::utility::log::Log::DefaultLog.WriteWarn("[%s] No edges in graph.", ClassName());
    } else {
        edges.at(0)->writeTsvHeaders(connections_file);
        connections_file << std::endl;
        for (auto const& edge : edges) {
            edge->writeTsv(connections_file);
            connections_file << std::endl;
        }
    }

    // Write meta data if available
    std::ofstream timesteps_file((dirname / "timesteps.tsv").generic_u8string().c_str());
    auto const& tsv_graph = dynamic_cast<TsvGenerator const*>(graph.get());
    if (tsv_graph != nullptr) {
        tsv_graph->writeTsvHeaders(timesteps_file);
        timesteps_file << std::endl;
        tsv_graph->writeTsv(timesteps_file);
    }

    return true;
}
