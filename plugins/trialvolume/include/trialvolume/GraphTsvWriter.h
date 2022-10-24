#pragma once

#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "trialvolume/GraphCall.h"

namespace megamol::trialvolume {

class GraphTsvWriter : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "GraphTsvWriter";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Writes a graph to a collection of tsv files.";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

    /** Ctor. */
    GraphTsvWriter();

    /** Dtor. */
    ~GraphTsvWriter() override = default;

protected:
    bool create() override;
    void release() override;

    /** Callback for when the button is pressed
     *
     * @param slot The slot that was pressed
     *
     * @return true if the write was successful
     */
    bool writeButtonCallback(core::param::ParamSlot& slot);

    /** Callback for when the lhs is queried
     *
     * @param call The call that was made
     *
     * @return true if the write was successful
     */
    bool getDataCallback(core::Call& call);

    /** Writes the tsv files
     *
     * @return true if the write was successful
     */
    bool writeTsvFiles(GraphCall* graph_call);

private:
    /** The slot asking for data. */
    core::CallerSlot graph_query_slot_;

    /** The slot receiving data */
    core::CalleeSlot graph_receiving_slot_;

    /** The slot containing the resulting filename. */
    core::param::ParamSlot filename_slot_;
    /** The slot containing the button to force the write. */
    core::param::ParamSlot write_button_slot_;
};

} // namespace megamol::trialvolume
