#pragma once

#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

namespace megamol::trialvolume {

class GraphDotWriter : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "GraphDotWriter";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Writes a graph to a dot file.";
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
    GraphDotWriter();

    /** Dtor. */
    ~GraphDotWriter() override = default;

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

    /** Writes the dot file
     *
     * @return true if the write was successful
     */
    bool writeDotFile();

private:
    /** The slot asking for data. */
    core::CallerSlot graph_slot_;

    /** The slot containing the resulting filename. */
    core::param::ParamSlot filename_slot_;
    /** The slot containing the button to force the write. */
    core::param::ParamSlot write_button_slot_;
};

} // namespace megamol::trialvolume
