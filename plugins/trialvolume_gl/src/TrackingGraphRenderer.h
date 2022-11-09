/**
 * MegaMol
 * Copyright (c) 2018, MegaMol Dev Team
 * All rights reserved.
 */

#pragma once

#include <memory>

#include <glowl/glowl.h>

#include "mmcore/CallerSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmstd_gl/renderer/Renderer3DModuleGL.h"

namespace megamol::trialvolume_gl {

/**
 * Renders incoming tracking graph to screen.
 */
class TrackingGraphRenderer : public mmstd_gl::Renderer3DModuleGL {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "TrackingGraphRenderer";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Renders the connections of the incoming tracking graph";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

    /** Constructor. */
    TrackingGraphRenderer();

    /** Destructor. */
    ~TrackingGraphRenderer() override;

protected:
    /**
     * Implementation of 'Create'.
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() override;

    /**
     * Implementation of 'Release'.
     */
    void release() override;

private:
    /**
     * The get extents callback. The module should set the members of
     * 'call' to tell the caller the extents of its data (bounding boxes
     * and times).
     *
     * TUTORIAL: This method computes the extents of the rendered data, namely the bounding box and other relevant
     *  values, and writes them into the calling call.
     *
     * @param call The calling call.
     * @return The return value of the function.
     */
    bool GetExtents(mmstd_gl::CallRender3DGL& call) override;

    /**
     * The Open GL Render callback.
     *
     * TUTORIAL: Mandatory method for each renderer. It performs the main OpenGL rendering.
     *
     * @param call The calling call.
     * @return The return value of the function.
     */
    bool Render(mmstd_gl::CallRender3DGL& call) override;

    /** The input data slot. */
    core::CallerSlot in_graph_data_slot_;

    /** The vertex buffer object for the rendered vertices. */
    GLuint vbo;

    /** The vertex array for the rendered vertices. */
    GLuint va;

    /** The index buffer object for the rendered vertices. */
    GLuint ibo;

    /** The number of line indices to render */
    size_t num_indices_;

    /** The data hash of the most recent rendered data */
    std::size_t last_data_hash_;

    /** The simple shader to draw the lines between the nodes */
    std::unique_ptr<glowl::GLSLProgram> line_shader_;

    /** Slot for the line shader enable */
    core::param::ParamSlot draw_connections_slot_;

    /** The shader to draw the bboxes of the clusters in the current timestep*/
    std::unique_ptr<glowl::GLSLProgram> bbox_shader_;

    /** Slot for the bbox shader enable */
    core::param::ParamSlot draw_bboxes_slot_;

    /** Slot for the scaling factor of the line width*/
    core::param::ParamSlot line_width_slot_;

    /** Bounding box */
    vislib::math::Cuboid<float> bbox_;

    /** Last frame encountered */
    int last_frame_;
};

} // namespace megamol::trialvolume_gl
