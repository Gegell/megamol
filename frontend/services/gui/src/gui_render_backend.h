/*
 * gui_render_backend.h
 *
 * Copyright (C) 2021 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOL_GUI_GUIRENDERBACKEND_H_INCLUDED
#define MEGAMOL_GUI_GUIRENDERBACKEND_H_INCLUDED
#pragma once


#include "ImageWrapper.h"
#include "imgui.h"
#include "imgui_impl_generic.h"
#include "imgui_sw.h"
#include "mmcore/view/CPUFramebuffer.h"
#include <glm/glm.hpp>
#include <memory>

#ifdef WITH_GL
#include <glowl/FramebufferObject.hpp>
#endif // WITH_GL


namespace megamol {
namespace gui {


enum class GUIRenderBackend { NONE, CPU, OPEN_GL };


/** ************************************************************************
 * Managing available ImGui backends
 */
class gui_render_backend {
public:
    /**
     * CTOR.
     */
    gui_render_backend();

    /**
     * DTOR.
     */
    ~gui_render_backend();

    bool IsBackendInitialized() {
        return (this->initialized_backend != GUIRenderBackend::NONE);
    }

    bool CheckPrerequisites(GUIRenderBackend backend);

    bool Init(GUIRenderBackend backend);

    void ClearFrame();

    void NewFrame(glm::vec2 framebuffer_size, glm::vec2 window_size);

    bool EnableRendering(unsigned int width, unsigned int height);

    bool Render(ImDrawData* draw_data);

    bool ShutdownBackend();

    bool CreateFontsTexture();

    megamol::frontend_resources::ImageWrapper GetImage();

private:
    // VARIABLES --------------------------------------------------------------

    GUIRenderBackend initialized_backend;

    GenericWindow cpu_window;
    GenericMonitor cpu_monitor;

#ifdef WITH_GL
    std::shared_ptr<glowl::FramebufferObject> ogl_fbo = nullptr;
#endif // WITH_GL
    std::shared_ptr<megamol::core::view::CPUFramebuffer> cpu_fbo = nullptr;

    // FUNCTIONS --------------------------------------------------------------

    bool createCPUFramebuffer(unsigned int width, unsigned int height);
    bool createOGLFramebuffer(unsigned int width, unsigned int height);
};


} // namespace gui
} // namespace megamol

#endif // MEGAMOL_GUI_GUIRENDERBACKEND_H_INCLUDED
