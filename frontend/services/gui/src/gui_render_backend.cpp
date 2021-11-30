/*
 * gui_render_backend.cpp
 *
 * Copyright (C) 2021 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#include "gui_render_backend.h"
#include "mmcore/utility/log/Log.h"

#ifdef WITH_GL
#include "glad/glx.h"
#include "imgui_impl_opengl3.h"
#endif

using namespace megamol::gui;


gui_render_backend::gui_render_backend()
        : initialized_backend(GUIRenderBackend::NONE)
        , sw_window({1280, 720, 0, 0})
        , sw_monitor({1920, 1080}) {}


gui_render_backend::~gui_render_backend() {
#ifdef WITH_GL
    this->ogl_fbo.reset();
#endif // WITH_GL
    this->cpu_fbo.reset();
}


bool megamol::gui::gui_render_backend::CheckPrerequisites(GUIRenderBackend backend) {

    switch (backend) {
    case (GUIRenderBackend::OPEN_GL): {
#ifdef WITH_GL
        bool prerequisities_given = true;

#ifdef _WIN32 // Windows
        HDC ogl_current_display = ::wglGetCurrentDC();
        HGLRC ogl_current_context = ::wglGetCurrentContext();
        if (ogl_current_display == nullptr) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] There is no OpenGL rendering context available.");
            prerequisities_given = false;
        }
        if (ogl_current_context == nullptr) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] There is no current OpenGL rendering context available from the calling thread.");
            prerequisities_given = false;
        }
#else
        // LINUX
        /// XXX The following throws segfault if OpenGL is not loaded yet:
        // Display* gl_current_display = ::glXGetCurrentDisplay();
        // GLXContext ogl_current_context = ::glXGetCurrentContext();
        /// XXX Is there a better way to check existing OpenGL context?
        if (glXGetCurrentDisplay == nullptr) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] There is no OpenGL rendering context available.");
            prerequisities_given = false;
        }
        if (glXGetCurrentContext == nullptr) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] There is no current OpenGL rendering context available from the calling thread.");
            prerequisities_given = false;
        }
#endif // _WIN32

        if (!prerequisities_given) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] Missing prerequisities to initialize render backend OpenGL. [%s, %s, line %d]", __FILE__,
                __FUNCTION__, __LINE__);
            return false;
        }
#else
        megamol::core::utility::log::Log::DefaultLog.WriteError("[GUI] Render backend OpenGL is not available.");
        return false;
#endif // WITH_GL
    } break;
    case (GUIRenderBackend::CPU): {
        /// nothing to check for CPU rendering ...
    } break;
    default: {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Unknown render backend... [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    } break;
    }

    return true;
}


bool gui_render_backend::Init(GUIRenderBackend backend) {

    if (this->IsBackendInitialized()) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Render backend is already initialized. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }

    switch (backend) {
    case (GUIRenderBackend::OPEN_GL): {
#ifdef WITH_GL
        if (ImGui_ImplOpenGL3_Init(nullptr)) { // "#version 130"
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("[GUI] Created ImGui render backend for Open GL.");
        } else {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] Unable to initialize OpenGL render backend for ImGui. [%s, %s, line %d]\n", __FILE__,
                __FUNCTION__, __LINE__);
            return false;
        }
#endif
    } break;
    case (GUIRenderBackend::CPU): {
        if (ImGui_ImplGeneric_Init(&this->sw_window)) { /// XXX How is sw_window used?
            imgui_sw::bind_imgui_painting();
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("[GUI] Created ImGui render Backend for CPU.");
        } else {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] Unable to initialize CPU render backend for ImGui. [%s, %s, line %d]\n", __FILE__, __FUNCTION__,
                __LINE__);
            return false;
        }
    } break;
    default: {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Unknown render backend... [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    } break;
    }

    this->initialized_backend = backend;
    return true;
}


void gui_render_backend::NewFrame(glm::vec2 framebuffer_size, glm::vec2 window_size) {

    switch (this->initialized_backend) {
    case (GUIRenderBackend::OPEN_GL): {
#ifdef WITH_GL
        ImGui_ImplOpenGL3_NewFrame();
#endif
    } break;
    case (GUIRenderBackend::CPU): {
        this->sw_window.width = static_cast<int>(window_size.x);
        this->sw_window.height = static_cast<int>(window_size.y);
        this->sw_window.x = 0;
        this->sw_window.y = 0;
        this->sw_monitor.res_x = static_cast<int>(framebuffer_size.x);
        this->sw_monitor.res_y = static_cast<int>(framebuffer_size.y);
        ImGui_ImplGeneric_NewFrame(&this->sw_window, &this->sw_monitor);
    } break;
    default: {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Unknown render backend... [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
    } break;
    }
}


bool gui_render_backend::EnableRendering(unsigned int width, unsigned int height) {

    switch (this->initialized_backend) {
    case (GUIRenderBackend::OPEN_GL): {
#ifdef WITH_GL
        auto width_i = static_cast<int>(width);
        auto height_i = static_cast<int>(height);
        bool create_fbo = false;
        if (this->ogl_fbo == nullptr) {
            create_fbo = true;
        } else if (((this->ogl_fbo->getWidth() != width_i) || (this->ogl_fbo->getHeight() != height_i)) &&
                   (width_i != 0) && (height_i != 0)) {
            create_fbo = true;
        }
        if (create_fbo) {
            try {
                glBindFramebuffer(GL_FRAMEBUFFER, 0);
                this->ogl_fbo.reset();
                this->ogl_fbo = std::make_shared<glowl::FramebufferObject>(width_i, height_i);
                this->ogl_fbo->createColorAttachment(GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE);
                // TODO: check completness and throw if not?
            } catch (glowl::FramebufferObjectException const& exc) {
                megamol::core::utility::log::Log::DefaultLog.WriteError(
                    "[GUI] Unable to create framebuffer object: %s [%s, %s, line %d]\n", exc.what(), __FILE__,
                    __FUNCTION__, __LINE__);
            }
        }
        if (this->ogl_fbo == nullptr) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[GUI] Unable to create valid FBO. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
            return false;
        }

        this->ogl_fbo->bind();
        // glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        // glClearDepth(1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0, 0, width, height);
        glEnable(GL_DEPTH_TEST);
        return true;
#endif
    } break;
    case (GUIRenderBackend::CPU): {
        bool create_fbo = false;
        if (this->cpu_fbo == nullptr) {
            create_fbo = true;
        } else if (((this->cpu_fbo->getWidth() != width) || (this->cpu_fbo->getHeight() != height)) && (width != 0) &&
                   (height != 0)) {
            create_fbo = true;
        }
        if (create_fbo) {
            this->cpu_fbo.reset();
            this->cpu_fbo = std::make_shared<megamol::core::view::CPUFramebuffer>();
            if (this->cpu_fbo == nullptr) {
                megamol::core::utility::log::Log::DefaultLog.WriteError(
                    "[GUI] Unable to create valid FBO. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
                return false;
            }
            this->cpu_fbo->colorBuffer = std::vector<uint32_t>(width * height, 0);
            this->cpu_fbo->width = width;
            this->cpu_fbo->height = height;
        }
        return true;
    } break;
    default: {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Unknown render backend... [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
    } break;
    }
    return false;
}


bool gui_render_backend::Render(ImDrawData* draw_data) {

    switch (this->initialized_backend) {
    case (GUIRenderBackend::OPEN_GL): {
#ifdef WITH_GL
        ImGui_ImplOpenGL3_RenderDrawData(draw_data);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
#endif
    } break;
    case (GUIRenderBackend::CPU): {
        std::fill_n(this->cpu_fbo->colorBuffer.data(), this->cpu_fbo->colorBuffer.size(), 0x19191919u);
        imgui_sw::paint_imgui(this->cpu_fbo->colorBuffer.data(), static_cast<int>(this->cpu_fbo->getWidth()),
            static_cast<int>(this->cpu_fbo->getHeight()));
    } break;
    default: {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Unknown render backend... [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
    } break;
    }
}


bool megamol::gui::gui_render_backend::ShutdownBackend() {

    switch (this->initialized_backend) {
    case (GUIRenderBackend::OPEN_GL): {
#ifdef WITH_GL
        ImGui_ImplOpenGL3_Shutdown();
#endif
    } break;
    case (GUIRenderBackend::CPU): {
        ImGui_ImplGeneric_Shutdown();
    } break;
    default:
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Unknown render backend... [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        break;
    }
    this->initialized_backend = GUIRenderBackend::NONE;
    return true;
}


bool gui_render_backend::CreateFont() {

    switch (this->initialized_backend) {
    case (GUIRenderBackend::NONE): {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Fonts can only be loaded after render backend was initialized. [%s, %s, line %d]\n", __FILE__,
            __FUNCTION__, __LINE__);
    } break;
    case (GUIRenderBackend::OPEN_GL): {
#ifdef WITH_GL
        return ImGui_ImplOpenGL3_CreateFontsTexture();
#endif
    } break;
    case (GUIRenderBackend::CPU): {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] CPU rendering does not support painting with any other texture than the default font texture. [%s, "
            "%s, line %d]\n",
            __FILE__, __FUNCTION__, __LINE__);
    } break;
    default: {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[GUI] Unknown render backend... [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
    } break;
    }
    return false;
}
