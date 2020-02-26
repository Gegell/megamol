/*
 * ClusterRenderer.h
 *
 * Copyright (C) 2017 by VISUS (Universitaet Stuttgart).
 * Alle Rechte vorbehalten.
 */

#ifndef MOLSURFMAPCLUSTER_CLUSTERHIERARCHIERENDERER_INCLUDED
#define MOLSURFMAPCLUSTER_CLUSTERHIERARCHIERENDERER_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#    pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "mmcore/CallerSlot.h"
#include "mmcore/utility/SDFFont.h"
#include "mmcore/view/CallRender2D.h"
#include "mmcore/view/Renderer2DModule.h"

#include "vislib/graphics/gl/GLSLShader.h"
#include "vislib/graphics/gl/ShaderSource.h"

#include "vislib/Array.h"

#include "CallCluster.h"
#include "ClusterRenderer.h"

namespace megamol {
namespace MolSurfMapCluster {

/**
 * Mesh-based renderer for b�zier curve tubes
 */
class ClusterHierarchieRenderer : public core::view::Renderer2DModule {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "ClusterHierarchieRenderer"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) { return "Render the Hierarchie of the given Clustering"; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) { return true; }

    /** Ctor. */
    ClusterHierarchieRenderer(void);

    /** Dtor. */
    virtual ~ClusterHierarchieRenderer(void);

    /** The mouse button pressed/released callback. */
    virtual bool OnMouseButton(megamol::core::view::MouseButton button, megamol::core::view::MouseButtonAction action,
        megamol::core::view::Modifiers mods) override;

    /** The mouse movement callback. */
    virtual bool OnMouseMove(double x, double y) override;

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    virtual bool create(void);

    /**
     * Implementation of 'Release'.
     */
    virtual void release(void);

    /**
     * The get extents callback. The module should set the members of
     * 'call' to tell the caller the extents of its data (bounding boxes
     * and times).
     */
    virtual bool GetExtents(core::view::CallRender2D& call);

    /**
     * The render callback.
     */
    virtual bool Render(core::view::CallRender2D& call);

    /**
     * The GetData and GetExtents Callback for Position
     */
    virtual bool GetPositionExtents(core::Call& call);
    virtual bool GetPositionData(core::Call& call);

private:
    /**********************************************************************
     * variables
     **********************************************************************/
    /**The Viewport*/
    vislib::math::Vector<float, 2> viewport;

    // Font rendering
    megamol::core::utility::SDFFont theFont;
    float fontSize;

    /** DataHash*/
    SIZE_T lastHashClustering;
    SIZE_T lastHashPosition;
    SIZE_T DataHashPosition;

    /** Clustering*/
    HierarchicalClustering* clustering;
    HierarchicalClustering::CLUSTERNODE* position;
    HierarchicalClustering::CLUSTERNODE* root;

    std::vector<HierarchicalClustering::CLUSTERNODE*>* cluster;
    std::vector<std::tuple<HierarchicalClustering::CLUSTERNODE*, ClusterRenderer::RGBCOLOR*>*>* colors;

    bool rendered;
    unsigned int counter;

    bool newcolor;
    SIZE_T hashoffset;
    SIZE_T colorhash;

    vislib::graphics::gl::GLSLShader textureShader;
    std::unique_ptr<glowl::BufferObject> texBuffer;

    GLuint texVa;

    // Popup
    HierarchicalClustering::CLUSTERNODE* popup;
    int x;
    int y;


    /* Mouse Variablen*/
    /** The current mouse coordinates */
    float mouseX;
    float mouseY;

    /** The last mouse coordinates */
    float lastMouseX;
    float lastMouseY;

    core::view::MouseButton mouseButton;
    core::view::MouseButtonAction mouseAction;

    bool actionavailable;


    /*** INPUT ********************************************************/


    /**********************************************************************
     * functions
     **********************************************************************/
    double drawTree(HierarchicalClustering::CLUSTERNODE*, glm::mat4, double, double, double, double,
        std::vector<std::tuple<HierarchicalClustering::CLUSTERNODE*, ClusterRenderer::RGBCOLOR*>*>*);
    void renderPopup(glm::mat4);
    double checkposition(HierarchicalClustering::CLUSTERNODE*, float, float, double, double, double, double);


    /**********************************************************************
     * callback stuff
     **********************************************************************/

    /** The input data slot. */
    core::CallerSlot clusterDataSlot;
    core::CallerSlot positionDataSlot;

    core::CalleeSlot positionoutslot;
    bool newposition;
};

} // namespace MolSurfMapCluster
} // namespace megamol

#endif /*MOLSURFMAPCLUSTER_CLUSTERHIERARCHIERENDERER_INCLUDED*/
