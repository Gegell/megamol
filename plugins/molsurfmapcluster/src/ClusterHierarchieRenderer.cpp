/*
 * ClusterRenderer.cpp
 *
 * Copyright (C) 2017 by VISUS (Universitaet Stuttgart).
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include <tuple>

#include "mmcore/param/BoolParam.h"
#include "mmcore/param/ColorParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "mmcore/view/Renderer2DModule.h"

#include "vislib/sys/Log.h"

#include "CallClusterPosition.h"
#include "ClusterHierarchieRenderer.h"
#include "DBScanClusteringProvider.h"
#include "EnzymeClassProvider.h"
#include "TextureLoader.h"

#define VIEWPORT_WIDTH 2560
#define VIEWPORT_HEIGHT 1440

using namespace megamol;
using namespace megamol::core;
using namespace megamol::core::view;
using namespace megamol::molsurfmapcluster;


/*
 * ClusterRenderer::ClusterRenderer
 */
ClusterHierarchieRenderer::ClusterHierarchieRenderer(void)
    : view::Renderer2DModule()
    , clusterDataSlot("inData", "The input data slot for sphere data.")
    , positionDataSlot("position", "The inoput data slot for the aktual position")
    , positionoutslot("getposition", "Returns the aktual Rendered-Root-Node from clustering")
    , showEnzymeClassesParam("showEnzymeClasses", "Display the Enzyme classes alongside with the popup renders")
    , showPDBIdsParam("showPDBIds", "Display the PDB Ids alongside with the popup renders")
    , fontSizeParam("fontSize", "Size of the rendered font")
    , useDistanceColors("useDistanceColors", "Use the distance colors instead of the coloring by cluster")
    , minColorParam("color::minColor", "the minimum color for interpolation")
    , midColorParam("color::midColor", "the mid color for interpolation")
    , maxColorParam("color::maxColor", "the maximum color for interpolation")
    , failColorParam("color::failColor", "color used for failed comparisons")
    , windowHeightParam("window::height", "height of the displayed window")
    , windowWidthParam("window::width", "width of the displayed window")
    , theFont(megamol::core::utility::SDFFont::FontName::ROBOTO_SANS)
    , texVa(0) {

    // Callee Slot
    this->positionoutslot.SetCallback(CallClusterPosition::ClassName(), CallClusterPosition::FunctionName(1),
        &ClusterHierarchieRenderer::GetPositionExtents);
    this->positionoutslot.SetCallback(CallClusterPosition::ClassName(), CallClusterPosition::FunctionName(0),
        &ClusterHierarchieRenderer::GetPositionData);
    this->MakeSlotAvailable(&this->positionoutslot);

    // CallerSlot
    this->clusterDataSlot.SetCompatibleCall<CallClusteringDescription>();
    this->MakeSlotAvailable(&this->clusterDataSlot);

    this->positionDataSlot.SetCompatibleCall<CallClusterPositionDescription>();
    this->MakeSlotAvailable(&this->positionDataSlot);

    // ParamSlot
    this->showEnzymeClassesParam.SetParameter(new core::param::BoolParam(false));
    this->MakeSlotAvailable(&this->showEnzymeClassesParam);

    this->showPDBIdsParam.SetParameter(new core::param::BoolParam(true));
    this->MakeSlotAvailable(&this->showPDBIdsParam);

    this->fontSizeParam.SetParameter(new core::param::FloatParam(22.0f, 5.0f, 300.0f));
    this->MakeSlotAvailable(&this->fontSizeParam);

    this->useDistanceColors.SetParameter(new param::BoolParam(false));
    this->MakeSlotAvailable(&this->useDistanceColors);

    this->minColorParam.SetParameter(new param::ColorParam(0.266666f, 0.05098f, 0.32941f, 1.0f));
    this->MakeSlotAvailable(&this->minColorParam);

    this->midColorParam.SetParameter(new param::ColorParam(0.127255f, 0.54117f, 0.55294f, 1.0f));
    this->MakeSlotAvailable(&this->midColorParam);

    this->maxColorParam.SetParameter(new param::ColorParam(0.992157f, 0.90588f, 0.14509f, 1.0f));
    this->MakeSlotAvailable(&this->maxColorParam);

    this->failColorParam.SetParameter(new param::ColorParam(1.0f, 0.0f, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->failColorParam);

    this->windowHeightParam.SetParameter(new param::IntParam(VIEWPORT_HEIGHT, 1000, 20000));
    this->MakeSlotAvailable(&this->windowHeightParam);

    this->windowWidthParam.SetParameter(new param::IntParam(VIEWPORT_WIDTH, 1000, 20000));
    this->MakeSlotAvailable(&this->windowWidthParam);

    // Variablen
    this->lastHashClustering = 0;
    this->lastHashPosition = 0;
    this->DataHashPosition = 0;
    this->clustering = nullptr;
    this->rendered = false;
    this->newposition = false;
    this->position = nullptr;

    this->newcolor = false;
    this->hashoffset = 0;
    this->colorhash = 0;

    this->actionavailable = true;
    this->dbscanclustercolor = false;

    this->popup = nullptr;
    this->x = 0;
    this->y = 0;
}


/*
 * ClusterRenderer::~ClusterRenderer
 */
ClusterHierarchieRenderer::~ClusterHierarchieRenderer(void) { this->Release(); }


/*
 * ClusterRenderer::create
 */
bool ClusterHierarchieRenderer::create(void) {

    // Initialise font
    if (!this->theFont.Initialise(this->GetCoreInstance())) {
        vislib::sys::Log::DefaultLog.WriteError("Couldn't initialize the font.");
        return false;
    }

    vislib::graphics::gl::ShaderSource texVertShader;
    vislib::graphics::gl::ShaderSource texFragShader;

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("molsurfTexture::vertex", texVertShader)) {
        vislib::sys::Log::DefaultLog.WriteMsg(
            vislib::sys::Log::LEVEL_ERROR, "Unable to load vertex shader source for texture Vertex Shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("molsurfTexture::fragment", texFragShader)) {
        vislib::sys::Log::DefaultLog.WriteMsg(
            vislib::sys::Log::LEVEL_ERROR, "Unable to load fragment shader source for texture Fragment Shader");
        return false;
    }

    try {
        if (!this->textureShader.Create(
                texVertShader.Code(), texVertShader.Count(), texFragShader.Code(), texFragShader.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        vislib::sys::Log::DefaultLog.WriteError("Unable to create shader: %s\n", e.GetMsgA());
        return false;
    }

    texVertShader.Clear();
    texFragShader.Clear();

    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("molsurfPassthrough::vertex", texVertShader)) {
        vislib::sys::Log::DefaultLog.WriteMsg(
            vislib::sys::Log::LEVEL_ERROR, "Unable to load vertex shader source for passthrough Vertex Shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource(
            "molsurfPassthrough::fragment", texFragShader)) {
        vislib::sys::Log::DefaultLog.WriteMsg(
            vislib::sys::Log::LEVEL_ERROR, "Unable to load fragment shader source for passthrough Fragment Shader");
        return false;
    }

    try {
        if (!this->passthroughShader.Create(
                texVertShader.Code(), texVertShader.Count(), texFragShader.Code(), texFragShader.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        vislib::sys::Log::DefaultLog.WriteError("Unable to create shader: %s\n", e.GetMsgA());
        return false;
    }

    const float size = 1.0f;
    std::vector<float> texVerts = {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, size, 0.0f, 0.0f, 0.0f, size, 0.0f, 0.0f, 1.0f,
        1.0f, size, size, 0.0f, 1.0f, 0.0f};

    this->texBuffer = std::make_unique<glowl::BufferObject>(GL_ARRAY_BUFFER, texVerts, GL_STATIC_DRAW);
    this->geometrySSBO = std::make_unique<glowl::BufferObject>(
        GL_SHADER_STORAGE_BUFFER, std::vector<glm::vec3>(), GL_DYNAMIC_DRAW); // this is filled later

    glGenVertexArrays(1, &this->dummyVa);
    // it is only a dummy, so we do not fill it

    glGenVertexArrays(1, &this->texVa);
    glBindVertexArray(this->texVa);

    this->texBuffer->bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), 0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);

    return true;
}


/*
 * ClusterRenderer::release
 */
void ClusterHierarchieRenderer::release(void) {
    if (this->texVa != 0) {
        glDeleteVertexArrays(1, &this->texVa);
    }
}


/*
 * ClusterRenderer::GetExtents
 */
bool ClusterHierarchieRenderer::GetExtents(view::CallRender2D& call) {

    // Incoming Call
    core::view::CallRender2D* cr = dynamic_cast<core::view::CallRender2D*>(&call);
    if (cr == nullptr) return false;

    this->windowMeasurements = cr->GetViewport();

    vislib::math::Vector<float, 2> currentViewport;
    currentViewport.SetX(static_cast<float>(this->windowWidthParam.Param<param::IntParam>()->Value()));
    currentViewport.SetY(static_cast<float>(this->windowHeightParam.Param<param::IntParam>()->Value()));

    cr->SetBoundingBox(0, 0, currentViewport.GetX(), currentViewport.GetY());

    // Check for new Data in clustering
    CallClustering* cc = this->clusterDataSlot.CallAs<CallClustering>();
    if (cc == nullptr) return false;

    if (!(*cc)(CallClustering::CallForGetExtent)) return false;

    // Check for new Position to render
    CallClusterPosition* ccp = this->positionDataSlot.CallAs<CallClusterPosition>();
    if (ccp == nullptr) return false;

    if (!(*ccp)(CallClusterPosition::CallForGetExtent)) return false;

    // if viewport changes ....
    if (currentViewport != this->viewport) {
        this->viewport = currentViewport;
    }

    return true;
}

double ClusterHierarchieRenderer::drawTree(HierarchicalClustering::CLUSTERNODE* node, glm::mat4 mvp, double minheight,
    double minwidth, double spacey, double spacex,
    std::vector<std::tuple<HierarchicalClustering::CLUSTERNODE*, ClusterRenderer::RGBCOLOR*>*>* colors) {

    double posx = 0;
    double posy = 0;
    double posLeft = 0;
    double posRight = 0;

    // draw child node
    if (node->level == 0) {
        posx = minwidth + (counter * spacex);
        posy = minheight + (node->level * spacey);
        this->counter++;

    } else {
        posLeft = drawTree(node->left, mvp, minheight, minwidth, spacey, spacex, colors);
        posRight = drawTree(node->right, mvp, minheight, minwidth, spacey, spacex, colors);

        posx = (posLeft + posRight) / 2;
        posy = minheight + (node->level * spacey);
    }

    // Select Color
    bool clusternode = false;
    glm::vec4 currentcolor;
    if (this->colors != nullptr) {
        for (std::tuple<HierarchicalClustering::CLUSTERNODE*, ClusterRenderer::RGBCOLOR*>* colortuple : *colors) {
            if (this->clustering->parentIs(node, std::get<0>(*colortuple))) {
                ClusterRenderer::RGBCOLOR* color = std::get<1>(*colortuple);
                double r = (255 - color->r) / 255;
                double g = (255 - color->g) / 255;
                double b = (255 - color->b) / 255;
                currentcolor = glm::vec4(r, g, b, 1.0f);
                clusternode = true;
            }
        }
    }

    if (this->useDistanceColors.Param<param::BoolParam>()->Value()) {
        std::string leftpdb = node->left == nullptr ? node->pic->pdbid : node->left->pic->pdbid;
        std::string rightpdb = node->right == nullptr ? node->pic->pdbid : node->right->pic->pdbid;
        auto leftvec = EnzymeClassProvider::RetrieveClassesForPdbId(leftpdb, *this->GetCoreInstance());
        auto rightvec = EnzymeClassProvider::RetrieveClassesForPdbId(rightpdb, *this->GetCoreInstance());

        float mindist = 10.0f;
        for (const auto& l : leftvec) {
            for (const auto& r : rightvec) {
                auto d = this->enzymeClassDistance(l, r);
                if (d < mindist) mindist = d;
            }
        }

        auto minColor = glm::make_vec4(this->minColorParam.Param<param::ColorParam>()->Value().data());
        auto midColor = glm::make_vec4(this->midColorParam.Param<param::ColorParam>()->Value().data());
        auto maxColor = glm::make_vec4(this->maxColorParam.Param<param::ColorParam>()->Value().data());
        auto failColor = glm::make_vec4(this->failColorParam.Param<param::ColorParam>()->Value().data());

        if (mindist > 5.0f) {
            currentcolor = failColor;
        } else {
            if (mindist <= 2.0f) {
                currentcolor = glm::mix(maxColor, midColor, mindist / 2.0f);
            } else {
                currentcolor = glm::mix(midColor, minColor, (mindist - 2.0f) / 2.0f);
            }
        }
    }

    if (this->dbscanclustercolor) {
        if (this->dbscancluster.count(node->pic->pdbid) > 0) {
            currentcolor = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
        } else {
            currentcolor = glm::vec4(0.3f, 0.3f, 0.3f, 1.0f);
        }
    }

    // Draw Point
    glPointSize(10);
    glBindVertexArray(this->dummyVa);
    this->passthroughShader.Enable();

    if (!clusternode) {
        if (this->clustering->parentIs(node, this->position)) {
            currentcolor = glm::vec4(0.5f, 0.5f, 0.5f, 1.0f);
        } else {
            currentcolor = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
        }
    }
    glUniform4f(this->passthroughShader.ParameterLocation("color"), currentcolor.x, currentcolor.y, currentcolor.z,
        currentcolor.w);
    glUniformMatrix4fv(this->passthroughShader.ParameterLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

    this->geometrySSBO->rebuffer(std::vector<glm::vec4>{glm::vec4(posx, posy, 0.0f, 1.0f)});
    this->geometrySSBO->bind(11);
    glDrawArrays(GL_POINTS, 0, 1);

    if (node->level != 0) {
        // Connect the Nodes
        double posLeftY = minheight + (node->left->level * spacey);
        double posRightY = minheight + (node->right->level * spacey);
        glLineWidth(2);
        std::vector<glm::vec4> data(6);
        data[0] = glm::vec4(posRight, posy, 0.0f, 1.0f);
        data[1] = glm::vec4(posLeft, posy, 0.0f, 1.0f);
        data[2] = glm::vec4(posRight, posy, 0.0f, 1.0f);
        data[3] = glm::vec4(posRight, posRightY, 0.0f, 1.0f);
        data[4] = glm::vec4(posLeft, posy, 0.0f, 1.0f);
        data[5] = glm::vec4(posLeft, posLeftY, 0.0f, 1.0f);
        this->geometrySSBO->rebuffer(data);
        glUniformMatrix4fv(this->passthroughShader.ParameterLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
        glDrawArrays(GL_LINES, 0, 6);
    }

    this->passthroughShader.Disable();
    glBindVertexArray(0);

    return posx;
}

/*
 * ClusterRenderer::Render
 */
bool ClusterHierarchieRenderer::Render(view::CallRender2D& call) {

    core::view::CallRender2D* cr = dynamic_cast<core::view::CallRender2D*>(&call);
    if (cr == nullptr) return false;

    this->windowMeasurements = cr->GetViewport();

    // Update data Clustering
    CallClustering* ccc = this->clusterDataSlot.CallAs<CallClustering>();
    if (!ccc) return false;
    // Updated data from cinematic camera call
    if (!(*ccc)(CallClustering::CallForGetData)) return false;

    // read matrices (old bullshit)
    GLfloat viewMatrixColumn[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, viewMatrixColumn);
    glm::mat4 view = glm::make_mat4(viewMatrixColumn);
    GLfloat projMatrixColumn[16];
    glGetFloatv(GL_PROJECTION_MATRIX, projMatrixColumn);
    glm::mat4 proj = glm::make_mat4(projMatrixColumn);
    glm::mat4 mvp = proj * view;
    this->zoomFactor = view[0][0];

    // Update data Position
    CallClusterPosition* ccp = this->positionDataSlot.CallAs<CallClusterPosition>();
    if (!ccp) return false;
    // Updated data from cinematic camera call
    if (!(*ccp)(CallClusterPosition::CallForGetData)) return false;

    if (ccc->DataHash() != this->lastHashClustering) {
        // update Clustering to work with
        this->clustering = ccc->getClustering();
        this->lastHashClustering = ccc->DataHash();
        this->root = this->clustering->getRoot();
    }

    if (ccp->DataHash() != this->lastHashPosition && !this->newposition) {
        // update Clustering to work with
        this->position = ccp->getPosition();
        this->lastHashPosition = ccp->DataHash();

        this->cluster = this->clustering->getClusterNodesOfNode(this->position);
        this->colors = ccp->getClusterColors();
    }

    if (this->clustering != nullptr) {
        if (this->clustering->finished()) {

            this->counter = 0;

            double height = this->viewport.GetY() * 0.9;
            double width = this->viewport.GetX() * 0.9;

            double minheight = this->viewport.GetY() * 0.05;
            double minwidth = this->viewport.GetX() * 0.05;

            double spacey = height / (this->root->level);
            double spacex = width / (this->clustering->getLeaves()->size() - 1);

            drawTree(this->root, mvp, minheight, minwidth, spacey, spacex, colors);

            // Render Popup
            renderPopup(mvp);
        }
    }
    return true;
}

/*
 * ClusterHierarchieRenderer::OnMouseButton
 */
bool ClusterHierarchieRenderer::OnMouseButton(megamol::core::view::MouseButton button,
    megamol::core::view::MouseButtonAction action, megamol::core::view::Modifiers mods) {

    auto down = action == MouseButtonAction::PRESS;
    auto shiftmod = mods.test(Modifier::SHIFT);
    this->mouseAction = action;
    this->mouseButton = button;

    if (actionavailable) {
        // Wenn mouse-click auf cluster => change position ...
        // Check position
        if (!shiftmod) {

            this->counter = 0;

            double height = this->viewport.GetY() * 0.9;
            double width = this->viewport.GetX() * 0.9;

            double minheight = this->viewport.GetY() * 0.05;
            double minwidth = this->viewport.GetX() * 0.05;

            double spacey = height / (this->root->level);
            double spacex = width / (this->clustering->getLeaves()->size() - 1);

            double distanceX = 30.0 / (windowMeasurements.Width() * 2.0 * this->zoomFactor);
            double distanceY = 30.0 / (windowMeasurements.Height() * 2.0 * this->zoomFactor);

            if (checkposition(this->root, this->mouseX, this->mouseY, minheight, minwidth, spacey, spacex, distanceX,
                    distanceY) == -1) {
                this->position = this->popup;
            }

            // Todo Clusterparent neu berechnen wenn nicht gesetzt...
            auto parent = this->root;
            bool change = false;
            while (this->position->clusterparent == nullptr) {
                change = false;
                auto tmpcluster = this->clustering->getClusterNodesOfNode(parent);
                for (HierarchicalClustering::CLUSTERNODE* node : *tmpcluster) {
                    if (this->position == node) {
                        this->position->clusterparent = parent;
                        change = true;
                        break;
                    } else if (this->clustering->parentIs(this->position, node)) {
                        parent = node;
                        change = true;
                        break;
                    }
                }

                if (!change) {
                    this->position->clusterparent = parent;
                }
            }
        } else {
            this->counter = 0;

            double height = this->viewport.GetY() * 0.9;
            double width = this->viewport.GetX() * 0.9;

            double minheight = this->viewport.GetY() * 0.05;
            double minwidth = this->viewport.GetX() * 0.05;

            double spacey = height / (this->root->level);
            double spacex = width / (this->clustering->getLeaves()->size() - 1);

            double distanceX = 30.0 / (windowMeasurements.Width() * 2.0 * this->zoomFactor);
            double distanceY = 30.0 / (windowMeasurements.Height() * 2.0 * this->zoomFactor);

            if (checkposition(this->root, this->mouseX, this->mouseY, minheight, minwidth, spacey, spacex, distanceX,
                    distanceY) == -1) {
                auto pdbid = this->popup->pic->pdbid;
                auto clusterids = DBScanClusteringProvider::RetrieveClustersForPdbId(pdbid, *this->GetCoreInstance());
                this->dbscancluster.clear();
                for (const auto& id : clusterids) {
                    if (id >= 0) {
                        auto ids = DBScanClusteringProvider::RetrievePDBIdsForCluster(id, *this->GetCoreInstance());
                        for (const auto val : ids) {
                            this->dbscancluster.insert(val);
                        }
                    }
                }
                this->dbscanclustercolor = true;
            } else {
                this->dbscancluster.clear();
                this->dbscanclustercolor = false;
            }
        }
    } else {
        if (action == MouseButtonAction::RELEASE) {
            actionavailable = true;
        }
        return true;
    }

    return false;
}

/*
 * ClusterHierarchieRenderer::OnMouseMove
 */
bool ClusterHierarchieRenderer::OnMouseMove(double x, double y) {
    // Only save actual Mouse Position

    this->mouseX = (float)static_cast<int>(x);
    this->mouseY = (float)static_cast<int>(y);

    // Check mouse position => if position is on cluster => render popup
    // first delete old popup boolean
    for (HierarchicalClustering::CLUSTERNODE* node : *this->clustering->getAllNodes()) {
        node->pic->popup = false;
    }
    this->popup = nullptr;

    // Check position
    this->counter = 0;

    double height = this->viewport.GetY() * 0.9;
    double width = this->viewport.GetX() * 0.9;

    double minheight = this->viewport.GetY() * 0.05;
    double minwidth = this->viewport.GetX() * 0.05;

    double spacey = height / (this->root->level);
    double spacex = width / (this->clustering->getLeaves()->size() - 1);

    double distanceX = 30.0 / (windowMeasurements.Width() * 2.0 * this->zoomFactor);
    double distanceY = 30.0 / (windowMeasurements.Height() * 2.0 * this->zoomFactor);

    checkposition(this->root, this->mouseX, this->mouseY, minheight, minwidth, spacey, spacex, distanceX, distanceY);

    this->x = this->mouseX;
    this->y = this->mouseY;

    return false;
}

void ClusterHierarchieRenderer::renderPopup(glm::mat4 mvp) {
    if (this->popup != nullptr) {
        // Load texture if not loaded
        glDisable(GL_CULL_FACE);
        glEnable(GL_TEXTURE_2D);
        this->popup->pic->popup = true;
        TextureLoader::loadTexturesToRender(this->clustering);

        // Render Texture
        // Position -> Mouse Position? Cluster Center?
        int width = static_cast<int>(1.0f / this->zoomFactor);
        if (width % 2 == 1) width++;
        int height = width / 2;

        // Berechne verschiebung
        int shiftx = 0;
        int shifty = 0;

        if (this->x + width > this->viewport.GetX()) {
            shiftx = this->viewport.GetX() - (this->x + width);
        }

        if (this->y + height > this->viewport.GetY()) {
            shifty = this->viewport.GetY() - (this->y + height);
        }

        this->popup->pic->texture->bindTexture();

        glBindVertexArray(this->texVa);
        this->textureShader.Enable();

        glUniform2f(this->textureShader.ParameterLocation("lowerleft"), this->x + shiftx, this->y + shifty);
        glUniform2f(
            this->textureShader.ParameterLocation("upperright"), this->x + shiftx + width, this->y + shifty + height);
        glUniformMatrix4fv(this->textureShader.ParameterLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
        glUniform1i(this->textureShader.ParameterLocation("tex"), 0);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

        this->textureShader.Disable();
        glBindVertexArray(0);
        glDisable(GL_TEXTURE_2D);

        auto fontSize = this->fontSizeParam.Param<param::FloatParam>()->Value();
        auto lineHeight = theFont.LineHeight(fontSize);

        // render fonts if necessary
        if (this->showPDBIdsParam.Param<param::BoolParam>()->Value()) {
            auto stringToDraw = this->popup->pic->pdbid.c_str();
            auto lineWidth = theFont.LineWidth(fontSize, stringToDraw);

            std::array<float, 4> color = {0.0f, 0.0f, 0.0f, 1.0f};
            this->theFont.DrawString(color.data(), this->x + shiftx + width * 0.5f, this->y + shifty + height * 0.5f,
                fontSize, false, stringToDraw, core::utility::AbstractFont::Alignment::ALIGN_CENTER_MIDDLE);
        }
        if (this->showEnzymeClassesParam.Param<param::BoolParam>()->Value()) {
            auto classes =
                EnzymeClassProvider::RetrieveClassesForPdbId(this->popup->pic->pdbid, *this->GetCoreInstance());
            auto numClasses = static_cast<uint32_t>(classes.size());
            float margin = 5.0f;

            if (numClasses > 0) {
                std::string text = "";
                uint32_t idx = 0;
                for (const auto& v : classes) {
                    text += (v[0] != -1 ? std::to_string(v[0]) : "") + (v[1] != -1 ? "." + std::to_string(v[1]) : "") +
                            (v[2] != -1 ? "." + std::to_string(v[2]) : "") +
                            (v[3] != -1 ? "." + std::to_string(v[3]) : "");
                    if (idx < numClasses - 1) {
                        text += "\n";
                    }
                    ++idx;
                }

                std::array<float, 4> color = {1.0f, 1.0f, 1.0f, 1.0f};
                this->theFont.DrawString(color.data(), this->x + shiftx + width + margin,
                    this->y + shifty + height * 0.5f, fontSize, false, text.c_str(),
                    core::utility::AbstractFont::Alignment::ALIGN_LEFT_MIDDLE);
            }
        }
        glEnable(GL_CULL_FACE);
    }
}

double ClusterHierarchieRenderer::checkposition(HierarchicalClustering::CLUSTERNODE* node, float x, float y,
    double minheight, double minwidth, double spacey, double spacex, double distanceX, double distanceY) {

    double posx = 0;
    double posy = 0;
    double posLeft = 0;
    double posRight = 0;

    // draw child node
    if (node->level == 0) {
        posx = minwidth + (counter * spacex);
        posy = minheight + (node->level * spacey);
        this->counter++;

        // Check position => if found return -1;
        if (x > posx - distanceX && x < posx + distanceX && y > posy - distanceY && y < posy + distanceY) {
            this->popup = node;
            return -1;
        }

    } else {
        posLeft = checkposition(node->left, x, y, minheight, minwidth, spacey, spacex, distanceX, distanceY);
        posRight = checkposition(node->right, x, y, minheight, minwidth, spacey, spacex, distanceX, distanceY);

        if (posLeft == -1 || posRight == -1) {
            return -1;
        } else {
            // Check position
            // Check position => if found return -1;

            posx = (posLeft + posRight) / 2;
            posy = minheight + (node->level * spacey);

            if (x > posx - distanceX && x < posx + distanceY && y > posy - distanceX && y < posy + distanceY) {
                this->popup = node;
                return -1;
            }
        }
    }
    return posx;
}


bool ClusterHierarchieRenderer::GetPositionExtents(Call& call) {
    CallClusterPosition* ccp = dynamic_cast<CallClusterPosition*>(&call);
    if (ccp == nullptr) return false;

    // Wenn neuer root node
    if (ccp->getPosition() != this->position) {
        this->hashoffset++;
        this->newposition = true;
        ccp->SetDataHash(this->DataHashPosition + this->hashoffset);
    }
    return true;
}

bool ClusterHierarchieRenderer::GetPositionData(Call& call) {
    CallClusterPosition* ccp = dynamic_cast<CallClusterPosition*>(&call);
    if (ccp == nullptr) return false;
    if (this->newposition) {
        ccp->setPosition(this->position);
        ccp->setClusterColors(nullptr);
        this->newposition = false;
        this->DataHashPosition += this->hashoffset;
    }
    return true;
}

float ClusterHierarchieRenderer::enzymeClassDistance(const std::array<int, 4>& arr1, const std::array<int, 4>& arr2) {
    if (arr1[0] < 0 || arr2[0] < 0) return 10.0f;
    if (arr1[0] == arr2[0]) {
        if (arr1[1] == arr2[1] || arr1[2] == arr2[2]) {
            if (arr1[1] == arr2[1] && arr1[2] == arr2[2]) {
                if (arr1[3] == arr2[3]) {
                    return 0.0f;
                }
                return 1.0f;
            }
            return 2.0f;
        }
        return 3.0f;
    }
    return 4.0f;
}
