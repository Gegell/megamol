/*
 * Call.h
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOL_GUI_GRAPH_CALL_H_INCLUDED
#define MEGAMOL_GUI_GRAPH_CALL_H_INCLUDED


#include "vislib/sys/Log.h"

#include <map>
#include <memory>
#include <vector>

#include "CallSlot.h"
#include "GUIUtils.h"


namespace megamol {
namespace gui {
namespace configurator {


// Forward declaration
class Call;
class CallSlot;
class Module;

// Pointer types to classes
typedef std::shared_ptr<Call> CallPtrType;
typedef std::shared_ptr<CallSlot> CallSlotPtrType;
typedef std::shared_ptr<Module> ModulePtrType;

/**
 * Defines call data structure for graph.
 */
class Call {
public:
    struct StockCall {
        std::string class_name;
        std::string description;
        std::string plugin_name;
        std::vector<std::string> functions;
    };

    enum Presentations : size_t { DEFAULT = 0, _COUNT_ = 1 };

    Call(ImGuiID uid);
    ~Call();

    const ImGuiID uid;

    // Init when adding call from stock
    std::string class_name;
    std::string description;
    std::string plugin_name;
    std::vector<std::string> functions;

    bool IsConnected(void);
    bool ConnectCallSlots(CallSlotPtrType call_slot_1, CallSlotPtrType call_slot_2);
    bool DisConnectCallSlots(void);
    const CallSlotPtrType GetCallSlot(CallSlot::CallSlotType type);

    // GUI Presentation -------------------------------------------------------

    // Returns uid if the call is selected.
    ImGuiID GUI_Present(const CanvasType& in_canvas, HotKeyArrayType& inout_hotkeys) {
        return this->present.Present(*this, in_canvas, inout_hotkeys);
    }

    void GUI_SetLabelVisibility(bool visible) { this->present.label_visible = visible; }
    void GUI_SetPresentation(Call::Presentations present) { this->present.presentations = present; }

private:
    std::map<CallSlot::CallSlotType, CallSlotPtrType> connected_call_slots;

    /**
     * Defines GUI call presentation.
     */
    class Presentation {
    public:
        Presentation(void);

        ~Presentation(void);

        ImGuiID Present(Call& inout_call, const CanvasType& in_canvas, HotKeyArrayType& inout_hotkeys);

        Call::Presentations presentations;
        bool label_visible;

    private:
        GUIUtils utils;
        bool selected;

    } present;
};


} // namespace configurator
} // namespace gui
} // namespace megamol

#endif // MEGAMOL_GUI_GRAPH_CALL_H_INCLUDED