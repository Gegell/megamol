/*
 * InterfaceSlot.h
 *
 * Copyright (C) 2020 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOL_GUI_GRAPH_INTERFACESLOT_H_INCLUDED
#define MEGAMOL_GUI_GRAPH_INTERFACESLOT_H_INCLUDED


#include "GUIUtils.h"
#include "widgets/HoverToolTip.h"


namespace megamol {
namespace gui {


// Forward declarations
class InterfaceSlot;
class CallSlot;
#ifndef _CALL_SLOT_TYPE_
enum CallSlotType { CALLEE, CALLER };
#    define _CALL_SLOT_TYPE_
#endif
typedef std::vector<CallSlotPtrType> CallSlotPtrVectorType;

// Types
typedef std::shared_ptr<InterfaceSlot> InterfaceSlotPtrType;
typedef std::vector<InterfaceSlotPtrType> InterfaceSlotPtrVectorType;
typedef std::map<CallSlotType, InterfaceSlotPtrVectorType> InterfaceSlotPtrMapType;


/** ************************************************************************
 * Defines GUI call slot presentation.
 */
class InterfaceSlotPresentation {
public:
    friend class InterfaceSlot;

    struct GroupState {
        ImGuiID uid;
        bool collapsed_view;
    };

    // VARIABLES --------------------------------------------------------------

    GroupState group;
    bool label_visible;


    // Widgets
    HoverToolTip tooltip;

    // FUNCTIONS --------------------------------------------------------------

    InterfaceSlotPresentation(void);
    ~InterfaceSlotPresentation(void);

    std::string GetLabel(void) { return this->label; }
    ImVec2 GetPosition(InterfaceSlot& inout_interfaceslot);
    inline bool IsGroupViewCollapsed(void) { return this->group.collapsed_view; }

    void SetPosition(ImVec2 pos) { this->position = pos; }

private:
    // VARIABLES --------------------------------------------------------------

    bool selected;
    std::string label;
    ImGuiID last_compat_callslot_uid;
    ImGuiID last_compat_interface_uid;
    bool compatible;
    // Absolute position including canvas offset and zooming
    ImVec2 position;

    // FUNCTIONS --------------------------------------------------------------

    void Present(megamol::gui::PresentPhase phase, InterfaceSlot& inout_interfaceslot, GraphItemsStateType& state);
};


/** ************************************************************************
 * Defines group interface slots bundling and redirecting calls of compatible call slots.
 */
class InterfaceSlot {
public:
    // VARIABLES --------------------------------------------------------------

    const ImGuiID uid;
    InterfaceSlotPresentation present;


    // FUNCTIONS --------------------------------------------------------------

    InterfaceSlot(ImGuiID uid, bool auto_create);
    ~InterfaceSlot();

    bool AddCallSlot(const CallSlotPtrType& callslot_ptr, const InterfaceSlotPtrType& parent_interfaceslot_ptr);
    bool RemoveCallSlot(ImGuiID callslot_uid);
    bool ContainsCallSlot(ImGuiID callslot_uid);
    bool IsConnectionValid(CallSlot& callslot);
    bool IsConnectionValid(InterfaceSlot& interfaceslot);
    bool GetCompatibleCallSlot(CallSlotPtrType& out_callslot_ptr);
    CallSlotPtrVectorType& GetCallSlots(void) { return this->callslots; }
    bool IsConnected(void);
    CallSlotType GetCallSlotType(void);
    bool IsEmpty(void);
    bool IsAutoCreated(void) { return this->auto_created; }

    // Presentation ----------------------------------------------------

    inline void PresentGUI(megamol::gui::PresentPhase phase, GraphItemsStateType& state) {
        this->present.Present(phase, *this, state);
    }
    inline ImVec2 GetGUIPosition(void) { return this->present.GetPosition(*this); }

private:
    // VARIABLES --------------------------------------------------------------

    bool auto_created;
    CallSlotPtrVectorType callslots;

    // FUNCTIONS --------------------------------------------------------------

    bool is_callslot_compatible(CallSlot& callslot);
};


} // namespace gui
} // namespace megamol

#endif // MEGAMOL_GUI_GRAPH_INTERFACESLOT_H_INCLUDED