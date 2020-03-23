/*
 * Module.cpp
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"

#include "Call.h"
#include "CallSlot.h"
#include "Module.h"


using namespace megamol;
using namespace megamol::gui;
using namespace megamol::gui::configurator;


megamol::gui::configurator::Module::Module(ImGuiID uid) : uid(uid), present() {

    this->call_slots.clear();
    this->call_slots.emplace(
        megamol::gui::configurator::CallSlot::CallSlotType::CALLER, std::vector<CallSlotPtrType>());
    this->call_slots.emplace(
        megamol::gui::configurator::CallSlot::CallSlotType::CALLEE, std::vector<CallSlotPtrType>());
}


megamol::gui::configurator::Module::~Module() { this->RemoveAllCallSlots(); }


bool megamol::gui::configurator::Module::AddCallSlot(megamol::gui::configurator::CallSlotPtrType call_slot) {

    if (call_slot == nullptr) {
        vislib::sys::Log::DefaultLog.WriteWarn(
            "Pointer to given call slot is nullptr. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }
    auto type = call_slot->type;
    for (auto& call_slot_ptr : this->call_slots[type]) {
        if (call_slot_ptr == call_slot) {
            throw std::invalid_argument("Pointer to call slot already registered in modules call slot list.");
        }
    }
    this->call_slots[type].emplace_back(call_slot);
    return true;
}


bool megamol::gui::configurator::Module::RemoveAllCallSlots(void) {

    try {
        for (auto& call_slots_map : this->call_slots) {
            for (auto& call_slot_ptr : call_slots_map.second) {
                if (call_slot_ptr == nullptr) {
                    // vislib::sys::Log::DefaultLog.WriteWarn(
                    //     "Call slot is already disconnected. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
                } else {
                    call_slot_ptr->DisConnectCalls();
                    call_slot_ptr->DisConnectParentModule();

                    // vislib::sys::Log::DefaultLog.WriteWarn(
                    //      "Found %i references pointing to call slot. [%s, %s, line %d]\n", call_slot_ptr.use_count(),
                    //      __FILE__, __FUNCTION__, __LINE__);
                    assert(call_slot_ptr.use_count() == 1);

                    call_slot_ptr.reset();
                }
            }
            call_slots_map.second.clear();
        }
    } catch (std::exception e) {
        vislib::sys::Log::DefaultLog.WriteError(
            "Error: %s [%s, %s, line %d]\n", e.what(), __FILE__, __FUNCTION__, __LINE__);
        return false;
    } catch (...) {
        vislib::sys::Log::DefaultLog.WriteError("Unknown Error. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }
    return true;
}


const CallSlotPtrType megamol::gui::configurator::Module::GetCallSlot(ImGuiID call_slot_uid) {

    if (call_slot_uid != GUI_INVALID_ID) {
        for (auto& call_slot_map : this->GetCallSlots()) {
            for (auto& call_slot : call_slot_map.second) {
                if (call_slot->uid == call_slot_uid) {
                    return call_slot;
                }
            }
        }
    }

    return nullptr;
}


const std::vector<megamol::gui::configurator::CallSlotPtrType>& megamol::gui::configurator::Module::GetCallSlots(
    megamol::gui::configurator::CallSlot::CallSlotType type) {

    // if (this->call_slots[type].empty()) {
    //    vislib::sys::Log::DefaultLog.WriteWarn(
    //        "Returned call slot list is empty. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
    //}
    return this->call_slots[type];
}


const std::map<megamol::gui::configurator::CallSlot::CallSlotType,
    std::vector<megamol::gui::configurator::CallSlotPtrType>>&
megamol::gui::configurator::Module::GetCallSlots(void) {

    return this->call_slots;
}


// MODULE PRESENTATION ####################################################

megamol::gui::configurator::Module::Presentation::Presentation(void)
    : presentations(Module::Presentations::DEFAULT)
    , label_visible(true)
    , position(ImVec2(FLT_MAX, FLT_MAX))
    , size(ImVec2(0.0f, 0.0f))
    , class_label()
    , name_label()
    , utils()
    , selected(false)
    , module_updated(true) {}


megamol::gui::configurator::Module::Presentation::~Presentation(void) {}


ImGuiID megamol::gui::configurator::Module::Presentation::Present(megamol::gui::configurator::Module& inout_mod,
    const CanvasType& in_canvas, megamol::gui::HotKeyArrayType& inout_hotkeys,
    megamol::gui::configurator::CallSlot::InteractType& inout_slot_interact) {

    if (ImGui::GetCurrentContext() == nullptr) {
        vislib::sys::Log::DefaultLog.WriteError(
            "No ImGui context available. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }
    ImGuiID retval_id = GUI_INVALID_ID;
    bool popup_rename = false;
    ImGuiStyle& style = ImGui::GetStyle();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    assert(draw_list != nullptr);

    try {
        // Condition for initialization position (if position is not set yet via tag in project file)
        if ((this->position.x == FLT_MAX) && (this->position.y == FLT_MAX)) {
            this->position = ImVec2(10.0f, 10.0f) + (ImGui::GetWindowPos() - in_canvas.offset) / in_canvas.zooming;
        }

        // Add module update state to general canvas update state
        auto canvas_update_state = in_canvas;
        canvas_update_state.updated = (canvas_update_state.updated || this->module_updated);
        this->module_updated = false;

        // Trigger only when canvas was updated
        // Always update position before clipping -> calls need updated slot positions.
        if (canvas_update_state.updated || (this->size.x == 0.0f) || (this->size.y == 0.0f)) {
            this->UpdateSize(inout_mod, in_canvas.zooming);
        }

        // Draw call slots ----------------------------------------------------
        /// Draw call slots prior to modules to catch mouse clicks on slot area lying over module box.
        ImGuiID module_slot_hovered_uid = inout_slot_interact.out_hovered_uid;
        for (auto& slot_pair : inout_mod.GetCallSlots()) {
            for (auto& slot : slot_pair.second) {
                slot->GUI_Present(canvas_update_state, inout_slot_interact);
            }
        }
        // Register hovering of call slots only being part of current module
        bool module_slot_hovered = (module_slot_hovered_uid != inout_slot_interact.out_hovered_uid);

        // Draw module --------------------------------------------------------
        ImVec2 module_size = this->size * in_canvas.zooming;
        ImVec2 module_rect_min = in_canvas.offset + this->position * in_canvas.zooming;
        ImVec2 module_rect_max = module_rect_min + module_size;
        ImVec2 module_center = module_rect_min + ImVec2(module_size.x / 2.0f, module_size.y / 2.0f);

        // Clip module if lying ouside the canvas
        ImVec2 canvas_rect_min = in_canvas.position;
        ImVec2 canvas_rect_max = in_canvas.position + in_canvas.size;
        if (!((canvas_rect_min.x < module_rect_max.x) && (canvas_rect_max.x > module_rect_min.x) &&
                (canvas_rect_min.y < module_rect_max.y) && (canvas_rect_max.y > module_rect_min.y))) {
            retval_id = GUI_INVALID_ID;
            if (ImGui::GetIO().MouseClicked[0]) {
                this->selected = false;
            }
            if (this->selected) {
                retval_id = inout_mod.uid;
            }
            return retval_id;
        }

        ImGui::PushID(inout_mod.uid);

        ImVec4 tmpcol = style.Colors[ImGuiCol_Button];
        // tmpcol = ImVec4(tmpcol.x * tmpcol.w, tmpcol.y * tmpcol.w, tmpcol.z * tmpcol.w, 1.0f);
        const ImU32 COLOR_MODULE_BACKGROUND = ImGui::ColorConvertFloat4ToU32(tmpcol);
        const ImU32 COLOR_MODULE_HIGHTLIGHT = ImGui::ColorConvertFloat4ToU32(style.Colors[ImGuiCol_ButtonActive]);
        const ImU32 COLOR_MODULE_BORDER = ImGui::ColorConvertFloat4ToU32(style.Colors[ImGuiCol_PopupBg]);

        // Draw text
        std::string label;
        if (this->label_visible) {
            draw_list->ChannelsSetCurrent(1); // Foreground

            ImGui::BeginGroup();

            float line_offset = 0.0f;
            if (inout_mod.is_view_instance) {
                line_offset = -0.5f * ImGui::GetTextLineHeightWithSpacing();
            }

            label = this->class_label;
            float name_width = this->utils.TextWidgetWidth(label);
            ImGui::SetCursorScreenPos(
                module_center + ImVec2(-(name_width / 2.0f), line_offset - ImGui::GetTextLineHeightWithSpacing()));
            ImGui::Text(label.c_str());

            label = this->name_label;
            name_width = this->utils.TextWidgetWidth(label);
            ImGui::SetCursorScreenPos(module_center + ImVec2(-(name_width / 2.0f), line_offset));
            ImGui::Text(label.c_str());

            if (inout_mod.is_view_instance) {
                label = "[Main View]";
                name_width = this->utils.TextWidgetWidth(label);
                ImGui::SetCursorScreenPos(
                    module_center + ImVec2(-(name_width / 2.0f), line_offset + ImGui::GetTextLineHeightWithSpacing()));
                ImGui::Text(label.c_str());
            }
            ImGui::EndGroup();
        }

        // Draw box
        draw_list->ChannelsSetCurrent(0); // Background

        ImGui::SetCursorScreenPos(module_rect_min);
        label = "module_" + inout_mod.name;
        ImGui::InvisibleButton(label.c_str(), module_size);
        ImGui::SetItemAllowOverlap();
        bool hovered = ImGui::IsItemHovered() && (!module_slot_hovered);
        bool mouse_clicked = ImGui::GetIO().MouseClicked[0];
        if (mouse_clicked && (!hovered || (module_slot_hovered))) {
            this->selected = false;
        }
        // Gives slots which overlap modules priority for ToolTip and Context Menu.
        if (!module_slot_hovered) {
            std::string hover_text = inout_mod.description;
            if (!this->label_visible) {
                hover_text = "[" + inout_mod.name + "] " + hover_text;
            }
            this->utils.HoverToolTip(hover_text.c_str(), ImGui::GetID(label.c_str()), 0.5f, 5.0f);
            // Context menu
            if (ImGui::BeginPopupContextItem("invisible_button_context")) {
                if (ImGui::MenuItem(
                        "Delete", std::get<0>(inout_hotkeys[HotkeyIndex::DELETE_GRAPH_ITEM]).ToString().c_str())) {
                    std::get<1>(inout_hotkeys[HotkeyIndex::DELETE_GRAPH_ITEM]) = true;
                    retval_id = inout_mod.uid;
                }
                if (ImGui::MenuItem("Rename")) {
                    popup_rename = true;
                }
                ImGui::EndPopup();
            }
            bool active = ImGui::IsItemActive();
            if (active) {
                this->selected = true;
                if (ImGui::IsMouseDragging(0)) {
                    this->position =
                        ((module_rect_min - in_canvas.offset) + ImGui::GetIO().MouseDelta) / in_canvas.zooming;
                    this->module_updated = true;
                }
            }
        }
        ImU32 module_bg_color = (hovered || this->selected) ? COLOR_MODULE_HIGHTLIGHT : COLOR_MODULE_BACKGROUND;
        draw_list->AddRectFilled(module_rect_min, module_rect_max, module_bg_color, 5.0f);
        draw_list->AddRect(module_rect_min, module_rect_max, COLOR_MODULE_BORDER, 5.0f);

        /// XXX Use ImGui::ArrowButton to show/hide parameters inside module box.

        // Rename pop-up
        if (this->utils.RenamePopUp("Rename Project", popup_rename, inout_mod.name)) {
            this->module_updated = true;
        }

        if (this->selected) {
            retval_id = inout_mod.uid;
        }

        ImGui::PopID();

    } catch (std::exception e) {
        vislib::sys::Log::DefaultLog.WriteError(
            "Error: %s [%s, %s, line %d]\n", e.what(), __FILE__, __FUNCTION__, __LINE__);
        return GUI_INVALID_ID;
    } catch (...) {
        vislib::sys::Log::DefaultLog.WriteError("Unknown Error. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return GUI_INVALID_ID;
    }

    return retval_id;
}


void megamol::gui::configurator::Module::Presentation::UpdateSize(
    megamol::gui::configurator::Module& mod, float canvas_zooming) {

    float max_label_length = 0.0f;
    if (this->label_visible) {
        this->class_label = " Class: " + mod.class_name + " ";
        float class_name_length = this->utils.TextWidgetWidth(this->class_label);
        this->name_label = " Name: " + mod.FullName() + " ";
        float name_length = this->utils.TextWidgetWidth(mod.present.name_label);
        max_label_length = std::max(class_name_length, name_length);
    }
    max_label_length /= canvas_zooming;

    float max_slot_name_length = 0.0f;
    for (auto& call_slot_type_list : mod.GetCallSlots()) {
        for (auto& call_slot : call_slot_type_list.second) {
            if (call_slot->GUI_GetLabelVisibility()) {
                max_slot_name_length = std::max(this->utils.TextWidgetWidth(call_slot->name), max_slot_name_length);
            }
        }
    }
    if (max_slot_name_length != 0.0f) {
        max_slot_name_length = (2.0f * max_slot_name_length / canvas_zooming) + (4.0f * GUI_CALL_SLOT_RADIUS);
    }

    float module_width = (max_label_length + max_slot_name_length) + (2.0f * GUI_CALL_SLOT_RADIUS);

    auto max_slot_count = std::max(mod.GetCallSlots(CallSlot::CallSlotType::CALLEE).size(),
        mod.GetCallSlots(CallSlot::CallSlotType::CALLER).size());
    float module_slot_height =
        (static_cast<float>(max_slot_count) * (GUI_CALL_SLOT_RADIUS * 2.0f) * 1.5f) + GUI_CALL_SLOT_RADIUS;

    float module_height = std::max(module_slot_height,
        (1.0f / canvas_zooming) * (ImGui::GetTextLineHeightWithSpacing() * ((mod.is_view_instance) ? (4.0f) : (3.0f))));

    // Clamp to minimum size
    this->size = ImVec2(std::max(module_width, 75.0f), std::max(module_height, 25.0f));
}
