/*
 * MinMax.cpp
 *
 * Copyright (C) 2019 by VISUS (Universitaet Stuttgart)
 * All rights reserved.
 */

#include "stdafx.h"
#include "MinMax.h"

#include "vislib/sys/Log.h"

/*
 * megamol::astro::MinMax::MinMax
 */
megamol::astro::MinMax::MinMax(void) : Module(),
        slotVolumetricDataIn("volumetricDataIn", "Input slot for volumetric data"),
        slotVolumetricDataOut("volumetricDataOut", "Output slot for volumetric data") {
    // Publish the slots.
    this->slotVolumetricDataIn.SetCompatibleCall<megamol::core::misc::VolumetricDataCallDescription>();
    this->MakeSlotAvailable(&this->slotVolumetricDataIn);

    this->slotVolumetricDataOut.SetCallback(megamol::core::misc::VolumetricDataCall::ClassName(),
            megamol::core::misc::VolumetricDataCall::FunctionName(megamol::core::misc::VolumetricDataCall::IDX_GET_DATA), &MinMax::onGetData);
    this->slotVolumetricDataOut.SetCallback(megamol::core::misc::VolumetricDataCall::ClassName(),
            megamol::core::misc::VolumetricDataCall::FunctionName(megamol::core::misc::VolumetricDataCall::IDX_GET_EXTENTS), &MinMax::onGetExtents);
    this->slotVolumetricDataOut.SetCallback(megamol::core::misc::VolumetricDataCall::ClassName(),
            megamol::core::misc::VolumetricDataCall::FunctionName(megamol::core::misc::VolumetricDataCall::IDX_GET_METADATA), &MinMax::onGetMetadata);
    this->slotVolumetricDataOut.SetCallback(megamol::core::misc::VolumetricDataCall::ClassName(),
            megamol::core::misc::VolumetricDataCall::FunctionName(megamol::core::misc::VolumetricDataCall::IDX_START_ASYNC), &MinMax::onUnsupportedCallback);
    this->slotVolumetricDataOut.SetCallback(megamol::core::misc::VolumetricDataCall::ClassName(),
            megamol::core::misc::VolumetricDataCall::FunctionName(megamol::core::misc::VolumetricDataCall::IDX_STOP_ASYNC), &MinMax::onUnsupportedCallback);
    this->slotVolumetricDataOut.SetCallback(megamol::core::misc::VolumetricDataCall::ClassName(),
            megamol::core::misc::VolumetricDataCall::FunctionName(megamol::core::misc::VolumetricDataCall::IDX_TRY_GET_DATA), &MinMax::onUnsupportedCallback);
    this->MakeSlotAvailable(&this->slotVolumetricDataOut);
}

/*
 * megamol::astro::MinMax::~MinMax
 */
megamol::astro::MinMax::~MinMax(void) {
    this->Release();
}

/*
 * megamol::astro::MinMax::create
 */
bool megamol::astro::MinMax::create(void) {
    return true;
}

/*
 * megamol::astro::MinMax::release
 */
void megamol::astro::MinMax::release(void) { }

bool megamol::astro::MinMax::onGetData(megamol::core::Call &call) {
    return pipeVolumetricDataCall(call, megamol::core::misc::VolumetricDataCall::IDX_GET_DATA);
}

bool megamol::astro::MinMax::onGetExtents(megamol::core::Call &call) {
    return pipeVolumetricDataCall(call, megamol::core::misc::VolumetricDataCall::IDX_GET_EXTENTS);
}

bool megamol::astro::MinMax::onGetMetadata(megamol::core::Call &call) {
    return pipeVolumetricDataCall(call, megamol::core::misc::VolumetricDataCall::IDX_GET_METADATA);
}

bool megamol::astro::MinMax::onUnsupportedCallback(megamol::core::Call &call) {
    return false;
}

bool megamol::astro::MinMax::pipeVolumetricDataCall(megamol::core::Call &call, unsigned int funcIdx) {
    using megamol::core::misc::VolumetricDataCall;
    using vislib::sys::Log;

    auto dst = dynamic_cast<VolumetricDataCall*>(&call);
    auto src = this->slotVolumetricDataIn.CallAs<VolumetricDataCall>();

    if (dst == nullptr) {
        Log::DefaultLog.WriteError(L"Call %hs of %hs received a wrong request.",
                                   VolumetricDataCall::FunctionName(funcIdx),
                                   MinMax::ClassName());
        return false;
    }
    if (src == nullptr) {
        Log::DefaultLog.WriteError(L"Call %hs of %hs has a wrong source.",
                                   VolumetricDataCall::FunctionName(funcIdx),
                                   MinMax::ClassName());
        return false;
    }

    *src = *dst;
    if (!(*src)(funcIdx)) {
        Log::DefaultLog.WriteError(L"%hs failed to call %hs.",
                                   MinMax::ClassName(),
                                   VolumetricDataCall::FunctionName(funcIdx));

        return false;
    }
    *dst = *src;

    return true;
}