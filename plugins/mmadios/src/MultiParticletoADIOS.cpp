/*
 * MultiParticletoADIOS.h
 *
 * Copyright (C) 2018 by Universitaet Stuttgart (VISUS).
 * Alle Rechte vorbehalten.
 */

#include "MultiParticletoADIOS.h"
#include "geometry_calls/MultiParticleDataCall.h"
#include "mmadios/CallADIOSData.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/utility/log/Log.h"
#include <algorithm>

namespace megamol {
namespace adios {

MultiParticletoADIOS::MultiParticletoADIOS(void)
        : core::Module()
        , mpSlot("mpSlot", "Slot to request multi particle data.")
        , adiosSlot("adiosSlot", "Slot to send ADIOS IO")
        , orderSlot("order", "Sets xyz parameter order") {

    this->adiosSlot.SetCallback(
        CallADIOSData::ClassName(), CallADIOSData::FunctionName(0), &MultiParticletoADIOS::getDataCallback);
    this->adiosSlot.SetCallback(
        CallADIOSData::ClassName(), CallADIOSData::FunctionName(1), &MultiParticletoADIOS::getHeaderCallback);
    this->MakeSlotAvailable(&this->adiosSlot);

    auto tmpEnum = new core::param::EnumParam(0);
    tmpEnum->SetTypePair(0, "separated");
    tmpEnum->SetTypePair(1, "interleaved");
    this->orderSlot << tmpEnum;
    this->MakeSlotAvailable(&this->orderSlot);


    this->mpSlot.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    this->MakeSlotAvailable(&this->mpSlot);
}

MultiParticletoADIOS::~MultiParticletoADIOS(void) {
    this->Release();
}

bool MultiParticletoADIOS::create(void) {
    return true;
}

void MultiParticletoADIOS::release(void) {}

bool MultiParticletoADIOS::getDataCallback(core::Call& call) {
    CallADIOSData* cad = dynamic_cast<CallADIOSData*>(&call);
    if (cad == nullptr)
        return false;

    geocalls::MultiParticleDataCall* mpdc = this->mpSlot.CallAs<geocalls::MultiParticleDataCall>();
    if (mpdc == nullptr)
        return false;

    if (!(*mpdc)(1))
        return false;

    // set frame to load from view
    mpdc->SetFrameID(cad->getFrameIDtoLoad());

    if (!(*mpdc)(0))
        return false;

    auto availVars = cad->getAvailableVars();

    // get types
    auto& list0 = mpdc->AccessParticles(0);
    size_t pCount = 0;

    if (list0.GetVertexDataType() == geocalls::MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZ) {
        if (orderSlot.Param<core::param::EnumParam>()->ValueString() == "separated") {
            auto xCont = std::make_shared<FloatContainer>(FloatContainer());
            auto yCont = std::make_shared<FloatContainer>(FloatContainer());
            auto zCont = std::make_shared<FloatContainer>(FloatContainer());
            auto lboxCont = std::make_shared<FloatContainer>(FloatContainer());
            std::vector<float>& tmp_x = xCont->getVec();
            std::vector<float>& tmp_y = yCont->getVec();
            std::vector<float>& tmp_z = zCont->getVec();
            std::vector<float>& tmp_lbox = lboxCont->getVec();
            tmp_lbox.reserve(mpdc->GetParticleListCount() * 6);
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();
                tmp_x.reserve(tmp_x.size() + num);
                tmp_y.reserve(tmp_y.size() + num);
                tmp_z.reserve(tmp_z.size() + num);
                pCount += num;
                float minx, maxx, miny, maxy, minz, maxz = 0;
                for (auto j = 0; j < num; j++) {
                    tmp_x.push_back(parts.GetParticleStore().GetXAcc()->Get_f(j));
                    tmp_y.push_back(parts.GetParticleStore().GetYAcc()->Get_f(j));
                    tmp_z.push_back(parts.GetParticleStore().GetZAcc()->Get_f(j));
                    minx = std::min(minx, tmp_x.back());
                    miny = std::min(miny, tmp_y.back());
                    minz = std::min(minz, tmp_z.back());
                    maxx = std::max(maxx, tmp_x.back());
                    maxy = std::max(maxy, tmp_y.back());
                    maxz = std::max(maxz, tmp_z.back());
                }
                tmp_lbox.insert(tmp_lbox.end(), {minx, miny, minz, maxx, maxy, maxz});
            }
            dataMap["x"] = std::move(xCont);
            dataMap["y"] = std::move(yCont);
            dataMap["z"] = std::move(zCont);
            lboxCont->shape = {mpdc->GetParticleListCount(), 6};
            dataMap["list_box"] = std::move(lboxCont);
        } else {
            auto mixCont = std::make_shared<FloatContainer>(FloatContainer());
            auto lboxCont = std::make_shared<FloatContainer>(FloatContainer());
            std::vector<float>& tmp_mix = mixCont->getVec();
            std::vector<float>& tmp_lbox = lboxCont->getVec();
            tmp_lbox.reserve(mpdc->GetParticleListCount() * 6);
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();
                tmp_mix.reserve(tmp_mix.size() + 3 * num);
                pCount += num;

                float minx, maxx, miny, maxy, minz, maxz = 0;
                for (auto j = 0; j < num; j++) {
                    auto x = parts.GetParticleStore().GetXAcc()->Get_f(j);
                    auto y = parts.GetParticleStore().GetYAcc()->Get_f(j);
                    auto z = parts.GetParticleStore().GetZAcc()->Get_f(j);
                    tmp_mix.push_back(x);
                    tmp_mix.push_back(y);
                    tmp_mix.push_back(z);
                    minx = std::min(minx, x);
                    miny = std::min(miny, y);
                    minz = std::min(minz, z);
                    maxx = std::max(maxx, x);
                    maxy = std::max(maxy, y);
                    maxz = std::max(maxz, z);
                }
                tmp_lbox.insert(tmp_lbox.end(), {minx, miny, minz, maxx, maxy, maxz});
            }
            dataMap["xyz"] = std::move(mixCont);
            lboxCont->shape = {mpdc->GetParticleListCount(), 6};
            dataMap["list_box"] = std::move(lboxCont);
        }
    } else if (list0.GetVertexDataType() == geocalls::MultiParticleDataCall::Particles::VERTDATA_DOUBLE_XYZ) {
        if (orderSlot.Param<core::param::EnumParam>()->ValueString() == "separated") {
            auto xCont = std::make_shared<DoubleContainer>(DoubleContainer());
            auto yCont = std::make_shared<DoubleContainer>(DoubleContainer());
            auto zCont = std::make_shared<DoubleContainer>(DoubleContainer());
            auto lboxCont = std::make_shared<DoubleContainer>(DoubleContainer());
            std::vector<double>& tmp_x = xCont->getVec();
            std::vector<double>& tmp_y = yCont->getVec();
            std::vector<double>& tmp_z = zCont->getVec();
            std::vector<double>& tmp_lbox = lboxCont->getVec();
            tmp_lbox.reserve(mpdc->GetParticleListCount() * 6);
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();
                tmp_x.reserve(tmp_x.size() + num);
                tmp_y.reserve(tmp_y.size() + num);
                tmp_z.reserve(tmp_z.size() + num);
                pCount += num;

                double minx, maxx, miny, maxy, minz, maxz = 0;
                for (auto j = 0; j < num; j++) {
                    tmp_x.push_back(parts.GetParticleStore().GetXAcc()->Get_d(j));
                    tmp_y.push_back(parts.GetParticleStore().GetYAcc()->Get_d(j));
                    tmp_z.push_back(parts.GetParticleStore().GetZAcc()->Get_d(j));
                    minx = std::min(minx, tmp_x.back());
                    miny = std::min(miny, tmp_y.back());
                    minz = std::min(minz, tmp_z.back());
                    maxx = std::max(maxx, tmp_x.back());
                    maxy = std::max(maxy, tmp_y.back());
                    maxz = std::max(maxz, tmp_z.back());
                }
                tmp_lbox.insert(tmp_lbox.end(), {minx, miny, minz, maxx, maxy, maxz});
            }
            dataMap["x"] = std::move(xCont);
            dataMap["y"] = std::move(yCont);
            dataMap["z"] = std::move(zCont);
            lboxCont->shape = {mpdc->GetParticleListCount(), 6};
            dataMap["list_box"] = std::move(lboxCont);
        } else {
            auto mixCont = std::make_shared<DoubleContainer>(DoubleContainer());
            auto lboxCont = std::make_shared<DoubleContainer>(DoubleContainer());
            std::vector<double>& tmp_mix = mixCont->getVec();
            std::vector<double>& tmp_lbox = lboxCont->getVec();
            tmp_lbox.reserve(mpdc->GetParticleListCount() * 6);
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();
                tmp_mix.reserve(tmp_mix.size() + 3 * num);
                pCount += num;

                double minx, maxx, miny, maxy, minz, maxz = 0;
                for (auto j = 0; j < num; j++) {
                    auto x = parts.GetParticleStore().GetXAcc()->Get_d(j);
                    auto y = parts.GetParticleStore().GetYAcc()->Get_d(j);
                    auto z = parts.GetParticleStore().GetZAcc()->Get_d(j);
                    tmp_mix.push_back(x);
                    tmp_mix.push_back(y);
                    tmp_mix.push_back(z);
                    minx = std::min(minx, x);
                    miny = std::min(miny, y);
                    minz = std::min(minz, z);
                    maxx = std::max(maxx, x);
                    maxy = std::max(maxy, y);
                    maxz = std::max(maxz, z);
                }
                tmp_lbox.insert(tmp_lbox.end(), {minx, miny, minz, maxx, maxy, maxz});
            }
            dataMap["xyz"] = std::move(mixCont);
            lboxCont->shape = {mpdc->GetParticleListCount(), 6};
            dataMap["list_box"] = std::move(lboxCont);
        }
    } else if (list0.GetVertexDataType() == geocalls::MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZR) {
        if (orderSlot.Param<core::param::EnumParam>()->ValueString() == "separated") {
            auto xCont = std::make_shared<FloatContainer>(FloatContainer());
            auto yCont = std::make_shared<FloatContainer>(FloatContainer());
            auto zCont = std::make_shared<FloatContainer>(FloatContainer());
            auto radiusCont = std::make_shared<FloatContainer>(FloatContainer());
            auto lboxCont = std::make_shared<FloatContainer>(FloatContainer());
            std::vector<float>& tmp_x = xCont->getVec();
            std::vector<float>& tmp_y = yCont->getVec();
            std::vector<float>& tmp_z = zCont->getVec();
            std::vector<float>& tmp_radius = radiusCont->getVec();
            std::vector<float>& tmp_lbox = lboxCont->getVec();
            tmp_lbox.reserve(mpdc->GetParticleListCount() * 6);
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();
                tmp_x.reserve(tmp_x.size() + num);
                tmp_y.reserve(tmp_y.size() + num);
                tmp_z.reserve(tmp_z.size() + num);
                tmp_radius.reserve(tmp_radius.size() + num);
                pCount += num;

                float minx, maxx, miny, maxy, minz, maxz = 0;
                for (auto j = 0; j < num; j++) {
                    tmp_x.push_back(parts.GetParticleStore().GetXAcc()->Get_f(j));
                    tmp_y.push_back(parts.GetParticleStore().GetYAcc()->Get_f(j));
                    tmp_z.push_back(parts.GetParticleStore().GetZAcc()->Get_f(j));
                    tmp_radius.push_back(parts.GetParticleStore().GetRAcc()->Get_f(j));
                    minx = std::min(minx, tmp_x.back());
                    miny = std::min(miny, tmp_y.back());
                    minz = std::min(minz, tmp_z.back());
                    maxx = std::max(maxx, tmp_x.back());
                    maxy = std::max(maxy, tmp_y.back());
                    maxz = std::max(maxz, tmp_z.back());
                }
                tmp_lbox.insert(tmp_lbox.end(), {minx, miny, minz, maxx, maxy, maxz});
            }
            dataMap["x"] = std::move(xCont);
            dataMap["y"] = std::move(yCont);
            dataMap["z"] = std::move(zCont);
            dataMap["radius"] = std::move(radiusCont);
            lboxCont->shape = {mpdc->GetParticleListCount(), 6};
            dataMap["list_box"] = std::move(lboxCont);
        } else {
            auto mixCont = std::make_shared<FloatContainer>(FloatContainer());
            auto radiusCont = std::make_shared<FloatContainer>(FloatContainer());
            auto lboxCont = std::make_shared<FloatContainer>(FloatContainer());
            std::vector<float>& tmp_mix = mixCont->getVec();
            std::vector<float>& tmp_radius = radiusCont->getVec();
            std::vector<float>& tmp_lbox = lboxCont->getVec();
            tmp_lbox.reserve(mpdc->GetParticleListCount() * 6);
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();
                tmp_mix.reserve(tmp_mix.size() + 3 * num);
                tmp_radius.reserve(tmp_radius.size() + num);
                pCount += num;

                float minx, maxx, miny, maxy, minz, maxz = 0;
                for (auto j = 0; j < num; j++) {
                    auto x = parts.GetParticleStore().GetXAcc()->Get_f(j);
                    auto y = parts.GetParticleStore().GetYAcc()->Get_f(j);
                    auto z = parts.GetParticleStore().GetZAcc()->Get_f(j);
                    tmp_radius.push_back(parts.GetParticleStore().GetRAcc()->Get_f(j));
                    tmp_mix.push_back(x);
                    tmp_mix.push_back(y);
                    tmp_mix.push_back(z);
                    minx = std::min(minx, x);
                    miny = std::min(miny, y);
                    minz = std::min(minz, z);
                    maxx = std::max(maxx, x);
                    maxy = std::max(maxy, y);
                    maxz = std::max(maxz, z);
                }
                tmp_lbox.insert(tmp_lbox.end(), {minx, miny, minz, maxx, maxy, maxz});
            }
            dataMap["xyz"] = std::move(mixCont);
            dataMap["radius"] = std::move(radiusCont);
            lboxCont->shape = {mpdc->GetParticleListCount(), 6};
            dataMap["list_box"] = std::move(lboxCont);
        }
    }
    if (cad->isInVars("r")) {
        if (list0.GetColourDataType() == geocalls::MultiParticleDataCall::Particles::COLDATA_UINT8_RGBA) {
            auto rCont = std::make_shared<UCharContainer>(UCharContainer());
            auto gCont = std::make_shared<UCharContainer>(UCharContainer());
            auto bCont = std::make_shared<UCharContainer>(UCharContainer());
            auto aCont = std::make_shared<UCharContainer>(UCharContainer());
            std::vector<unsigned char>& tmp_r = rCont->getVec();
            std::vector<unsigned char>& tmp_g = gCont->getVec();
            std::vector<unsigned char>& tmp_b = bCont->getVec();
            std::vector<unsigned char>& tmp_a = aCont->getVec();

            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();
                tmp_r.reserve(tmp_r.size() + num);
                tmp_g.reserve(tmp_g.size() + num);
                tmp_b.reserve(tmp_b.size() + num);
                tmp_a.reserve(tmp_a.size() + num);
                for (auto j = 0; j < num; j++) {
                    tmp_r.push_back(parts.GetParticleStore().GetCRAcc()->Get_u8(j));
                    tmp_g.push_back(parts.GetParticleStore().GetCGAcc()->Get_u8(j));
                    tmp_b.push_back(parts.GetParticleStore().GetCBAcc()->Get_u8(j));
                    tmp_a.push_back(parts.GetParticleStore().GetCAAcc()->Get_u8(j));
                }
            }
            dataMap["r"] = std::move(rCont);
            dataMap["g"] = std::move(gCont);
            dataMap["b"] = std::move(bCont);
            dataMap["a"] = std::move(aCont);

        } else if (list0.GetColourDataType() == geocalls::MultiParticleDataCall::Particles::COLDATA_FLOAT_RGB) {

            auto rCont = std::make_shared<FloatContainer>(FloatContainer());
            auto gCont = std::make_shared<FloatContainer>(FloatContainer());
            auto bCont = std::make_shared<FloatContainer>(FloatContainer());
            auto aCont = std::make_shared<FloatContainer>(FloatContainer());
            std::vector<float>& tmp_r = rCont->getVec();
            std::vector<float>& tmp_g = gCont->getVec();
            std::vector<float>& tmp_b = bCont->getVec();
            std::vector<float>& tmp_a = aCont->getVec();

            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();

                tmp_r.reserve(tmp_r.size() + num);
                tmp_g.reserve(tmp_g.size() + num);
                tmp_b.reserve(tmp_b.size() + num);
                tmp_a.reserve(tmp_a.size() + num);
                for (auto j = 0; j < num; j++) {
                    tmp_r.push_back(parts.GetParticleStore().GetCRAcc()->Get_f(j));
                    tmp_g.push_back(parts.GetParticleStore().GetCGAcc()->Get_f(j));
                    tmp_b.push_back(parts.GetParticleStore().GetCBAcc()->Get_f(j));
                    tmp_a.push_back(1.0f);
                }
            }
            dataMap["r"] = std::move(rCont);
            dataMap["g"] = std::move(gCont);
            dataMap["b"] = std::move(bCont);
            dataMap["a"] = std::move(aCont);
        } else if (list0.GetColourDataType() == geocalls::MultiParticleDataCall::Particles::COLDATA_FLOAT_RGBA) {

            auto rCont = std::make_shared<FloatContainer>(FloatContainer());
            auto gCont = std::make_shared<FloatContainer>(FloatContainer());
            auto bCont = std::make_shared<FloatContainer>(FloatContainer());
            auto aCont = std::make_shared<FloatContainer>(FloatContainer());
            std::vector<float>& tmp_r = rCont->getVec();
            std::vector<float>& tmp_g = gCont->getVec();
            std::vector<float>& tmp_b = bCont->getVec();
            std::vector<float>& tmp_a = aCont->getVec();

            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();

                tmp_r.reserve(tmp_r.size() + num);
                tmp_g.reserve(tmp_g.size() + num);
                tmp_b.reserve(tmp_b.size() + num);
                tmp_a.reserve(tmp_a.size() + num);
                for (auto j = 0; j < num; j++) {
                    tmp_r.push_back(parts.GetParticleStore().GetCRAcc()->Get_f(j));
                    tmp_g.push_back(parts.GetParticleStore().GetCGAcc()->Get_f(j));
                    tmp_b.push_back(parts.GetParticleStore().GetCBAcc()->Get_f(j));
                    tmp_a.push_back(parts.GetParticleStore().GetCAAcc()->Get_f(j));
                }
            }
            dataMap["r"] = std::move(rCont);
            dataMap["g"] = std::move(gCont);
            dataMap["b"] = std::move(bCont);
            dataMap["a"] = std::move(aCont);
        }
    } else if (cad->isInVars("global_r")) {
        auto rCont = std::make_shared<FloatContainer>(FloatContainer());
        auto gCont = std::make_shared<FloatContainer>(FloatContainer());
        auto bCont = std::make_shared<FloatContainer>(FloatContainer());
        auto aCont = std::make_shared<FloatContainer>(FloatContainer());
        std::vector<float>& tmp_r = rCont->getVec();
        std::vector<float>& tmp_g = gCont->getVec();
        std::vector<float>& tmp_b = bCont->getVec();
        std::vector<float>& tmp_a = aCont->getVec();
        tmp_r.push_back(list0.GetParticleStore().GetCRAcc()->Get_f(0));
        tmp_g.push_back(list0.GetParticleStore().GetCGAcc()->Get_f(0));
        tmp_b.push_back(list0.GetParticleStore().GetCBAcc()->Get_f(0));
        tmp_a.push_back(list0.GetParticleStore().GetCAAcc()->Get_f(0));

        rCont->singleValue = true;
        gCont->singleValue = true;
        bCont->singleValue = true;
        aCont->singleValue = true;

        dataMap["global_r"] = std::move(rCont);
        dataMap["global_g"] = std::move(gCont);
        dataMap["global_b"] = std::move(bCont);
        dataMap["global_a"] = std::move(aCont);
    } else if (cad->isInVars("list_r")) {
        auto rCont = std::make_shared<FloatContainer>(FloatContainer());
        auto gCont = std::make_shared<FloatContainer>(FloatContainer());
        auto bCont = std::make_shared<FloatContainer>(FloatContainer());
        auto aCont = std::make_shared<FloatContainer>(FloatContainer());
        std::vector<float>& tmp_r = rCont->getVec();
        std::vector<float>& tmp_g = gCont->getVec();
        std::vector<float>& tmp_b = bCont->getVec();
        std::vector<float>& tmp_a = aCont->getVec();
        tmp_r.reserve(mpdc->GetParticleListCount());
        tmp_g.reserve(mpdc->GetParticleListCount());
        tmp_b.reserve(mpdc->GetParticleListCount());
        tmp_a.reserve(mpdc->GetParticleListCount());

        for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
            geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);

            const unsigned char* rgba = parts.GetGlobalColour();
            tmp_r.push_back(static_cast<float>(rgba[0]) / 255.0f);
            tmp_g.push_back(static_cast<float>(rgba[1]) / 255.0f);
            tmp_b.push_back(static_cast<float>(rgba[2]) / 255.0f);
            tmp_a.push_back(static_cast<float>(rgba[3]) / 255.0f);
        }
        dataMap["list_r"] = std::move(rCont);
        dataMap["list_g"] = std::move(gCont);
        dataMap["list_b"] = std::move(bCont);
        dataMap["list_a"] = std::move(aCont);
    } else if (cad->isInVars("i")) {
        if (list0.GetColourDataType() == geocalls::MultiParticleDataCall::Particles::COLDATA_FLOAT_I) {

            auto iCont = std::make_shared<FloatContainer>(FloatContainer());
            std::vector<float>& tmp_i = iCont->getVec();

            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();

                tmp_i.reserve(tmp_i.size() + num);

                for (auto j = 0; j < num; j++) {
                    tmp_i.push_back(parts.GetParticleStore().GetCRAcc()->Get_f(j));
                }
            }
            dataMap["i"] = std::move(iCont);
        } else if (list0.GetColourDataType() == geocalls::MultiParticleDataCall::Particles::COLDATA_DOUBLE_I) {

            auto iCont = std::make_shared<DoubleContainer>(DoubleContainer());
            std::vector<double>& tmp_i = iCont->getVec();

            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();

                tmp_i.reserve(tmp_i.size() + num);
                for (auto j = 0; j < num; j++) {
                    tmp_i.push_back(parts.GetParticleStore().GetCRAcc()->Get_d(j));
                }
            }
            dataMap["i"] = std::move(iCont);
        }
    }

    if (list0.HasID()) {
        if (list0.GetIDDataType() == geocalls::MultiParticleDataCall::Particles::IDDATA_UINT64) {
            auto idCont = std::make_shared<UInt64Container>(UInt64Container());
            std::vector<uint64_t>& tmp_id = idCont->getVec();
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();

                tmp_id.reserve(tmp_id.size() + num);
                for (auto j = 0; j < num; j++) {
                    tmp_id.push_back(parts.GetParticleStore().GetIDAcc()->Get_u64(j));
                }
            }
            dataMap["id"] = std::move(idCont);
        } else if (list0.GetIDDataType() == geocalls::MultiParticleDataCall::Particles::IDDATA_UINT32) {
            auto idCont = std::make_shared<UInt32Container>(UInt32Container());
            std::vector<uint32_t>& tmp_id = idCont->getVec();
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                const size_t num = parts.GetCount();

                tmp_id.reserve(tmp_id.size() + num);
                for (auto j = 0; j < num; j++) {
                    tmp_id.push_back(parts.GetParticleStore().GetIDAcc()->Get_u32(j));
                }
            }
            dataMap["id"] = std::move(idCont);
        }
    }
    if (cad->isInVars("list_offset")) {
        auto listCont = std::make_shared<UInt64Container>(UInt64Container());
        std::vector<uint64_t>& tmp_list = listCont->getVec();
        tmp_list.reserve(mpdc->GetParticleListCount());
        size_t old_num = 0;

        for (uint64_t i = 0; i < mpdc->GetParticleListCount(); i++) {
            geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
            const size_t num = parts.GetCount();

            tmp_list.push_back(old_num);
            old_num = num;
        }
        dataMap["list_offset"] = std::move(listCont);
    }

    // Particle Count
    auto ic_pcount = std::make_shared<UInt64Container>(UInt64Container());
    std::vector<uint64_t>& tmp_pcount = ic_pcount->getVec();
    tmp_pcount.resize(1);
    tmp_pcount[0] = pCount;
    dataMap["count"] = std::move(ic_pcount);

    // Bounding Box
    auto fc_box = std::make_shared<FloatContainer>(FloatContainer());
    std::vector<float>& tmp_box = fc_box->getVec();
    vislib::math::Cuboid<float> bbox = mpdc->GetBoundingBoxes().ObjectSpaceBBox();
    tmp_box.resize(6);
    tmp_box = {bbox.GetLeft(), bbox.GetBottom(), bbox.GetBack(), bbox.GetRight(), bbox.GetTop(), bbox.GetFront()};
    fc_box->singleValue = true;
    dataMap["global_box"] = std::move(fc_box);


    // global radius
    if (cad->isInVars("global_radius")) {
        auto fc_radius = std::make_shared<FloatContainer>(FloatContainer());
        std::vector<float>& tmp_radius = fc_radius->getVec();
        tmp_radius.push_back(list0.GetParticleStore().GetRAcc()->Get_f(0));
        fc_radius->singleValue = true;
        dataMap["global_radius"] = std::move(fc_radius);
    } else if (cad->isInVars("list_radius")) {
        auto radiusCont = std::make_shared<FloatContainer>(FloatContainer());
        std::vector<float>& tmp_radius = radiusCont->getVec();
        tmp_radius.resize(mpdc->GetParticleListCount());
        for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
            geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);

            tmp_radius[i] = parts.GetGlobalRadius();
        }
        dataMap["list_radius"] = std::move(radiusCont);
    } else if (cad->isInVars("radius")) {
        auto radiusCont = std::make_shared<FloatContainer>(FloatContainer());
        std::vector<float>& tmp_radius = radiusCont->getVec();

        for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
            geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
            const size_t num = parts.GetCount();

            tmp_radius.reserve(tmp_radius.size() + num);
            for (auto j = 0; j < num; j++) {
                tmp_radius.push_back(parts.GetGlobalRadius());
            }
        }
        dataMap["radius"] = std::move(radiusCont);
    }


    // set stuff in call
    cad->setData(std::make_shared<adiosDataMap>(dataMap));
    cad->setDataHash(mpdc->DataHash());

    return true;
}

bool MultiParticletoADIOS::getHeaderCallback(core::Call& call) {

    CallADIOSData* cad = dynamic_cast<CallADIOSData*>(&call);
    if (cad == nullptr)
        return false;

    geocalls::MultiParticleDataCall* mpdc = this->mpSlot.CallAs<geocalls::MultiParticleDataCall>();
    if (mpdc == nullptr)
        return false;

    if (!(*mpdc)(1))
        return false;
    if (!(*mpdc)(0))
        return false;

    geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(0);

    // get total frame cound from data source
    cad->setFrameCount(mpdc->FrameCount());
    // DEBUG
    // cad->setFrameCount(1);

    // set available vars
    std::vector<std::string> availVars = {"global_box", "list_box", "count"};

    if (this->orderSlot.Param<core::param::EnumParam>()->Value() == 0) {
        availVars.push_back("x");
        availVars.push_back("y");
        availVars.push_back("z");
    } else {
        availVars.push_back("xyz");
    }

    if (parts.GetColourDataType() == geocalls::SimpleSphericalParticles::COLDATA_NONE) {
        if (mpdc->GetParticleListCount() > 1) {
            availVars.push_back("list_r");
            availVars.push_back("list_g");
            availVars.push_back("list_b");
            availVars.push_back("list_a");
        } else {
            availVars.push_back("global_r");
            availVars.push_back("global_g");
            availVars.push_back("global_b");
            availVars.push_back("global_a");
        }
    } else if ((parts.GetColourDataType() == geocalls::SimpleSphericalParticles::COLDATA_FLOAT_I) ||
               (parts.GetColourDataType() == geocalls::SimpleSphericalParticles::COLDATA_DOUBLE_I)) {
        availVars.push_back("i");
    } else {
        availVars.push_back("r");
        availVars.push_back("g");
        availVars.push_back("b");
        availVars.push_back("a");
    }

    if (parts.GetVertexDataType() == geocalls::SimpleSphericalParticles::VERTDATA_FLOAT_XYZR) {
        availVars.push_back("radius");
    } else {
        if (mpdc->GetParticleListCount() > 1) {
            std::vector<float> list_radius(mpdc->GetParticleListCount());
            for (auto i = 0; i < mpdc->GetParticleListCount(); i++) {
                geocalls::MultiParticleDataCall::Particles& parts = mpdc->AccessParticles(i);
                list_radius[i] = parts.GetGlobalRadius();
            }
            if (std::equal(list_radius.begin() + 1, list_radius.end(), list_radius.begin())) {
                availVars.push_back("global_radius");
            } else {
                availVars.push_back("list_radius");
            }
        } else {
            availVars.push_back("global_radius");
        }
    }

    availVars.push_back("list_offset");

    if (parts.HasID()) {
        availVars.push_back("id");
    }

    cad->setAvailableVars(availVars);

    return true;
}

} // end namespace adios
} // end namespace megamol
