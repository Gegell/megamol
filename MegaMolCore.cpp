/*
 * MegaMolCore.cpp
 *
 * Copyright (C) 2006 - 2008 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "api/MegaMolCore.h"

#define _LOG_CORE_HASH_INFO 1
#define _SEND_CORE_HASH_INFO 1

#include "mmd3d.h"
#include "ApiHandle.h"
#include "CoreInstance.h"
#include "JobDescription.h"
#include "JobInstance.h"
#include "api/MegaMolCore.inl"
#include "ObjectDescription.h"
#include "ObjectDescriptionManager.h"
#include "versioninfo.h"
#include "param/ParamHandle.h"
#include "utility/Configuration.h"
#include "ViewDescription.h"
#include "ViewInstance.h"
#include "view/AbstractTileView.h"
#include "view/AbstractView.h"
#include "view/ViewDirect3D.h"
#include "job/AbstractJob.h"
#include "ModuleDescriptionManager.h"
#include "CallDescriptionManager.h"
#include "CallerSlot.h"

#include "vislib/assert.h"
#include "vislib/CriticalSection.h"
#include "vislib/Console.h"
#include "vislib/File.h"
#include "vislib/functioncast.h"
#include "vislib/Log.h"
#include "vislib/mathfunctions.h"
#include "vislib/MD5HashProvider.h"
#include "vislib/SHA1HashProvider.h"
#include "vislib/Path.h"
#include "vislib/String.h"
#include "vislib/StringConverter.h"
#include "vislib/SystemInformation.h"
#include "vislib/ThreadSafeStackTrace.h"
#include "vislib/Trace.h"
#include "vislib/Socket.h"


#ifdef _WIN32
/* windows dll entry point */
#ifdef _MANAGED
#pragma managed(push, off)
#endif /* _MANAGED */

HMODULE mmCoreModuleHandle;

BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call,
        LPVOID lpReserved) {
    mmCoreModuleHandle = hModule;
    switch (ul_reason_for_call) {
        case DLL_PROCESS_ATTACH:
        case DLL_THREAD_ATTACH:
        case DLL_THREAD_DETACH:
        case DLL_PROCESS_DETACH:
            break;
    }
    return TRUE;
}

#ifdef _MANAGED
#pragma managed(pop)
#endif /* _MANAGED */

#else /* _WIN32 */
/* linux shared object main */

extern "C" {

const char interp[] __attribute__((section(".interp"))) = 
"/lib/ld-linux.so.2";

void mmCoreMain(int argc, char *argv[]) {
    printf("Horst!\n");
    //printf("argc = %i (%u)\nargv = %p\n", argc, argc, argv);
    //printf("*argv = %s\n", *argv);
    exit(0);
}

}

#endif /* _WIN32 */


/*
 * mmvGetVersionInfo
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcGetVersionInfo(
        unsigned short *outMajorVersion, unsigned short *outMinorVersion,
        unsigned short *outMajorRevision, unsigned short *outMinorRevision,
        mmcOSys *outSys, mmcHArch *outArch, unsigned int *outFlags,
        char *outNameStr, unsigned int *inOutNameSize,
        char *outCopyrightStr, unsigned int *inOutCopyrightSize,
        char *outCommentStr, unsigned int *inOutCommentSize) {
    VLSTACKTRACE("mmvGetVersionInfo", __FILE__, __LINE__);

    // Set version data
    if (outMajorVersion != NULL) *outMajorVersion = MEGAMOL_CORE_MAJOR_VER;
    if (outMinorVersion != NULL) *outMinorVersion = MEGAMOL_CORE_MINOR_VER;
    if (outMajorRevision != NULL) *outMajorRevision = MEGAMOL_CORE_MAJOR_REV;
    if (outMinorRevision != NULL) *outMinorRevision = MEGAMOL_CORE_MINOR_REV;

    // Set system architecture information
    if (outSys != NULL) {
        *outSys = MMC_OSYSTEM_UNKNOWN;
#ifdef _WIN32
#if defined(WINVER)
#if (WINVER >= 0x0501)
        *outSys = MMC_OSYSTEM_WINDOWS;
#endif /* (WINVER >= 0x0501) */
#endif /* defined(WINVER) */
#else /* _WIN32 */
        *outSys = MMC_OSYSTEM_LINUX;
#endif /* _WIN32 */
    }
    if (outArch != NULL) {
        *outArch = MMC_HARCH_UNKNOWN;
#if defined(_WIN64) || defined(_LIN64)
        *outArch = MMC_HARCH_X64;
#else /* defined(_WIN64) || defined(_LIN64) */
        *outArch = MMC_HARCH_I86;
#endif /* defined(_WIN64) || defined(_LIN64) */
    }

    // Set build flags
    if (outFlags != NULL) {
        *outFlags = 0
#if defined(_DEBUG) || defined(DEBUG)
            | MMC_BFLAG_DEBUG
#endif /* defined(_DEBUG) || defined(DEBUG) */
#ifdef MEGAMOL_CORE_ISDIRTY
            | MMC_BFLAG_DIRTY
#endif /* MEGAMOL_GLUT_ISDIRTY */
            ;
    }

    // Set library name
    if (inOutNameSize != NULL) {
        SIZE_T length = vislib::CharTraitsA::SafeStringLength(MEGAMOL_CORE_NAME);
        if (outNameStr != NULL) {
            if (*inOutNameSize < static_cast<unsigned int>(length)) {
                length = static_cast<SIZE_T>(*inOutNameSize);
            }
            ::memcpy(outNameStr, MEGAMOL_CORE_NAME, length);
        }
        *inOutNameSize = static_cast<unsigned int>(length);
    }

    // Set library copyright
    if (inOutCopyrightSize != NULL) {
        SIZE_T length = vislib::CharTraitsA::SafeStringLength(MEGAMOL_CORE_COPYRIGHT);
        if (outCopyrightStr != NULL) {
            if (*inOutCopyrightSize < static_cast<unsigned int>(length)) {
                length = static_cast<SIZE_T>(*inOutCopyrightSize);
            }
            ::memcpy(outCopyrightStr, MEGAMOL_CORE_COPYRIGHT, length);
        }
        *inOutCopyrightSize = static_cast<unsigned int>(length);
    }

    // Set library comments
    if (inOutCommentSize != NULL) {
        SIZE_T length = vislib::CharTraitsA::SafeStringLength(MEGAMOL_CORE_COMMENTS);
        if (outCommentStr != NULL) {
            if (*inOutCommentSize < static_cast<unsigned int>(length)) {
                length = static_cast<SIZE_T>(*inOutCommentSize);
            }
            ::memcpy(outCommentStr, MEGAMOL_CORE_COMMENTS, length);
        }
        *inOutCommentSize = static_cast<unsigned int>(length);
    }
}


/*
 * mmcGetHandleSize
 */
MEGAMOLCORE_API unsigned int MEGAMOLCORE_CALL mmcGetHandleSize(void) {
    VLSTACKTRACE("mmcGetHandleSize", __FILE__, __LINE__);
    return megamol::core::ApiHandle::GetHandleSize();
}


/*
 * mmcDisposeHandle
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcDisposeHandle(void *hndl) {
    VLSTACKTRACE("mmcDisposeHandle", __FILE__, __LINE__);
    megamol::core::ApiHandle::DestroyHandle(hndl);
}


/*
 * mmcIsHandleValid
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcIsHandleValid(void *hndl) {
    VLSTACKTRACE("mmcIsHandleValid", __FILE__, __LINE__);
    return (megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ApiHandle>(hndl) != NULL);
}


/*
 * mmcGetHandleType
 */
MEGAMOLCORE_API mmcHandleType MEGAMOLCORE_CALL mmcGetHandleType(void *hndl) {
    VLSTACKTRACE("mmcGetHandleType", __FILE__, __LINE__);

    if (megamol::core::ApiHandle::InterpretHandle<
            megamol::core::ApiHandle>(hndl) == NULL) {
        return MMC_HTYPE_INVALID;

    } else if (megamol::core::ApiHandle::InterpretHandle<
            megamol::core::CoreInstance>(hndl) == NULL) {
        return MMC_HTYPE_COREINSTANCE;

    } else if (megamol::core::ApiHandle::InterpretHandle<
            megamol::core::ViewInstance>(hndl) == NULL) {
        return MMC_HTYPE_VIEWINSTANCE;

    } else if (megamol::core::ApiHandle::InterpretHandle<
            megamol::core::JobInstance>(hndl) == NULL) {
        return MMC_HTYPE_JOBINSTANCE;

    } else if (megamol::core::ApiHandle::InterpretHandle<
            megamol::core::param::ParamHandle>(hndl) == NULL) {
        return MMC_HTYPE_PARAMETER;

    } else {
        return MMC_HTYPE_UNKNOWN;
    }
}


/*
 * mmcCreateCoreInstance
 */
MEGAMOLCORE_API mmcErrorCode MEGAMOLCORE_CALL mmcCreateCore(void *hCore) {
    VLSTACKTRACE("mmcCreateCore", __FILE__, __LINE__);
    if (mmcIsHandleValid(hCore) != 0) {
        return MMC_ERR_HANDLE; // handle was already valid.
    }
    if (*static_cast<unsigned char*>(hCore) != 0) {
        return MMC_ERR_MEMORY; // memory pointer seams to be invalid.
    }

//    { // self test for licencing
//        vislib::MD5HashProvider hash;
//        void *apifunctions[] = {
//__ALL_API_FUNCS
//        };
//        void *rp = function_cast<void*>(mmcCreateCore);
//        SIZE_T r = reinterpret_cast<SIZE_T>(rp);
//        SIZE_T d = UINT_MAX;
//        for (unsigned int i = 0; apifunctions[i] != NULL; i++) {
//            SIZE_T t = r - reinterpret_cast<SIZE_T>(apifunctions[i]);
//            if ((t > 0) && (t < d)) d = t;
//        }
//#ifdef _LOG_CORE_HASH_INFO
//        UINT logLevel = vislib::sys::Log::DefaultLog.GetLevel();
//        vislib::sys::Log::DefaultLog.SetLevel(800);
//        vislib::sys::Log::DefaultLog.WriteMsg(800, "Calculating core Hash using %d bytes\n", d);
//        vislib::sys::Log::DefaultLog.SetLevel(logLevel);
//#endif /* _LOG_CORE_HASH_INFO */
//        hash.Initialise();
//        hash.ComputeHash(NULL, r, static_cast<BYTE*>(rp), d);
//        BYTE *hashVal = new BYTE[r];
//        SIZE_T otherr = r;
//        hash.GetHashValue(hashVal, otherr);
//        ASSERT(otherr == r);
//
//        // Test against manifest hash value
//
//        vislib::StringA s, tmp;
//        for (d = 0; d < r; d++) {
//            if ((d % 4) == 0) tmp.Append("-");
//            s.Format("%02x", /*(int)*/hashVal[d]);
//            tmp.Append(s);
//        }
//        delete[] hashVal;
//        s.Format("%s-%s%d%s%s",
//            vislib::sys::SystemInformation::ComputerNameA().PeekBuffer(),
//#ifdef _WIN32
//#if defined(WINVER)
//#if (WINVER >= 0x0501)
//            "Win",
//#endif /* (WINVER >= 0x0501) */
//#endif /* defined(WINVER) */
//#else /* _WIN32 */
//            "Lin",
//#endif /* _WIN32 */
//#if defined(_WIN64) || defined(_LIN64)
//            64,
//#else /* defined(_WIN64) || defined(_LIN64) */
//            32,
//#endif /* defined(_WIN64) || defined(_LIN64) */
//#if defined(_DEBUG) || defined(DEBUG)
//            "d",
//#else /* defined(_DEBUG) || defined(DEBUG) */
//            "",
//#endif /* defined(_DEBUG) || defined(DEBUG) */
//            tmp.PeekBuffer());
//        s.ToLowerCase();
//        tmp.ToLowerCase();
//
//#ifdef _LOG_CORE_HASH_INFO
//        logLevel = vislib::sys::Log::DefaultLog.GetLevel();
//        vislib::sys::Log::DefaultLog.SetLevel(800);
//        vislib::sys::Log::DefaultLog.WriteMsg(800, "Core Hash: %s\n", tmp.PeekBuffer() + 1);
//        vislib::sys::Log::DefaultLog.SetLevel(logLevel);
//#endif /* _LOG_CORE_HASH_INFO */
//
//#ifdef _SEND_CORE_HASH_INFO
//        // send infos
//        unsigned short s1, s2, s3, s4;
//        ::mmcGetVersion(&s1, &s2, &s3, &s4);
//
//        tmp.Format("-%d-%d-%d-%d", s1, s2, s3, s4);
//        s.Append(tmp);
//
//        vislib::net::Socket::Startup();
//        try {
//            vislib::net::Socket socket;
//            socket.Create(vislib::net::Socket::FAMILY_INET, vislib::net::Socket::TYPE_STREAM, vislib::net::Socket::PROTOCOL_TCP);
//            try {
//                vislib::net::IPEndPoint endPoint = vislib::net::IPEndPoint::CreateIPv4("www.vis.uni-stuttgart.de", 80);
//                socket.Connect(endPoint);
//                try {
//                    tmp = "GET /~grottel/megamol/corehashreg.php?hash=";
//                    tmp.Append(s);
//                    // socket.SetSndTimeo(5000); // Does not work under linux, why?
//                    socket.Send(tmp.PeekBuffer(), tmp.Length());
//                } catch(...) {
//                }
//                socket.Shutdown();
//            } catch(...) {
//            }
//            socket.Close();
//        } catch(...) {
//        }
//        vislib::net::Socket::Cleanup();
//
//#endif /* _SEND_CORE_HASH_INFO */
//
//    }

#if !(defined(DEBUG) || defined(_DEBUG))
    vislib::Trace::GetInstance().SetLevel(vislib::Trace::LEVEL_VL);
#endif /* !(defined(DEBUG) || defined(_DEBUG)) */
    megamol::core::CoreInstance *inst = new megamol::core::CoreInstance();
    if (inst == NULL) {
        return MMC_ERR_MEMORY; // out of memory or initialisation failed.
    }

    return megamol::core::ApiHandle::CreateHandle(hCore, inst) 
        ? MMC_ERR_NO_ERROR : MMC_ERR_UNKNOWN;
}


/*
 * mmcSetInitialisationValue
 */
MEGAMOLCORE_API mmcErrorCode
MEGAMOLCORE_CALL mmcSetInitialisationValue(void *hCore, mmcInitValue key, 
        mmcValueType type, const void* value) {
    VLSTACKTRACE("mmcSetInitialisationValue", __FILE__, __LINE__);
    megamol::core::CoreInstance *inst
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);

    if (key == MMC_INITVAL_VISLIB_STACKTRACEMANAGER) {
        return vislib::sys::ThreadSafeStackTrace::Initialise(
            *static_cast<const vislib::SmartPtr<vislib::StackTrace>*>(value), true)
            ? MMC_ERR_NO_ERROR : MMC_ERR_UNKNOWN;
    }

    if (inst == NULL) { return MMC_ERR_INVALID_HANDLE; }
    try {

        return inst->SetInitValue(key, type, value);

    } catch(vislib::IllegalStateException) {
        return MMC_ERR_STATE;
    } catch(...) {
        return MMC_ERR_UNKNOWN;
    }

    return MMC_ERR_UNKNOWN;
}


/*
 * mmcInitialiseCoreInstance
 */
MEGAMOLCORE_API mmcErrorCode
MEGAMOLCORE_CALL mmcInitialiseCoreInstance(void *hCore) {
    VLSTACKTRACE("mmcInitialiseCoreInstance", __FILE__, __LINE__);
    megamol::core::CoreInstance *inst
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (inst == NULL) { return MMC_ERR_INVALID_HANDLE; }
    try {
        inst->Initialise();
        return MMC_ERR_NO_ERROR;
    } catch(vislib::Exception ex) {
        VLTRACE(VISLIB_TRCELVL_ERROR,
            "Failed to initialise core instance: %s (%s; %i)\n",
            ex.GetMsgA(), ex.GetFile(), ex.GetLine());
    } catch (...) {
        VLTRACE(VISLIB_TRCELVL_ERROR,
            "Failed to initialise core instance: %s\n", "unknown exception");
    }
    return MMC_ERR_UNKNOWN;
}


/*
 * mmcGetConfigurationValueA
 */
MEGAMOLCORE_API const void * MEGAMOLCORE_CALL mmcGetConfigurationValueA(
        void *hCore, mmcConfigID id, const char *name, 
        mmcValueType *outType) {
    VLSTACKTRACE("mmcGetConfigurationValueA", __FILE__, __LINE__);
    megamol::core::CoreInstance *inst
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (inst == NULL) { return NULL; }
    return inst->Configuration().GetValue(id, name, outType);
}


/*
 * mmcGetConfigurationValueW
 */
MEGAMOLCORE_API const void * MEGAMOLCORE_CALL mmcGetConfigurationValueW(
        void *hCore, mmcConfigID id, const wchar_t *name, 
        mmcValueType *outType) {
    VLSTACKTRACE("mmcGetConfigurationValueW", __FILE__, __LINE__);
    megamol::core::CoreInstance *inst
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (inst == NULL) { return NULL; }
    return inst->Configuration().GetValue(id, name, outType);
}


/*
 * mmcRequestInstanceA
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcRequestInstanceA(
        void *hCore, const char *name, const char *id) {
    VLSTACKTRACE("mmcRequestInstanceA", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;

    megamol::core::ViewDescription *vd = core->FindViewDescription(name);
    if (vd != NULL) {
        core->RequestViewInstantiation(vd, vislib::StringA(id));
        return;
    }
    megamol::core::JobDescription *jd = core->FindJobDescription(name);
    if (jd != NULL) {
        core->RequestJobInstantiation(jd, vislib::StringA(id));
        return;
    }

    vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_WARN,
        "Unable to queue instantiation of \"%s\": "
        "Description \"%s\" has not been found.", id, name);
}


/*
 * mmcRequestInstanceW
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcRequestInstanceW(
        void *hCore, const wchar_t *name, const wchar_t *id) {
    VLSTACKTRACE("mmcRequestInstanceW", __FILE__, __LINE__);
    mmcRequestInstanceA(hCore, vislib::StringA(name).PeekBuffer(),
        vislib::StringA(id).PeekBuffer());
}


/*
 * mmcHasPendingViewInstantiationRequests
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcHasPendingViewInstantiationRequests(
        void *hCore) {
    VLSTACKTRACE("mmcHasPendingViewInstantiationRequests", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return false;
    return core->HasPendingViewInstantiationRequests();
}


/*
 * mmcInstantiatePendingView
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcInstantiatePendingView(void *hCore,
        void *hView) {
    VLSTACKTRACE("mmcInstantiatePendingView", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return false;

    if (megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView) != NULL) return false;

    megamol::core::ViewInstance *view = core->InstantiatePendingView();
    if (view == NULL) return false;

    if (megamol::core::ApiHandle::CreateHandle(hView, view)) {
        megamol::core::ApiHandle::SetDeallocator(hView, core,
            megamol::core::CoreInstance::ViewJobHandleDalloc);
        return true;
    }
    return false;
}


/*
 * mmcHasPendingJobInstantiationRequests
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcHasPendingJobInstantiationRequests(
        void *hCore) {
    VLSTACKTRACE("mmcHasPendingJobInstantiationRequests", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return false;
    return core->HasPendingJobInstantiationRequests();
}


/*
 * mmcInstantiatePendingJob
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcInstantiatePendingJob(void *hCore,
        void *hJob) {
    VLSTACKTRACE("mmcInstantiatePendingJob", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return false;

    if (megamol::core::ApiHandle::InterpretHandle<
        megamol::core::JobInstance>(hJob) != NULL) return false;

    megamol::core::JobInstance *job = core->InstantiatePendingJob();
    if (job == NULL) return false;

    if (megamol::core::ApiHandle::CreateHandle(hJob, job)) {
        megamol::core::ApiHandle::SetDeallocator(hJob, core,
            megamol::core::CoreInstance::ViewJobHandleDalloc);
        return true;
    }
    return false;
}


/*
 * mmcRenderView
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcRenderView(void *hView,
        mmcRenderViewContext *context) {
    VLSTACKTRACE("mmcRenderView", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    ASSERT(context != NULL);
    ASSERT(sizeof(mmcRenderViewContext) == context->Size);

    if (view != NULL) {
        view->ModuleGraphLock().LockExclusive();

#ifdef MEGAMOLCORE_WITH_DIRECT3D11
        /* Pass in the D3D device that we created in the Viewer DLL. */
        megamol::core::view::ViewDirect3D *vd3d 
            = dynamic_cast<megamol::core::view::ViewDirect3D *>(view->View());
        if (vd3d != NULL) {
            ASSERT(context->Direct3DDevice != NULL);
            vd3d->UpdateFromContext(context);
        }
#endif /* MEGAMOLCORE_WITH_DIRECT3D11 */

        megamol::core::view::AbstractTileView *atv
            = dynamic_cast<megamol::core::view::AbstractTileView *>(view->View());
        if (atv != NULL) {
            atv->AdjustTileFromContext(context);
        }

        if (view->View() != NULL) {
            double it = context->SynchronisedTime;

            if ((it < 0.0)or(vislib::math::IsEqual(it, 0.0))) {
                // If we did not get a time via the context, determine the time
                // by using the standard method.
                it = view->View()->GetCoreInstance()->GetCoreInstanceTime();

                if (context->SynchronisedTime != 0) {
                    // The viewer module wants to reuse this time until it
                    // resets 'SynchronisedTime' to -1.
                    context->SynchronisedTime = it;
                }
            }
            view->View()->Render(view->View()->DefaultTime(it), it);
            context->ContinuousRedraw = true; // TODO: Implement the real thing
        }
        view->ModuleGraphLock().UnlockExclusive();
    }
}


/*
 * mmcRegisterViewCloseRequestFunction
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcRegisterViewCloseRequestFunction(
        void *hView, mmcViewCloseRequestFunction func, void *data) {
    VLSTACKTRACE("mmcRegisterViewCloseRequestFunction", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    if (view != NULL) {
        view->SetCloseRequestCallback(func, data);
    }
}


/*
 * mmcResizeView
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcResizeView(void *hView,
        unsigned int width, unsigned int height) {
    VLSTACKTRACE("mmcResizeView", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    if ((view != NULL) && (view->View() != NULL)) {
        view->View()->Resize(width, height);
    }
}


/*
 * mmcSetInputModifier
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcSetInputModifier(void *hView,
        mmcInputModifier mod, bool down) {
    VLSTACKTRACE("mmcSetInputModifier", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    if ((view != NULL) && (view->View() != NULL)) {
        view->View()->SetInputModifier(mod, down);
    }
}


/*
 * mmcSet2DMouseButton
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcSet2DMouseButton(void *hView,
        unsigned int btn, bool down) {
    VLSTACKTRACE("mmcSet2DMouseButton", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    if ((view != NULL) && (view->View() != NULL)) {
        view->View()->SetCursor2DButtonState(btn, down);
    }
}


/*
 * mmcSet2DMousePosition
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcSet2DMousePosition(void *hView,
        float x, float y) {
    VLSTACKTRACE("mmcSet2DMousePosition", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    if ((view != NULL) && (view->View() != NULL)) {
        view->View()->SetCursor2DPosition(x, y);
    }
}


/*
 * mmcDesiredViewWindowConfig
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcDesiredViewWindowConfig(void *hView,
        int *x, int *y, int *w, int *h, bool *nd) {
    VLSTACKTRACE("mmcDesiredViewWindowConfig", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    if ((view != NULL) && (view->View() != NULL)) {
        return view->View()->DesiredWindowPosition(x, y, w, h, nd);
    }
    return false;
}


/*
 * mmcIsJobRunning
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcIsJobRunning(void *hJob) {
    VLSTACKTRACE("mmcIsJobRunning", __FILE__, __LINE__);
    megamol::core::JobInstance *job
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::JobInstance>(hJob);
    if ((job != NULL) && (job->Job() != NULL)) {
        return job->Job()->IsRunning();
    }
    return false;
}


/*
 * mmcIsViewRunning
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcIsViewRunning(void *hView) {
    VLSTACKTRACE("mmcIsViewRunning", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    return (view != NULL) && (view->View() != NULL);
}


/*
 * mmcStartJob
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcStartJob(void *hJob) {
    VLSTACKTRACE("mmcStartJob", __FILE__, __LINE__);
    megamol::core::JobInstance *job
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::JobInstance>(hJob);
    if ((job != NULL) && (job->Job() != NULL)) {
        return job->Job()->Start();
    }
    return false;
}


/*
 * mmcTerminateJob
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcTerminateJob(void *hJob) {
    VLSTACKTRACE("mmcTerminateJob", __FILE__, __LINE__);
    megamol::core::JobInstance *job
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::JobInstance>(hJob);
    if ((job != NULL) && (job->Job() != NULL)) {
        job->Job()->Terminate();
    }
}


/*
 * mmcSetParameterA
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcSetParameterValueA(void *hParam,
        const char *value) {
    VLSTACKTRACE("mmcSetParameterValueA", __FILE__, __LINE__);
    megamol::core::param::ParamHandle *param
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::param::ParamHandle>(hParam);
    if (param == NULL) return;
    if (param->GetParameter()->ParseValue(A2T(value))) {
        vislib::StringA name;
        param->GetIDString(name);
        // TODO: Change text if it is a button parameter
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_INFO,
            "Setting parameter \"%s\" to \"%s\"",
            name.PeekBuffer(), vislib::StringA(
            param->GetParameter()->ValueString()).PeekBuffer());
    } else {
        vislib::StringA name;
        param->GetIDString(name);
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_ERROR,
            "Unable to set parameter \"%s\": Failed to parse value \"%s\"",
            name.PeekBuffer(), value);
    }
}


/*
 * mmcSetParameterW
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcSetParameterValueW(void *hParam,
        const wchar_t *value) {
    VLSTACKTRACE("mmcSetParameterValueW", __FILE__, __LINE__);
    megamol::core::param::ParamHandle *param
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::param::ParamHandle>(hParam);
    if (param == NULL) return;
    if (param->GetParameter()->ParseValue(W2T(value))) {
        vislib::StringA name;
        param->GetIDString(name);
        // TODO: Change text if it is a button parameter
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_INFO,
            "Setting parameter \"%s\" to \"%s\"",
            name.PeekBuffer(), vislib::StringA(
            param->GetParameter()->ValueString()).PeekBuffer());
    } else {
        vislib::StringA name;
        param->GetIDString(name);
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_ERROR,
            "Unable to set parameter \"%s\": Failed to parse value \"%s\"",
            name.PeekBuffer(), vislib::StringA(value).PeekBuffer());
    }
}


/*
 * mmcLoadProjectA
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcLoadProjectA(void *hCore,
        const char *filename) {
    VLSTACKTRACE("mmcLoadProjectA", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;
    core->LoadProject(filename);
}


/*
 * mmcLoadProjectW
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcLoadProjectW(void *hCore,
        const wchar_t *filename) {
    VLSTACKTRACE("mmcLoadProjectW", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;
    core->LoadProject(filename);
}


/*
 * mmcGetParameterA
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcGetParameterA(void *hCore,
        const char *name, void *hParam, bool bCreate) {
    VLSTACKTRACE("mmcGetParameterA", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return false;

    vislib::SmartPtr<megamol::core::param::AbstractParam>
        param = core->FindParameter(name, false, bCreate);
    if (param.IsNull()) return false;

    if (mmcIsHandleValid(hParam) != 0) {
        return false; // handle was already valid.
    }
    if (*static_cast<unsigned char*>(hParam) != 0) {
        return false; // memory pointer seams to be invalid.
    }

    return megamol::core::ApiHandle::CreateHandle(hParam,
        new megamol::core::param::ParamHandle(*core, param));
}


/*
 * mmcGetParameterW
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcGetParameterW(void *hCore,
        const wchar_t *name, void *hParam, bool bCreate) {
    VLSTACKTRACE("mmcGetParameterW", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return false;

    vislib::SmartPtr<megamol::core::param::AbstractParam>
        param = core->FindParameter(name, false, bCreate);
    if (param.IsNull()) return false;

    if (mmcIsHandleValid(hParam) != 0) {
        return false; // handle was already valid.
    }
    if (*static_cast<unsigned char*>(hParam) != 0) {
        return false; // memory pointer seams to be invalid.
    }

    return megamol::core::ApiHandle::CreateHandle(hParam,
        new megamol::core::param::ParamHandle(*core, param));
}


/*
 * mmcGetParameterValueA
 */
MEGAMOLCORE_API const char * MEGAMOLCORE_CALL mmcGetParameterValueA(
        void *hParam) {
    static vislib::StringA retval;
    VLSTACKTRACE("mmcGetParameterValueA", __FILE__, __LINE__);
    megamol::core::param::ParamHandle *param
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::param::ParamHandle>(hParam);
    if (param == NULL) return NULL;
    retval = param->GetParameter()->ValueString();
    return retval.PeekBuffer();

}


/*
 * mmcGetParameterValueW
 */
MEGAMOLCORE_API const wchar_t * MEGAMOLCORE_CALL mmcGetParameterValueW(
        void *hParam) {
    static vislib::StringW retval;
    VLSTACKTRACE("mmcGetParameterValueW", __FILE__, __LINE__);
    megamol::core::param::ParamHandle *param
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::param::ParamHandle>(hParam);
    if (param == NULL) return NULL;
    retval = param->GetParameter()->ValueString();
    return retval.PeekBuffer();
}


/*
 * mmcEnumParametersA
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcEnumParametersA(void *hCore,
        mmcEnumStringAFunction func, void *data) {
    VLSTACKTRACE("mmcEnumParametersA", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;
    core->EnumParameters(func, data);
}


/**
 * Helper struct for unicode parameter enumeration
 */
typedef struct _EnumParamDataHelperW_t {
    mmcEnumStringWFunction func;
    void *data;
} EnumParamDataHelperW;


extern "C" {

/**
 * Helper function for unicode parameter enumeration
 *
 * @param name The parameter name
 * @param data The enumeration data
 */
static void
#ifdef _WIN32
__stdcall
#endif /* _WIN32 */
EnumParamsW(const char *name, void *data) {
    EnumParamDataHelperW *context
        = reinterpret_cast<EnumParamDataHelperW*>(data);
    context->func(vislib::StringW(name).PeekBuffer(), context->data);
}

}


/*
 * mmcEnumParametersW
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcEnumParametersW(void *hCore,
        mmcEnumStringWFunction func, void *data) {
    VLSTACKTRACE("mmcEnumParametersW", __FILE__, __LINE__);
    EnumParamDataHelperW context;
    context.func = func;
    context.data = data;
    ::mmcEnumParametersA(hCore, ::EnumParamsW, &context);
}


/*
 * mmcGetInstanceIDA
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcGetInstanceIDA(void *hInst,
        char *buf, unsigned int *len) {
    VLSTACKTRACE("mmcGetInstanceIDA", __FILE__, __LINE__);
    if (len == NULL) return;

    vislib::StringA id;

    megamol::core::ViewInstance *vi
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hInst);
    megamol::core::JobInstance *ji
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::JobInstance>(hInst);
    megamol::core::param::ParamHandle *ph
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::param::ParamHandle>(hInst);

    if (vi != NULL) {
        id = vi->Name();
    } else if (ji != NULL) {
        id = ji->Name();
    } else if (ph != NULL){
        ph->GetIDString(id);
    }

    if (buf == NULL) {
        *len = id.Length() + 1;
    } else {
        memcpy(buf, id.PeekBuffer(),
            vislib::math::Min<SIZE_T>(id.Length() + 1, *len));
    }
}


/*
 * mmcGetInstanceIDW
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcGetInstanceIDW(void *hInst,
        wchar_t *buf, unsigned int *len) {
    VLSTACKTRACE("mmcGetInstanceIDW", __FILE__, __LINE__);
    if (len == NULL) return;

    char *bufA = new char[*len];
    ::mmcGetInstanceIDA(hInst, bufA, len);
    vislib::StringW id(vislib::StringA(bufA, *len));
    delete[] bufA;

    if (buf == NULL) {
        *len = id.Length() + 1;
    } else {
        memcpy(buf, id.PeekBuffer(), sizeof(wchar_t)
            * vislib::math::Min<SIZE_T>(id.Length() + 1, *len));
    }
}


/*
 * mmcIsParameterRelevant
 */
MEGAMOLCORE_API bool MEGAMOLCORE_CALL mmcIsParameterRelevant(void *hInst,
        void *hParam) {
    VLSTACKTRACE("mmcIsParameterRelevant", __FILE__, __LINE__);
    megamol::core::param::ParamHandle *param
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::param::ParamHandle>(hParam);
    if (param == NULL) return false;

    megamol::core::JobInstance *ji
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::JobInstance>(hInst);
    if ((ji != NULL) && (ji->Job() != NULL)) {
        return ji->Job()->IsParamRelevant(param->GetParameter());
    }

    megamol::core::ViewInstance *vi
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hInst);
    if ((vi != NULL) && (vi->View() != NULL)) {
        return vi->View()->IsParamRelevant(param->GetParameter());
    }

    return false;
}


/*
 * mmcGetParameterTypeDescription
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcGetParameterTypeDescription(
        void *hParam, unsigned char *buf, unsigned int *len) {
    VLSTACKTRACE("mmcGetParameterTypeDescription", __FILE__, __LINE__);
    megamol::core::param::ParamHandle *param
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::param::ParamHandle>(hParam);

    if (len == NULL) return;
    if (param != NULL) {
        vislib::RawStorage rs;
        param->GetParameter()->Definition(rs);
        if (buf != NULL) {
            unsigned int s = vislib::math::Min<unsigned int>(
                static_cast<unsigned int>(rs.GetSize()), *len);
            ::memcpy(buf, rs.As<unsigned char>(), s);
            *len = s;
        } else {
            *len = static_cast<unsigned int>(rs.GetSize());
        }
    } else {
        *len = 0;
    }
}


/*
 * mmcFreezeOrUpdateView
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcFreezeOrUpdateView(
        void *hView, bool freeze) {
    VLSTACKTRACE("mmcFreezeOrUpdateView", __FILE__, __LINE__);
    megamol::core::ViewInstance *view
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::ViewInstance>(hView);
    if ((view != NULL) && (view->View() != NULL)) {
        view->View()->UpdateFreeze(freeze);
    }
}


/*
 * mmcQuickstartA
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcQuickstartA(void *hCore, const char *filename) {
    VLSTACKTRACE("mmcQuickstartA", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;
    core->Quickstart(A2T(filename));
}


/*
 * mmcQuickstartW
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcQuickstartW(void *hCore, const wchar_t *filename) {
    VLSTACKTRACE("mmcQuickstartW", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;
    core->Quickstart(W2T(filename));
}


/*
 * mmcQuickstartRegistryA
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcQuickstartRegistryA(void *hCore,
        const char *frontend, const char *feparams,
        const char *filetype, bool unreg, bool overwrite) {
    VLSTACKTRACE("mmcQuickstartRegistryA", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;
    core->QuickstartRegistry(A2T(frontend), A2T(feparams), A2T(filetype), unreg, overwrite);
}


/*
 * mmcQuickstartRegistryW
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcQuickstartRegistryW(void *hCore,
        const wchar_t *frontend, const wchar_t *feparams,
        const wchar_t *filetype, bool unreg, bool overwrite) {
    VLSTACKTRACE("mmcQuickstartRegistryW", __FILE__, __LINE__);
    megamol::core::CoreInstance *core
        = megamol::core::ApiHandle::InterpretHandle<
        megamol::core::CoreInstance>(hCore);
    if (core == NULL) return;
    core->QuickstartRegistry(W2T(frontend), W2T(feparams), W2T(filetype), unreg, overwrite);
}


/*
 * mmcModuleCount
 */
MEGAMOLCORE_API int MEGAMOLCORE_CALL mmcModuleCount(void) {
    int cnt = 0;
    megamol::core::ModuleDescriptionManager::DescriptionIterator it
        = megamol::core::ModuleDescriptionManager::Instance()->GetIterator();
    while (it.HasNext()) {
        it.Next();
        cnt++;
    }
    return cnt;
}


/*
 * mmcModuleDescription
 */
MEGAMOLCORE_API void* MEGAMOLCORE_CALL mmcModuleDescription(int idx) {
    if (idx < 0) return nullptr;
    megamol::core::ModuleDescriptionManager::DescriptionIterator it
        = megamol::core::ModuleDescriptionManager::Instance()->GetIterator();
    while (it.HasNext()) {
        megamol::core::ModuleDescription *d = it.Next();
        if (idx == 0) return d;
        idx--;
    }
    return nullptr;
}


/*
 * mmcCallCount
 */
MEGAMOLCORE_API int MEGAMOLCORE_CALL mmcCallCount(void) {
    int cnt = 0;
    megamol::core::CallDescriptionManager::DescriptionIterator it
        = megamol::core::CallDescriptionManager::Instance()->GetIterator();
    while (it.HasNext()) {
        it.Next();
        cnt++;
    }
    return cnt;
}


/*
 * mmcCallDescription
 */
MEGAMOLCORE_API void* MEGAMOLCORE_CALL mmcCallDescription(int idx) {
    if (idx < 0) return nullptr;
    megamol::core::CallDescriptionManager::DescriptionIterator it
        = megamol::core::CallDescriptionManager::Instance()->GetIterator();
    while (it.HasNext()) {
        megamol::core::CallDescription *d = it.Next();
        if (idx == 0) return d;
        idx--;
    }
    return nullptr;
}


/**
 * TODO: Document me
 */
bool operator==(const mmcParamSlotDescription lhs, mmcParamSlotDescription rhs) {
    return ::memcmp(&lhs, &rhs, sizeof(mmcParamSlotDescription)) == 0;
}


/**
 * TODO: Document me
 */
bool operator==(const mmcCalleeSlotDescription lhs, mmcCalleeSlotDescription rhs) {
    return ::memcmp(&lhs, &rhs, sizeof(mmcCalleeSlotDescription)) == 0;
}


/**
 * TODO: Document me
 */
bool operator==(const mmcCallerSlotDescription lhs, mmcCallerSlotDescription rhs) {
    return ::memcmp(&lhs, &rhs, sizeof(mmcCallerSlotDescription)) == 0;
}


/*
 * mmcGetModuleSlotDescriptions
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcGetModuleSlotDescriptions(void * desc, 
        unsigned int *outCntParamSlots, mmcParamSlotDescription **outParamSlots,
        unsigned int *outCntCalleeSlots, mmcCalleeSlotDescription **outCalleeSlots,
        unsigned int *outCntCallerSlots, mmcCallerSlotDescription **outCallerSlots) {
    ASSERT(desc != NULL);
    ASSERT(outCntParamSlots != NULL);
    ASSERT(outParamSlots != NULL);
    ASSERT(outCntCalleeSlots != NULL);
    ASSERT(outCalleeSlots != NULL);
    ASSERT(outCntCallerSlots != NULL);
    ASSERT(outCallerSlots != NULL);

    megamol::core::RootModuleNamespace rms;
    megamol::core::ModuleDescription *md = static_cast<megamol::core::ModuleDescription*>(desc);
    ASSERT(md != NULL);

    megamol::core::Module *m = md->CreateModule(NULL, NULL);
    if (m == NULL) {
        *outCntParamSlots = 0;
        *outParamSlots = NULL;
        *outCntCalleeSlots = 0;
        *outCalleeSlots = NULL;
        *outCntCallerSlots = 0;
        *outCallerSlots = NULL;
        return;
    }
    rms.AddChild(m);

    vislib::Array<mmcParamSlotDescription> pa;
    vislib::Array<mmcCalleeSlotDescription> cea;
    vislib::Array<mmcCallerSlotDescription> cra;

    vislib::Stack<megamol::core::Module::ChildList::Iterator> stack;
    stack.Push(m->GetChildIterator());

    while (!stack.IsEmpty()) {
        megamol::core::Module::ChildList::Iterator iter = stack.Pop();
        while (iter.HasNext()) {
            megamol::core::AbstractNamedObject *ano = iter.Next();

            megamol::core::param::ParamSlot *ps = dynamic_cast<megamol::core::param::ParamSlot*>(ano);
            megamol::core::CalleeSlot *ces = dynamic_cast<megamol::core::CalleeSlot*>(ano);
            megamol::core::CallerSlot *crs = dynamic_cast<megamol::core::CallerSlot*>(ano);

            if (ps != NULL) {
                SIZE_T i = pa.Count();
                pa.SetCount(i + 1);

                vislib::StringA str = ps->FullName();
                vislib::StringA str2 = m->FullName() + "::";
                if (str.StartsWith(str2)) str.Remove(0, str2.Length());

                pa[i].name = new char[str.Length() + 1];
                ::memcpy(const_cast<char*>(pa[i].name), str.PeekBuffer(), str.Length() + 1);

                str = ps->Description();
                pa[i].desc = new char[str.Length() + 1];
                ::memcpy(const_cast<char*>(pa[i].desc), str.PeekBuffer(), str.Length() + 1);

                vislib::RawStorage blob;
                ps->Param< ::megamol::core::param::AbstractParam>()->Definition(blob);

                pa[i].typeInfoSize = static_cast<unsigned int>(blob.GetSize());
                pa[i].typeInfo = new unsigned char[pa[i].typeInfoSize];
                ::memcpy(const_cast<unsigned char*>(pa[i].typeInfo), blob, pa[i].typeInfoSize);

                str = ps->Parameter()->ValueString();
                pa[i].defVal = new char[str.Length() + 1];
                ::memcpy(const_cast<char*>(pa[i].defVal), str.PeekBuffer(), str.Length() + 1);

            } else if (ces != NULL) {
                SIZE_T i = cea.Count();
                cea.SetCount(i + 1);

                vislib::StringA str = ces->FullName();
                vislib::StringA str2 = m->FullName() + "::";
                if (str.StartsWith(str2)) str.Remove(0, str2.Length());

                cea[i].name = new char[str.Length() + 1];
                ::memcpy(const_cast<char*>(cea[i].name), str.PeekBuffer(), str.Length() + 1);

                str = ces->Description();
                cea[i].desc = new char[str.Length() + 1];
                ::memcpy(const_cast<char*>(cea[i].desc), str.PeekBuffer(), str.Length() + 1);

                cea[i].cntCallbacks = static_cast<unsigned int>(ces->GetCallbackCount());
                cea[i].callbackCallType = new const char*[cea[i].cntCallbacks];
                cea[i].callbackFuncName = new const char*[cea[i].cntCallbacks];

                for (unsigned int j = 0; j < cea[i].cntCallbacks; j++) {
                    str = ces->GetCallbackCallName(j);
                    cea[i].callbackCallType[j] = new char[str.Length() + 1];
                    ::memcpy(const_cast<char*>(cea[i].callbackCallType[j]), str.PeekBuffer(), str.Length() + 1);

                    str = ces->GetCallbackFuncName(j);
                    cea[i].callbackFuncName[j] = new char[str.Length() + 1];
                    ::memcpy(const_cast<char*>(cea[i].callbackFuncName[j]), str.PeekBuffer(), str.Length() + 1);
                }

            } else if (crs != NULL) {
                SIZE_T i = cra.Count();
                cra.SetCount(i + 1);

                vislib::StringA str = crs->FullName();
                vislib::StringA str2 = m->FullName() + "::";
                if (str.StartsWith(str2)) str.Remove(0, str2.Length());

                cra[i].name = new char[str.Length() + 1];
                ::memcpy(const_cast<char*>(cra[i].name), str.PeekBuffer(), str.Length() + 1);

                str = crs->Description();
                cra[i].desc = new char[str.Length() + 1];
                ::memcpy(const_cast<char*>(cra[i].desc), str.PeekBuffer(), str.Length() + 1);

                cra[i].cntCompCalls = static_cast<unsigned int>(crs->GetCompCallCount());
                cra[i].compCalls = new const char *[cra[i].cntCompCalls];
                for (unsigned int j = 0; j < cra[i].cntCompCalls; j++) {
                    str = crs->GetCompCallClassName(j);
                    cra[i].compCalls[j] = new char[str.Length() + 1];
                    ::memcpy(const_cast<char*>(cra[i].compCalls[j]), str.PeekBuffer(), str.Length() + 1);
                }

            }
        }
    }

    *outCntParamSlots = static_cast<unsigned int>(pa.Count());
    *outParamSlots = new mmcParamSlotDescription[*outCntParamSlots];
    ::memcpy(*outParamSlots, pa.PeekElements(), sizeof(mmcParamSlotDescription) * *outCntParamSlots);
    *outCntCalleeSlots = static_cast<unsigned int>(cea.Count());
    *outCalleeSlots = new mmcCalleeSlotDescription[*outCntCalleeSlots];
    ::memcpy(*outCalleeSlots, cea.PeekElements(), sizeof(mmcCalleeSlotDescription) * *outCntCalleeSlots);
    *outCntCallerSlots = static_cast<unsigned int>(cra.Count());
    *outCallerSlots = new mmcCallerSlotDescription[*outCntCallerSlots];
    ::memcpy(*outCallerSlots, cra.PeekElements(), sizeof(mmcCallerSlotDescription) * *outCntCallerSlots);

    rms.RemoveChild(m);
    m->SetAllCleanupMarks();
    m->PerformCleanup();
    delete m;
}


/*
 * mmcReleaseModuleSlotDescriptions
 */
MEGAMOLCORE_API void MEGAMOLCORE_CALL mmcReleaseModuleSlotDescriptions(
        unsigned int outCntParamSlots, mmcParamSlotDescription **outParamSlots,
        unsigned int outCntCalleeSlots, mmcCalleeSlotDescription **outCalleeSlots,
        unsigned int outCntCallerSlots, mmcCallerSlotDescription **outCallerSlots) {
    ASSERT(outParamSlots != NULL);
    ASSERT(outCalleeSlots != NULL);
    ASSERT(outCallerSlots != NULL);

    for (unsigned int i = 0; i < outCntParamSlots; i++) {
        delete[] (*outParamSlots)[i].name;
        delete[] (*outParamSlots)[i].desc;
        delete[] (*outParamSlots)[i].typeInfo;
        delete[] (*outParamSlots)[i].defVal;
    }
    delete[] (*outParamSlots);
    *outParamSlots = NULL;

    for (unsigned int i = 0; i < outCntCalleeSlots; i++) {
        delete[] (*outCalleeSlots)[i].name;
        delete[] (*outCalleeSlots)[i].desc;
        for (unsigned int j = 0; j < (*outCalleeSlots)[i].cntCallbacks; j++) {
            delete[] (*outCalleeSlots)[i].callbackCallType[j];
            delete[] (*outCalleeSlots)[i].callbackFuncName[j];
        }
        delete[] (*outCalleeSlots)[i].callbackCallType;
        delete[] (*outCalleeSlots)[i].callbackFuncName;
    }
    delete[] (*outCalleeSlots);
    *outCalleeSlots = NULL;

    for (unsigned int i = 0; i < outCntCallerSlots; i++) {
        delete[] (*outCallerSlots)[i].name;
        delete[] (*outCallerSlots)[i].desc;
        for (unsigned int j = 0; j < (*outCallerSlots)[i].cntCompCalls; j++) {
            delete[] (*outCallerSlots)[i].compCalls[j];
        }
        delete[] (*outCallerSlots)[i].compCalls;
    }
    delete[] (*outCallerSlots);
    *outCallerSlots = NULL;

}


/*
 * mmcGetModuleDescriptionInfo
 */
MEGAMOLCORE_EXT_APICALL(mmcModuleDescriptionInfo*, mmcGetModuleDescriptionInfo)(void * desc) {
    megamol::core::ModuleDescription *md = static_cast<megamol::core::ModuleDescription*>(desc);
    ASSERT(md != NULL);
    mmcModuleDescriptionInfo *d = new mmcModuleDescriptionInfo();

    vislib::StringA str = md->ClassName();
    d->name = new char[str.Length() + 1];
    ::memcpy(const_cast<char*>(d->name), str.PeekBuffer(), str.Length() + 1);

    str = md->Description();
    d->desc = new char[str.Length() + 1];
    ::memcpy(const_cast<char*>(d->desc), str.PeekBuffer(), str.Length() + 1);

    return d;
}


/*
 * mmcReleaseModuleDescriptionInfo
 */
MEGAMOLCORE_EXT_APICALL(void, mmcReleaseModuleDescriptionInfo)(mmcModuleDescriptionInfo* desc) {
    delete[] desc->name;
    delete[] desc->desc;
    delete desc;
}


/*
 * mmcGetCallDescriptionInfo
 */
MEGAMOLCORE_EXT_APICALL(mmcCallDescriptionInfo*, mmcGetCallDescriptionInfo)(void * desc) {
    megamol::core::CallDescription *cd = static_cast<megamol::core::CallDescription*>(desc);
    ASSERT(cd != NULL);
    mmcCallDescriptionInfo *d = new mmcCallDescriptionInfo();

    vislib::StringA str = cd->ClassName();
    d->name = new char[str.Length() + 1];
    ::memcpy(const_cast<char*>(d->name), str.PeekBuffer(), str.Length() + 1);

    str = cd->Description();
    d->desc = new char[str.Length() + 1];
    ::memcpy(const_cast<char*>(d->desc), str.PeekBuffer(), str.Length() + 1);

    d->cntFunc = cd->FunctionCount();
    d->funcNames = new const char*[d->cntFunc];
    for (unsigned int i = 0; i < d->cntFunc; i++) {
        str = cd->FunctionName(i);
        d->funcNames[i] = new char[str.Length() + 1];
        ::memcpy(const_cast<char*>(d->funcNames[i]), str.PeekBuffer(), str.Length() + 1);
    }

    return d;
}


/*
 * mmcReleaseCallDescriptionInfo
 */
MEGAMOLCORE_EXT_APICALL(void, mmcReleaseCallDescriptionInfo)(mmcCallDescriptionInfo* desc) {
    delete[] desc->name;
    delete[] desc->desc;
    for (unsigned int i = 0; i < desc->cntFunc; i++) {
        delete[] desc->funcNames[i];
    }
    delete[] desc->funcNames;
    delete desc;
}
