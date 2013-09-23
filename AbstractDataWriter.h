/*
 * AbstractDataWriter.h
 *
 * Copyright (C) 2010 by Universitaet Stuttgart (VISUS)
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOLCORE_ABSTRACTDATAWRITER_H_INCLUDED
#define MEGAMOLCORE_ABSTRACTDATAWRITER_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "Module.h"
#include "CalleeSlot.h"
#include "DataWriterCtrlCall.h"


namespace megamol {
namespace core {

    /**
     * Abstract base class for data writer modules
     */
    class MEGAMOLCORE_API AbstractDataWriter : public Module {
    public:

        /** Ctor. */
        AbstractDataWriter(void);

        /** Dtor. */
        virtual ~AbstractDataWriter(void);

    protected:

        /**
         * The main function
         *
         * @return True on success
         */
        virtual bool run(void) = 0;

        /**
         * Function querying the writers capabilities
         *
         * @param call The call to receive the capabilities
         *
         * @return True on success
         */
        virtual bool getCapabilities(DataWriterCtrlCall& call) = 0;

        /**
         * Called to abort the run function.
         *
         * @return True on success
         */
        virtual bool abort(void);

    private:

        /**
         * Event handler for incoming run calls
         *
         * @param call The incoming call
         *
         * @return True on success
         */
        bool onCallRun(Call &call);

        /**
         * Event handler for incoming run calls
         *
         * @param call The incoming call
         *
         * @return True on success
         */
        bool onCallGetCapability(Call &call);

        /**
         * Event handler for incoming run calls
         *
         * @param call The incoming call
         *
         * @return True on success
         */
        bool onCallAbort(Call &call);

        /** The slot for incoming control commands */
        CalleeSlot controlSlot;

    };

} /* end namespace core */
} /* end namespace megamol */

#endif /* MEGAMOLCORE_ABSTRACTDATAWRITER_H_INCLUDED */
