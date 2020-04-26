/*
 * ProbesToTable.cpp
 * Copyright (C) 2020 by MegaMol Team
 * Alle Rechte vorbehalten.
 */

#include "mmcore/Module.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/CalleeSlot.h"
#include "mmstd_datatools/table/TableDataCall.h"

namespace megamol {
namespace probe {

class ProbeToTable : public core::Module {
public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName(void) { return "ProbeToTable"; }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description(void) { return "Converts Probes data to table data"; }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable(void) { return true; }

    /** Ctor. */
    ProbeToTable(void);

    /** Dtor. */
    virtual ~ProbeToTable(void);

protected:
    virtual bool create();
    virtual void release();

    core::CallerSlot _getDataSlot;
    core::CalleeSlot _deployTableSlot;

private:
    bool getMetaData(core::Call& call);
    bool getData(core::Call& call);

    bool InterfaceIsDirty();
    size_t _currentFrame;
    std::vector<float> _floatBlob;
    std::vector<stdplugin::datatools::table::TableDataCall::ColumnInfo> _colinfo;
    size_t _cols = 0;
    size_t _rows = 0;
    size_t _datahash = 0;

    };
} // namespace probe
    } // namespace megamol