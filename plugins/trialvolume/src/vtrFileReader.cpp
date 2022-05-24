#include "vtrFileReader.h"

#include "mmcore/param/BoolParam.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/IntParam.h"

#include <iostream>
#include <fstream>

using namespace megamol;

trialvolume::vtrFileReader::vtrFileReader() : core::Module()
        , getDataCalleeSlot_("getdata", "Slot to request data from this data source.")
        , filename_("filename", "The path to the vtr file to load.")
        , dataHash(0)
        , metadata()
        , vtrFilename("")
        , fileChanged_(true) {
    // Setup filename input slot
    this->filename_.SetParameter(new core::param::FilePathParam(""));
    this->filename_.SetUpdateCallback(&trialvolume::vtrFileReader::filenameChanged);
    this->MakeSlotAvailable(&this->filename_);

    // Setup volumetric data output slot
    this->getDataCalleeSlot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_DATA),
        &trialvolume::vtrFileReader::getDataCallback);
    this->getDataCalleeSlot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &trialvolume::vtrFileReader::getDataCallback);
    this->getDataCalleeSlot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &trialvolume::vtrFileReader::getDataCallback);
    this->getDataCalleeSlot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &trialvolume::vtrFileReader::dummyCallback);
    this->getDataCalleeSlot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &trialvolume::vtrFileReader::dummyCallback);
    this->getDataCalleeSlot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &trialvolume::vtrFileReader::dummyCallback);
    this->MakeSlotAvailable(&this->getDataCalleeSlot_);
}

trialvolume::vtrFileReader::~vtrFileReader() {
    this->Release();
}

bool trialvolume::vtrFileReader::create(void) {
    return true;
}

void trialvolume::vtrFileReader::release(void) {
}

bool trialvolume::vtrFileReader::dummyCallback(core::Call& caller) {
    return true;
}

// TODO This function is probably *super* unsafe. Replace with something better.
bool trialvolume::vtrFileReader::filenameChanged(core::param::ParamSlot& slot) {
    vtrFilename = this->filename_.Param<core::param::FilePathParam>()->ValueString();
    this->fileChanged_ = true;
    return this->loadFile();
}

bool trialvolume::vtrFileReader::loadFile() {
    if (vtrFilename.empty()) {
        core::utility::log::Log::DefaultLog.WriteError("Empty vtr file %s!", vtrFilename);
        return false;
    }

    // Set general metadata which is not directly given in the vtr file
    metadata.Components = 1;
    metadata.GridType = geocalls::GridType_t::RECTILINEAR;
    metadata.ScalarType = geocalls::ScalarType_t::FLOATING_POINT;
    metadata.ScalarLength = sizeof(float);

    metadata.NumberOfFrames = 1;

    // Open file
    std::ifstream vtrFile(vtrFilename.c_str(), std::ios::binary);
    if (!vtrFile.is_open()) {
        core::utility::log::Log::DefaultLog.WriteError("Could not open vtr file %s!", vtrFilename.c_str());
        return false;
    }
    // Read header info
    // structured like following: Name, Type, Amount, Offset, Size
    // with the following entries:
    //  - Name: string "[xyz]_coordinates" or "vof-function[-]"
    //  - Type: string Float64 (for now)
    //  - Amount: int (expected to be 1)
    //  - Offset: int offset in blob in bytes
    //  - Size: int amount of bytes total
    
    // Parse header line by line
    std::string line;
    std::string name;
    std::string type;
    size_t amount;
    size_t offset;
    size_t size;

    size_t offsets[4];

    while (std::getline(vtrFile, line)) {
        // Parse line
        std::stringstream ss(line);
        ss >> name >> type >> amount >> offset >> size;
        // Check if line decomposition failed, if so we encountered the end of the headers
        if (!ss) {
            break;
        }
        // Check if correct type is given
        if (type != "Float64") {
            core::utility::log::Log::DefaultLog.WriteError("Wrong type given in vtr file %s!", vtrFilename.c_str());
            return false;
        }
        // Check if expected amount is given
        if (amount != 1) {
            core::utility::log::Log::DefaultLog.WriteError("Wrong amount given in vtr file %s!", vtrFilename.c_str());
            return false;
        }
        // Check if name is given
        if (name == "x_coordinates") {
            // TODO: Actually use the vof as volumes of the cells not as values on the grid points 
            metadata.Resolution[0] = size / 8 - 1;
            offsets[0] = offset;
        } else if (name == "y_coordinates") {
            metadata.Resolution[1] = size / 8 - 1;
            offsets[1] = offset;
        } else if (name == "z_coordinates") {
            metadata.Resolution[2] = size / 8 - 1;
            offsets[2] = offset;
        } else if (name == "vof-function[-]") {
            offsets[3] = offset;
        } else {
            core::utility::log::Log::DefaultLog.WriteError("Unknown name %s given in vtr file %s!", name, vtrFilename.c_str());
            return false;
        }
    }

    auto const base_offset = vtrFile.tellg();

    // Read all the axis data
    for (size_t axisIdx = 0; axisIdx < 3; axisIdx++) {

        vtrFile.seekg(base_offset + std::streamoff(offsets[axisIdx]), std::ios::beg);
        metadata.SliceDists[axisIdx] = new float[metadata.Resolution[axisIdx] - 1];
        char buffer[sizeof(double)];
        vtrFile.read(buffer, sizeof(double));
        auto prev_val = *reinterpret_cast<double*>(buffer);
        auto min = prev_val;
        auto max = prev_val;
        for (size_t i = 0; i < metadata.Resolution[axisIdx] - 1; i++) {
            vtrFile.read(buffer, sizeof(double));
            auto const val = *reinterpret_cast<double*>(buffer);
            metadata.SliceDists[axisIdx][i] = static_cast<float>(val - prev_val);
            prev_val = val;
            if (val < min) {
                min = val;
            } else if (val > max) {
                max = val;
            }
        }
        metadata.IsUniform[axisIdx] = false;
        metadata.Extents[axisIdx] = static_cast<float>(max - min);
        metadata.Origin[axisIdx] = static_cast<float>(min);
    }

    // Read the vof data
    vtrFile.seekg(base_offset + std::streamoff(offsets[3]), std::ios::beg);
    volume.resize((metadata.Resolution[0]) * (metadata.Resolution[1]) * (metadata.Resolution[2]));
    auto min = std::numeric_limits<double>::max();
    auto max = std::numeric_limits<double>::min();
    char buffer[sizeof(double)];
    for (size_t i = 0; i < volume.size(); i++) {
        vtrFile.read(buffer, sizeof(double));
        auto val = *reinterpret_cast<double*>(buffer);
        volume[i] = static_cast<float>(val);
        if (val < min) {
            min = val;
        } else if (val > max) {
            max = val;
        }
    }
    this->dataHash++;
    this->fileChanged_ = false;

    metadata.MinValues = new double[1];
    metadata.MinValues[0] = min;
    metadata.MaxValues = new double[1];
    metadata.MaxValues[0] = max;

    // Save the bounding box
    bbox = vislib::math::Cuboid<float>(metadata.Origin[0], metadata.Origin[1], metadata.Origin[2],
        metadata.Origin[0] + metadata.Extents[0], metadata.Origin[1] + metadata.Extents[1], metadata.Origin[2] + metadata.Extents[2]);

    // Report some statistics
    core::utility::log::Log::DefaultLog.WriteInfo("Loaded vtr file %s", vtrFilename.c_str());
    core::utility::log::Log::DefaultLog.WriteInfo("  Resolution: %d x %d x %d", metadata.Resolution[0], metadata.Resolution[1], metadata.Resolution[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("  Extents: %f x %f x %f", metadata.Extents[0], metadata.Extents[1], metadata.Extents[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("  Origin: %f x %f x %f", metadata.Origin[0], metadata.Origin[1], metadata.Origin[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("  Min: %f", min);
    core::utility::log::Log::DefaultLog.WriteInfo("  Max: %f", max);
}

bool trialvolume::vtrFileReader::getDataCallback(core::Call& caller) {
    if (this->fileChanged_) {
        this->fileChanged_ = false;
        this->loadFile();
    }
    if (this->dataHash != 0) {
        auto& dataCall = dynamic_cast<geocalls::VolumetricDataCall&>(caller);
        dataCall.SetDataHash(this->dataHash);
        dataCall.SetFrameID(0);
        dataCall.SetData(this->volume.data());
        dataCall.SetMetadata(&metadata);

        dataCall.AccessBoundingBoxes().SetObjectSpaceBBox(this->bbox);
        auto clipbox = vislib::math::Cuboid<float>(this->bbox);
        clipbox.Grow(.1f);
        dataCall.AccessBoundingBoxes().SetObjectSpaceClipBox(clipbox);
        dataCall.AccessBoundingBoxes().MakeScaledWorld(1.0f);
        return true;
    }
    return false;
}

