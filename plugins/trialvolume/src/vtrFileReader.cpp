#include "VtrFileReader.h"

#include <iostream>
#include <fstream>

#include "mmcore/param/BoolParam.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/IntParam.h"

using namespace megamol;

trialvolume::VtrFileReader::VtrFileReader() : core::Module()
        , get_data_callee_slot_("getdata", "Slot to request data from this data source.")
        , filename_("filename", "The path to the vtr file to load.")
        , data_hash_(0)
        , metadata_()
        , vtr_filename_("")
        , file_changed_(true) {
    // Setup filename input slot
    filename_.SetParameter(new core::param::FilePathParam(""));
    filename_.SetUpdateCallback(&trialvolume::VtrFileReader::filenameChanged);
    MakeSlotAvailable(&filename_);

    // Setup volumetric data output slot
    get_data_callee_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_DATA),
        &trialvolume::VtrFileReader::getDataCallback);
    get_data_callee_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &trialvolume::VtrFileReader::getDataCallback);
    get_data_callee_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &trialvolume::VtrFileReader::getDataCallback);
    get_data_callee_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &trialvolume::VtrFileReader::dummyCallback);
    get_data_callee_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &trialvolume::VtrFileReader::dummyCallback);
    get_data_callee_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &trialvolume::VtrFileReader::dummyCallback);
    MakeSlotAvailable(&get_data_callee_slot_);
}

trialvolume::VtrFileReader::~VtrFileReader() {
    Release();
}

bool trialvolume::VtrFileReader::create(void) {
    return true;
}

void trialvolume::VtrFileReader::release(void) {
}

bool trialvolume::VtrFileReader::dummyCallback(core::Call& caller) {
    return true;
}

// TODO This function is probably *super* unsafe. Replace with something better.
bool trialvolume::VtrFileReader::filenameChanged(core::param::ParamSlot& slot) {
    vtr_filename_ = filename_.Param<core::param::FilePathParam>()->ValueString();
    file_changed_ = true;
    return loadFile();
}

bool trialvolume::VtrFileReader::loadFile() {
    if (vtr_filename_.empty()) {
        core::utility::log::Log::DefaultLog.WriteError("Empty vtr file %s!", vtr_filename_);
        return false;
    }

    // Set general metadata which is not directly given in the vtr file
    metadata_.Components = 1;
    metadata_.GridType = geocalls::GridType_t::RECTILINEAR;
    metadata_.ScalarType = geocalls::ScalarType_t::FLOATING_POINT;
    metadata_.ScalarLength = sizeof(float);

    metadata_.NumberOfFrames = 1;

    // Open file
    std::ifstream vtrFile(vtr_filename_.c_str(), std::ios::binary);
    if (!vtrFile.is_open()) {
        core::utility::log::Log::DefaultLog.WriteError("Could not open vtr file %s!", vtr_filename_.c_str());
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
            core::utility::log::Log::DefaultLog.WriteError("Wrong type given in vtr file %s!", vtr_filename_.c_str());
            return false;
        }
        // Check if expected amount is given
        if (amount != 1) {
            core::utility::log::Log::DefaultLog.WriteError("Wrong amount given in vtr file %s!", vtr_filename_.c_str());
            return false;
        }
        // Check if name is given
        if (name == "x_coordinates") {
            // TODO: Actually use the vof as volumes of the cells not as values on the grid points 
            metadata_.Resolution[0] = size / 8 - 1;
            offsets[0] = offset;
        } else if (name == "y_coordinates") {
            metadata_.Resolution[1] = size / 8 - 1;
            offsets[1] = offset;
        } else if (name == "z_coordinates") {
            metadata_.Resolution[2] = size / 8 - 1;
            offsets[2] = offset;
        } else if (name == "vof-function[-]") {
            offsets[3] = offset;
        } else {
            core::utility::log::Log::DefaultLog.WriteError("Unknown name %s given in vtr file %s!", name, vtr_filename_.c_str());
            return false;
        }
    }

    auto const base_offset = vtrFile.tellg();

    // Read all the axis data
    for (size_t axisIdx = 0; axisIdx < 3; axisIdx++) {

        vtrFile.seekg(base_offset + std::streamoff(offsets[axisIdx]), std::ios::beg);
        metadata_.SliceDists[axisIdx] = new float[metadata_.Resolution[axisIdx] - 1];
        char buffer[sizeof(double)];
        vtrFile.read(buffer, sizeof(double));
        auto prev_val = *reinterpret_cast<double*>(buffer);
        auto min = prev_val;
        auto max = prev_val;
        for (size_t i = 0; i < metadata_.Resolution[axisIdx] - 1; i++) {
            vtrFile.read(buffer, sizeof(double));
            auto const val = *reinterpret_cast<double*>(buffer);
            metadata_.SliceDists[axisIdx][i] = static_cast<float>(val - prev_val);
            prev_val = val;
            if (val < min) {
                min = val;
            } else if (val > max) {
                max = val;
            }
        }
        metadata_.IsUniform[axisIdx] = false;
        metadata_.Extents[axisIdx] = static_cast<float>(max - min);
        metadata_.Origin[axisIdx] = static_cast<float>(min);
    }

    // Read the vof data
    vtrFile.seekg(base_offset + std::streamoff(offsets[3]), std::ios::beg);
    volume_.resize((metadata_.Resolution[0]) * (metadata_.Resolution[1]) * (metadata_.Resolution[2]));
    auto min = std::numeric_limits<double>::max();
    auto max = std::numeric_limits<double>::min();
    char buffer[sizeof(double)];
    for (size_t i = 0; i < volume_.size(); i++) {
        vtrFile.read(buffer, sizeof(double));
        auto val = *reinterpret_cast<double*>(buffer);
        volume_[i] = static_cast<float>(val);
        if (val < min) {
            min = val;
        } else if (val > max) {
            max = val;
        }
    }
    data_hash_++;
    file_changed_ = false;

    metadata_.MinValues = new double[1];
    metadata_.MinValues[0] = min;
    metadata_.MaxValues = new double[1];
    metadata_.MaxValues[0] = max;

    // Save the bounding box
    bbox_ = vislib::math::Cuboid<float>(metadata_.Origin[0], metadata_.Origin[1], metadata_.Origin[2],
        metadata_.Origin[0] + metadata_.Extents[0], metadata_.Origin[1] + metadata_.Extents[1], metadata_.Origin[2] + metadata_.Extents[2]);

    // Report some statistics
    core::utility::log::Log::DefaultLog.WriteInfo("Loaded vtr file %s", vtr_filename_.c_str());
    core::utility::log::Log::DefaultLog.WriteInfo("  Resolution: %d x %d x %d", metadata_.Resolution[0], metadata_.Resolution[1], metadata_.Resolution[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("  Extents: %f x %f x %f", metadata_.Extents[0], metadata_.Extents[1], metadata_.Extents[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("  Origin: %f x %f x %f", metadata_.Origin[0], metadata_.Origin[1], metadata_.Origin[2]);
    core::utility::log::Log::DefaultLog.WriteInfo("  Min: %f", min);
    core::utility::log::Log::DefaultLog.WriteInfo("  Max: %f", max);
}

bool trialvolume::VtrFileReader::getDataCallback(core::Call& caller) {
    if (file_changed_) {
        file_changed_ = false;
        loadFile();
    }
    if (data_hash_ != 0) {
        auto& dataCall = dynamic_cast<geocalls::VolumetricDataCall&>(caller);
        dataCall.SetDataHash(data_hash_);
        dataCall.SetFrameID(0);
        dataCall.SetData(volume_.data());
        dataCall.SetMetadata(&metadata_);

        dataCall.AccessBoundingBoxes().SetObjectSpaceBBox(bbox_);
        auto clipbox = vislib::math::Cuboid<float>(bbox_);
        clipbox.Grow(.1f);
        dataCall.AccessBoundingBoxes().SetObjectSpaceClipBox(clipbox);
        dataCall.AccessBoundingBoxes().MakeScaledWorld(1.0f);
        return true;
    }
    return false;
}

