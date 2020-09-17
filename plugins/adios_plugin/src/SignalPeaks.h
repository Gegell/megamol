#pragma once

#include <numeric>
#include <memory>

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

#include "mmstd_datatools/table/TableDataCall.h"

namespace megamol {
namespace adios {

class SignalPeaks : public core::Module {
public:
    enum class PeakSelector : int { EQUIDISTANT, MAXVARIANCE, NTHHIGHEST };

    /** Return module class name */
    static const char* ClassName(void) { return "SignalPeaks"; }

    /** Return module class description */
    static const char* Description(void) { return "Extracts local peaks from a signal"; }

    /** Module is always available */
    static bool IsAvailable(void) { return true; }

    /** Ctor */
    SignalPeaks(void);

    /** Dtor */
    virtual ~SignalPeaks(void);

protected:
    /** Lazy initialization of the module */
    bool create(void) override;

    /** Resource release */
    void release(void) override;

private:
    bool getDataCallback(core::Call& c);
    bool getHeaderCallback(core::Call& c);

    bool isDirty() {
        return num_peaks_slot_.IsDirty() || peak_selector_slot_.IsDirty() || column_grouping_factor_slot_.IsDirty();
    }

    void resetDirty() {
        num_peaks_slot_.ResetDirty();
        peak_selector_slot_.ResetDirty();
        column_grouping_factor_slot_.ResetDirty();
    }

    std::vector<size_t> equidistant(
        float const* data, size_t num_cols, size_t num_rows, size_t num_peaks) {
        return std::vector<size_t>();
    }
    std::vector<size_t> maxvariance(float const* data, size_t num_cols, size_t num_rows, size_t num_peaks) {
        std::vector<std::pair<size_t, float>> row_variance;
        row_variance.reserve(num_rows);

        for (size_t row = 0; row < num_rows; ++row) {
            auto const minmax = std::minmax_element(data + (row * num_cols), data + ((row + 1) * num_cols));
            row_variance.push_back(std::make_pair(row, *minmax.second - *minmax.first));
        }

        std::sort(row_variance.begin(), row_variance.end(),
            [](auto const& lhs, auto const& rhs) { return lhs.second > rhs.second; });

        std::vector<size_t> ret(num_peaks);
        for (size_t i = 0; i < num_peaks; ++i) {
            ret[i] = row_variance[i].first;
        }

        return ret;
    }
    std::vector<size_t> nthhighest(float const* data, size_t num_cols, size_t num_rows, size_t num_peaks) {
        std::vector<std::pair<size_t, float>> row_max;
        row_max.reserve(num_rows);

        for (size_t row = 0; row < num_rows; ++row) {
            auto const max =
                std::max_element(data + (row * num_cols), data + ((row + 1) * num_cols));
            row_max.push_back(std::make_pair(row, *max));
        }

        std::sort(row_max.begin(), row_max.end(),
            [](auto const& lhs, auto const& rhs) { return lhs.second > rhs.second; });

        std::vector<size_t> ret(num_peaks);
        for (size_t i = 0; i < num_peaks; ++i) {
            ret[i] = row_max[i].first;
        }

        return ret;
    }

    std::vector<float> extract_peaks(
        float const* data, size_t num_cols, size_t num_rows, std::vector<size_t> indices) {
        auto const num_peaks = indices.size();
        std::vector<float> ret(num_cols * num_peaks);

        for (size_t i = 0; i < num_peaks; ++i) {
            auto const idx = indices[i];
            std::copy(
                data + (idx * num_cols), data + ((idx + 1) * num_cols), ret.begin() + (i * num_cols));
        }

        return ret;
    }

    std::vector<std::pair<float, float>> reevaluate_colums(
        std::vector<float> const& data, size_t num_cols, size_t num_rows) {
        std::vector<std::pair<float, float>> ret;
        ret.reserve(num_cols);

        for (size_t col = 0; col < num_cols; ++col) {
            float min_val = std::numeric_limits<float>::max();
            float max_val = std::numeric_limits<float>::lowest();

            for (size_t row = 0; row < num_rows; ++row) {
                auto const val = data[col + row * num_cols];
                if (val < min_val) min_val = val;
                if (val > max_val) max_val = val;
            }

            ret.push_back(std::make_pair(min_val, max_val));
        }

        return ret;
    }

    void populate_ranges(std::vector<stdplugin::datatools::table::TableDataCall::ColumnInfo>& infos,
        std::vector<std::pair<float, float>> const& ranges) {
        for (size_t col = 0; col < infos.size(); ++col) {
            infos[col].SetMinimumValue(ranges[col].first);
            infos[col].SetMaximumValue(ranges[col].second);
        }
    }

    core::CalleeSlot data_out_slot_;
    core::CallerSlot data_in_slot_;

    core::param::ParamSlot num_peaks_slot_;
    core::param::ParamSlot peak_selector_slot_;
    core::param::ParamSlot column_grouping_factor_slot_;

    std::vector<float> data_;
    std::vector<stdplugin::datatools::table::TableDataCall::ColumnInfo> infos_;
    size_t out_num_columns_;
    size_t out_num_rows_;

    size_t data_hash_;
}; // class SignalPeaks

} // namespace adios
} // namespace megamol
