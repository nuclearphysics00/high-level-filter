#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <iomanip>
#include <cmath>
#include "FilterTimeFrameSliceBySomething.h"
#include "FilterTimeFrameSliceABC.icxx"
#include "fairmq/runDevice.h"
#include "utility/MessageUtil.h"
#include "UnpackTdc.h"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"

#define DEBUG 0

using nestdaq::FilterTimeFrameSliceBySomething;
namespace bpo = boost::program_options;

FilterTimeFrameSliceBySomething::FilterTimeFrameSliceBySomething()
{
}

struct GroupInfo {
    uint64_t ip;
    int ch_r;
    int ch_l;
};

// LTOF {detectorID, {IP, r, l}}
std::map<int, GroupInfo> groupInfoMap = {
    {1, {0xc0a802aa, 16, 22}}, // LTOF
    {2, {0xc0a802aa, 19, 25}}, // LTOF
    {3, {0xc0a802aa, 17, 23}}, // LTOF
    {4, {0xc0a802aa, 20, 27}}, // LTOF
    {5, {0xc0a802aa, 18, 24}}, // LTOF
    {6, {0xc0a802aa, 21, 28}}, // LTOF
};

bool FilterTimeFrameSliceBySomething::ProcessSlice(TTF& tf)
{
    std::cout << "関数内部" << std::endl;

    // initialize
    std::vector<std::shared_ptr<int>> t1tof_r_hits;
    std::vector<std::shared_ptr<int>> t1tof_l_hits;
    std::vector<std::shared_ptr<int>> utof_r_hits;
    std::vector<std::shared_ptr<int>> utof_l_hits;
    std::vector<std::vector<std::shared_ptr<int>>> ltof_r_hits(6); // 6つのTOF検出器群の右ヒット
    std::vector<std::vector<std::shared_ptr<int>>> ltof_l_hits(6); // 6つのTOF検出器群の左ヒット

    bool t1tofCondition = false;
    bool ltofCondition = false;

    auto tfHeader = tf.GetHeader();

    for (auto& SubTimeFrame : tf) {
        auto header = SubTimeFrame->GetHeader();
        auto& hbf = SubTimeFrame->at(0); 
        uint64_t nData = hbf->GetNumData();

        for (int i = 0; i < nData; ++i) {
            if (header->femType == SubTimeFrame::TDC64H) {
                TDC64H::tdc64 tdc;
                TDC64H::Unpack(hbf->UncheckedAt(i), &tdc);

                if (header->femId == 0xc0a802aa) {
                    //If you need the process, write it down.
                }
            } else if (header->femType == SubTimeFrame::TDC64L) {
                TDC64L::tdc64 tdc;
                TDC64L::Unpack(hbf->UncheckedAt(i), &tdc);
                if (header->femId == 0xc0802aa9 && (tdc.ch == 10 || tdc.ch == 12)) {
                    //If you need the process, write it down.
                }
            } else if (header->femType == SubTimeFrame::TDC64H_V3) {
                TDC64H_V3::tdc64 tdc;
                TDC64H_V3::Unpack(hbf->UncheckedAt(i), &tdc);
                if (header->femId == 0xc0a802aa) {
                    if (tdc.ch == 12) {
                        // TOF process
                        t1tof_r_hits.push_back(std::make_shared<int>(tdc.tdc));
                    } else if (tdc.ch == 10) {
                        // TOF process
                        t1tof_l_hits.push_back(std::make_shared<int>(tdc.tdc));
                    }
                    // LTOF process
                    for (const auto& groupInfo : groupInfoMap) {
                        if (header->femId == groupInfo.second.ip) {
                            if (tdc.ch == groupInfo.second.ch_r) {
                                // TOF
                                ltof_r_hits[groupInfo.first - 1].push_back(std::make_shared<int>(tdc.tdc));
                            } else if (tdc.ch == groupInfo.second.ch_l) {
                                // TOF
                                ltof_l_hits[groupInfo.first - 1].push_back(std::make_shared<int>(tdc.tdc));
                            }
                        }
                    }
                } else if (header->femId == 0xc0a802a9) {
                    if (tdc.ch == 10) {
                        // TOF
                        utof_r_hits.push_back(std::make_shared<int>(tdc.tdc));
                    } else if (tdc.ch == 8) {
                        // TOF
                        utof_l_hits.push_back(std::make_shared<int>(tdc.tdc));
                    }
                }
            } else if (header->femType == SubTimeFrame::TDC64L_V3) {
                TDC64L::tdc64 tdc;
                TDC64L::Unpack(hbf->UncheckedAt(i), &tdc);
                // Add processing for TDC64L_V3 if needed
            }
        }
    }

    // for debug
    flt->PrintList(t1tof_r_hits, "t1tof_r_hits");
    flt->PrintList(t1tof_l_hits, "t1tof_l_hits");
    flt->PrintList(utof_r_hits, "utof_r_hits");
    flt->PrintList(utof_l_hits, "utof_l_hits");
    for (size_t i = 0; i < 6; ++i) {
        flt->PrintList(ltof_r_hits[i], "ltof_r_hits[" + std::to_string(i + 1) + "]");
        flt->PrintList(ltof_l_hits[i], "ltof_l_hits[" + std::to_string(i + 1) + "]");
    }

    // calculate averages R and L
    auto t1_averages = flt->CalculateAveragePairs(t1tof_r_hits, t1tof_l_hits);
    auto utof_averages = flt->CalculateAveragePairs(utof_r_hits, utof_l_hits);
    std::vector<std::vector<std::shared_ptr<int>>> ltof_averages(6);
    for (size_t i = 0; i < 6; ++i) {
        ltof_averages[i] = flt->CalculateAveragePairs(ltof_r_hits[i], ltof_l_hits[i]);
    }

    // calculate averages of TOF
    flt->CalculateAndPrintTOF(utof_averages, t1_averages);
    for (const auto& t1_avg : t1_averages) {
        for (const auto& utof_avg : utof_averages) {
            double tof_difference = *utof_avg - *t1_avg;
            std::cout << "utof_avg: " << *utof_avg << ", t1_avg: " << *t1_avg << ", TOF Difference: " << tof_difference << std::endl;

            if (flt->CheckAllTOFConditions(*t1_avg, *utof_avg, t1tof_r_hits, t1tof_l_hits, utof_r_hits, utof_l_hits, -130000, -125000)) {
                t1tofCondition = true;
                std::cout << "合格tof: " << *utof_avg - *t1_avg << std::endl;
                std::cout << "min: -130000" << " " << "max: -125000" << std::endl;
            }
        }
    }

    // Dump the result of subtracting the average of utof and ltof as a list
    for (const auto& utof_avg : utof_averages) {
        for (size_t i = 0; i < 6; ++i) {
            for (const auto& ltof_avg : ltof_averages[i]) {
                double tof_difference = *utof_avg - *ltof_avg;
                std::cout << "utof_avg: " << *utof_avg << ", ltof_avg[" << i + 1 << "]: " << *ltof_avg << ", TOF Difference: " << tof_difference << std::endl;

                if (flt->CheckAllTOFConditions(*ltof_avg, *utof_avg, ltof_r_hits[i], ltof_l_hits[i], utof_r_hits, utof_l_hits, -130000, -125000)) {
                    ltofCondition = true;
                }
            }
        }
    }

    return t1tofCondition && ltofCondition;
}

void addCustomOptions(bpo::options_description& options)
{
    using opt = FilterTimeFrameSliceBySomething::OptionKey;

    options.add_options()
        (opt::InputChannelName.data(),
         bpo::value<std::string>()->default_value("in"),
         "Name of the input channel")
        (opt::OutputChannelName.data(),
         bpo::value<std::string>()->default_value("out"),
         "Name of the output channel")
        (opt::DQMChannelName.data(),
         bpo::value<std::string>()->default_value("dqm"),
         "Name of the data quality monitoring channel")
        (opt::PollTimeout.data(),
         bpo::value<std::string>()->default_value("1"),
         "Timeout of polling (in msec)")
        (opt::SplitMethod.data(),
         bpo::value<std::string>()->default_value("1"),
         "STF split method");
}

std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
    return std::make_unique<FilterTimeFrameSliceBySomething>();
}
