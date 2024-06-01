#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
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

bool FilterTimeFrameSliceBySomething::ProcessSlice(TTF& tf)
{

    //std::cout << "関数内部" << std::endl;

    std::vector<std::unique_ptr<int>> t1tof_r_hits;
    std::vector<std::unique_ptr<int>> t1tof_l_hits;
    std::vector<std::unique_ptr<int>> utof_r_hits;
    std::vector<std::unique_ptr<int>> utof_l_hits;

    bool tofCondition = false;

    auto tfHeader = tf.GetHeader();
    //std::cout << "TimeFrameヘッダーサイズ: " << sizeof(*tfHeader) << " bytes" << std::endl;

    for (auto& SubTimeFrame : tf) {
        //std::cout << "SubTimeFrame内" << std::endl;
        auto header = SubTimeFrame->GetHeader();
        auto& hbf = SubTimeFrame->at(0); // SubTimeFrameの上から１つ目のHB

        //std::cout << "SubTimeFrameヘッダーサイズ: " << sizeof(*header) << " bytes" << std::endl;
        uint64_t nData = hbf->GetNumData();

        for (int i = 0; i < nData; ++i) {
            if (header->femType == SubTimeFrame::TDC64H) {
                TDC64H::tdc64 tdc;
                TDC64H::Unpack(hbf->UncheckedAt(i), &tdc);

                if (header->femId == 0xc0a802aa) {
                    //std::cout << "64H_utofv1: " << tdc.tdc << std::endl;
                    //std::cout << "64H_t1v1: " << tdc.tdc << std::endl;
                }
            } else if (header->femType == SubTimeFrame::TDC64L) {
                TDC64L::tdc64 tdc;
                TDC64L::Unpack(hbf->UncheckedAt(i), &tdc);
                if (header->femId == 0xc0802aa9 && (tdc.ch == 10 || tdc.ch == 12)) {
                    //std::cout << "64L_utof: " << tdc.tdc << std::endl;
                    //std::cout << "64L_t1: " << tdc.tdc << std::endl;
                }
            } else if (header->femType == SubTimeFrame::TDC64H_V3) {
                TDC64H_V3::tdc64 tdc;
                TDC64H_V3::Unpack(hbf->UncheckedAt(i), &tdc);
                if (header->femId == 0xc0a802aa) {
                    if (tdc.ch == 12) {
                        // TOF
                        t1tof_r_hits.push_back(std::make_unique<int>(tdc.tdc));
                    } else if (tdc.ch == 10) {
                        // TOF
                        t1tof_l_hits.push_back(std::make_unique<int>(tdc.tdc));
                    }
                } else if (header->femId == 0xc0a802a9) {
                    if (tdc.ch == 10) {
                        // TOF
                        utof_r_hits.push_back(std::make_unique<int>(tdc.tdc));
                    } else if (tdc.ch == 8) {
                        // TOF
                        utof_l_hits.push_back(std::make_unique<int>(tdc.tdc));
                    }
                }
            } else if (header->femType == SubTimeFrame::TDC64L_V3) {
                TDC64L_V3::tdc64 tdc;
                TDC64L_V3::Unpack(hbf->UncheckedAt(i), &tdc);
                // Add processing for TDC64L_V3 if needed
            }

            /* マルチプリシティ前段処理
            if (wireMap.find(header->femId) != wireMap.end()) {
                Flt::Wire_map& wireMapEntry = wireMap[header->femId];
                GeoIDs[header->femId].push_back(wireMapEntry.geo);
                std::cout << "Wire ID: " << wireMapEntry.id << std::endl;
            } */
        }
    }

    // デバッグ用
/*     flt->PrintList(t1tof_r_hits, "t1tof_r_hits");
    flt->PrintList(t1tof_l_hits, "t1tof_l_hits");
    flt->PrintList(utof_r_hits, "utof_r_hits");
    flt->PrintList(utof_l_hits, "utof_l_hits"); */

    // 面ごとの各ペアの平均値の計算と表示
    if ((!t1tof_r_hits.empty() || !t1tof_l_hits.empty()) && (!utof_r_hits.empty() || !utof_l_hits.empty())) {
        auto t1_averages = flt->CalculateAveragePairs(t1tof_r_hits, t1tof_l_hits);
        auto utof_averages = flt->CalculateAveragePairs(utof_r_hits, utof_l_hits);

        flt->CalculateAndPrintTOF(utof_averages, t1_averages);

        for (const auto& t1_avg : t1_averages) {
            for (const auto& utof_avg : utof_averages) {
                //最小値 最大値の順番
                if (flt->CheckAllTOFConditions(*t1_avg, *utof_avg, t1tof_r_hits, t1tof_l_hits, utof_r_hits, utof_l_hits, -130000, -125000)) {
                    tofCondition = true;
                    //std::cout << "合格tof: " << *utof_avg - *t1_avg << std::endl;
                    //std::cout << "min: -130000" << " " << "max: -125000" << std::endl;
                    break;
                }
           }
            if (tofCondition) break;
        }
    }

    if (tofCondition) {
        // TOF フィルターデバッグ用 TDC値分布を見る用
        /* std::ofstream outFile("tof_output.txt", std::ios::app);
        double tof = *utof_ave - *t1tof_ave;
        outFile << *utof_ave << " " << *t1tof_ave << " " << tof << "\n";
        outFile.close(); */

        //std::cout << "関数内tofbool値: true" << std::endl;
        return true;
    }

    return false;
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
