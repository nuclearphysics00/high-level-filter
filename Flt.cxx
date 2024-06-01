#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <string>
#include "UnpackTdc.h"
#include "SubTimeFrameHeader.h"

class Flt
{
public:
    Flt();
    virtual ~Flt();
    void InitParam();
    std::vector<std::vector<int>> clusterNumbers(const std::vector<int>& numbers);

    struct GeofieldClusterInfo {
        int clusterCount; // cluster size 
        std::vector<int> clusterSizes; // number of elements in cluster
    };

    std::map<int, GeofieldClusterInfo> analyzeGeofieldClusters(const std::map<int, std::vector<int>>& GeoIDs);
    std::map<int, GeofieldClusterInfo> geofieldClusterMap;    
    
    typedef struct {
        int catid = -1;
        int id = -1;
        int sh = -1;
        int fp = -1;
        int det = -1;
        uint64_t geo = -1;
        int ch = -1;
    } Wire_map;

    //int geoToIndex(uint64_t geo);
    //int findWirenumber(Flt::Wire_map *record, int size, uint64_t geo, int ch, int *foundid, int *foundGeo, int &Geofield);
    //bool allKeysHaveAtLeastOneCluster(const std::map<int, GeofieldClusterInfo>& geofieldClusterMap);
    //void Charge(int charge, int fem, int &ChargeSum, int &u2ChargeSum, int &x1ChargeSum, int &u1ChargeSum, int &x2ChargeSum, int &nu1, int &nu2);
    //bool Chargelogic(int ChargeSum, int u2ChargeSum, int nu1, int nu2);
    //bool TOFlogic(int tdc1, int tdc2);
    std::vector<std::unique_ptr<int>> CalculateAveragePairs(const std::vector<std::unique_ptr<int>>& r_hits, const std::vector<std::unique_ptr<int>>& l_hits);
    void CalculateAndPrintTOF(const std::vector<std::unique_ptr<int>>& utof_averages, const std::vector<std::unique_ptr<int>>& t1_averages);
    void PrintList(const std::vector<std::unique_ptr<int>>& list, const std::string& listName);
    bool CheckAllTOFConditions(int t1tof_ave, int utof_ave, const std::vector<std::unique_ptr<int>>& t1tof_r_hits, const std::vector<std::unique_ptr<int>>& t1tof_l_hits,
                               const std::vector<std::unique_ptr<int>>& utof_r_hits, const std::vector<std::unique_ptr<int>>& utof_l_hits, int min, int max);

    std::map<int, std::vector<int>> GeoIDs; // ワイヤーIDと面ID(上流側から数えた時の面番号)を保存してるマップ
};

Flt::Flt() {
    return;
}

Flt::~Flt() {

    return;
}

// TOF 特定の検出機のマルチヒットの中から全てのRとLのコンビネーションを探し、平均値を計算。2
std::vector<std::unique_ptr<int>> Flt::CalculateAveragePairs(const std::vector<std::unique_ptr<int>>& r_hits, const std::vector<std::unique_ptr<int>>& l_hits) {
    std::vector<std::unique_ptr<int>> averages;
    for (const auto& r : r_hits) {
        for (const auto& l : l_hits) {
            averages.push_back(std::make_unique<int>((*r + *l) / 2.0));
        }
    }
    return averages;
}

// TOF ２つの検出器間のTOFを算出、検出器1 検出器2 の平均値全てのコンビネーションを計算 3
void Flt::CalculateAndPrintTOF(const std::vector<std::unique_ptr<int>>& utof_averages, const std::vector<std::unique_ptr<int>>& t1_averages) {
    for (const auto& utof_ave : utof_averages) {
        for (const auto& t1_ave : t1_averages) {
            int tof = *utof_ave - *t1_ave;
            //std::cout << "UTOF Average: " << *utof_ave << ", T1 Average: " << *t1_ave << ", TOF: " << tof << std::endl;
        }
    }
}

// TOFデバック用 1
void Flt::PrintList(const std::vector<std::unique_ptr<int>>& list, const std::string& listName) {
    std::cout << listName << " [";
    for (const auto& value : list) {
        std::cout << *value << " ";
    }
    std::cout << "]" << std::endl;
}

// TOFフラグ
bool Flt::CheckAllTOFConditions(int t1tof_ave, int utof_ave, const std::vector<std::unique_ptr<int>>& t1tof_r_hits, const std::vector<std::unique_ptr<int>>& t1tof_l_hits,
                                const std::vector<std::unique_ptr<int>>& utof_r_hits, const std::vector<std::unique_ptr<int>>& utof_l_hits, int min, int max) {
    int avg_tof = utof_ave - t1tof_ave;
    return (avg_tof > min && avg_tof < max);
}

// マルチプリシティGeoと面番号を修正する必要がある。
/* int Flt::findWirenumber(Flt::Wire_map *record, int size, uint64_t geo, int ch, int *foundid, int *foundGeo, int &Geofield){    
    for(int i=0; i<size; i++){
        if( (geo == record[i].geo) && (ch == record[i].ch) ){
            *foundid = record[i].id;
            *foundGeo = record[i].geo;
            if((geo == 0xc0a802a1) || (geo == 0xc0a802a2) ){
                Geofield =1 ;//plane
            }
            else if((geo == 0xc0a802a3) || (geo == 0xc0a802a4) ){
                Geofield = 2 ;//plane2
            }
            else if((geo == 0xc0a802a5) || (geo == 0xc0a802a6) ){
                Geofield = 3 ;//plane3
            }
            else if((geo == 0xc0a802a7) || ((geo == 0xc0a802aa))){
                Geofield = 4 ;//plane4
            }
            return 1;
        } 
    }
    return 0;
} */

/* int Flt::geoToIndex(uint64_t geo) {
    if (geo >= 0xc0a802a1 && geo <= 0xc0a802a7) {
        return geo - 0xc0a802a1;
    } else if (geo == 0xc0a802aa) {
        return 7;
    }
    return -1; 
} */

/* std::vector<std::vector<int>> Flt::clusterNumbers(const std::vector<int>& numbers) {
    std::vector<std::vector<int>> clusters;
    std::vector<int> currentCluster;

    for (int number : numbers) {
        if (currentCluster.empty() || number == currentCluster.back() + 1) {
            currentCluster.push_back(number);
        } else {
            if (currentCluster.size() >= 3) {
                clusters.push_back(currentCluster);
            }
            currentCluster.clear();
            currentCluster.push_back(number);
        }
    }

    if (currentCluster.size() >= 3) {
        clusters.push_back(currentCluster);
    }

    return clusters;
} */

// マルチプリシティメイン
/* std::map<int, Flt::GeofieldClusterInfo> Flt::analyzeGeofieldClusters(const std::map<int, std::vector<int>>& GeoIDs) {
    std::map<int, GeofieldClusterInfo> geofieldClusterMap;

    for (int geofieldID = 1; geofieldID <= 4; ++geofieldID) {
        geofieldClusterMap[geofieldID] = GeofieldClusterInfo{0, std::vector<int>()};
    }

    for (const auto& pair : GeoIDs) {
        int geofieldID = pair.first;
        const std::vector<int>& wireIDs = pair.second;
        std::vector<std::vector<int>> clusters = clusterNumbers(wireIDs);

        GeofieldClusterInfo& info = geofieldClusterMap[geofieldID];
        info.clusterCount = clusters.size();
        for (const auto& cluster : clusters) {
            info.clusterSizes.push_back(cluster.size());
        }
    }

    return geofieldClusterMap;
} */

// マルチプリシティ flag
/* bool Flt::allKeysHaveAtLeastOneCluster(const std::map<int, Flt::GeofieldClusterInfo>& geofieldClusterMap) {
    for (const auto& pair : geofieldClusterMap) {
        if (pair.second.clusterCount < 1) {
            return false;
        }
    }
    return true;
} */

// 
/* void Flt::Charge(int charge,int fem, int &ChargeSum, int &u2ChargeSum, int &x1ChargeSum, int &u1ChargeSum, int &x2ChargeSum, int &nu1, int &nu2){
    if((fem == 0xc0a802a1) || (fem == 0xc0a802a2) || (fem == 0xc0a802a3) || (fem == 0xc0a802a4) || (fem == 0xc0a802a5) || (fem == 0xc0a802a6) || (fem == 0xc0a802a7) || (fem == 0xc0a802aa)){
        ChargeSum += charge;
        if((fem == 0xc0a802a7) || ((fem == 0xc0a802aa))){
            u2ChargeSum += charge;
            nu2 ++;
        }
        else if((fem == 0xc0a802a3) || (fem == 0xc0a802a4) ){
            u1ChargeSum += charge;
            nu1 ++;
        }
    }
} */

// 電荷 flag
/* bool Flt::Chargelogic(int ChargeSum, int u2ChargeSum, int nu1, int nu2){
    if((((double)(ChargeSum)/nu1) > 35)  && ((double)(u2ChargeSum)/nu2 > 45)){
        return 1;
    }
    return 0;
} */


/* v1シングルヒットの時のTOFコード
bool Flt::TOFlogic(int tdc1,int tdc2){
    int tof = tdc1 - tdc2;
    std::cout << "tof: " << tof << std::endl;
    std::cout << "abs(tof): " << abs(tof) << std::endl; 
    if(tof> -130000 && tof < -125000){
        return 1;
    }
    return 0;
} */