// System includes
#include <vector>
#include <iostream>
#include <array>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <set>
#include <memory>
#include <string>
#include <sstream>
#include <chrono>

// Project includes


template<unsigned int TDim>
class DataPoint
{
public:
    ///@name Type definitions
    ///@{

    using Pointer = std::shared_ptr<DataPoint>;

    ///@}
    ///@name Life cycle
    ///@{

    DataPoint(
        const std::size_t Id,
        const std::array<double, TDim>& rPosition)
        : mId(Id),
          mPosition(rPosition)
    {
    }

    ///@}
    ///@Public operations
    ///@{

    std::size_t GetId() const { return mId; }

    std::array<double, TDim> Coordinates() const { return mPosition; }

    std::string Info() const
    {
        std::stringstream info;
        info << mId << "[" << mPosition[0];
        for (std::size_t i = 1; i < TDim; ++i) {
            info << "," << mPosition[i];
        }
        info << "]";
        return info.str();
    }

    ///@}
    ///@name Static Position
    ///@{

    static std::array<double, TDim> GetPosition(const DataPoint<TDim>::Pointer& pDataPoint)
    {
        return pDataPoint->Coordinates();
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::size_t mId;

    const std::array<double, TDim> mPosition;

    ///@}
};

template<unsigned int TDim>
std::ostream& operator<<(
    std::ostream& rOStream,
    const DataPoint<TDim>& rThis)
{
    rOStream << rThis.Info();
    return rOStream;
}

template<unsigned int TDim>
std::ostream& operator<<(
    std::ostream& rOStream,
    const typename DataPoint<TDim>::Pointer& pThis)
{
    rOStream << pThis->Info();
    return rOStream;
}


#include "../includes/KDTree.h"

int main()
{
    constexpr unsigned int Dim = 2;
    // create the data points
    std::size_t N = 6e+5;
    std::vector<DataPoint<Dim>::Pointer> points;

    // test data
    const double distance = 1;
    const std::array<double, Dim> test_location = {3.5, 3.3};

    for (std::size_t i = 0; i < N; ++i) {
        std::array<double, Dim> position;
        for (std::size_t j = 0; j < Dim; ++j) {
            const auto rand = std::rand();
            position[j] = 0.0 + 10.0 * rand / RAND_MAX;
        }
        points.push_back(std::make_shared<DataPoint<Dim>>(i + 1, position));
    }

    using k_d_tree_type = KDTree<Dim, DataPoint<Dim>::Pointer,  DataPoint<Dim>>;
    k_d_tree_type kd_tree(points.begin(), points.end());
    // std::cout<< "---------------- KDTree: ----------------";
    // std::cout << kd_tree.Info() << "\n";

    // now find the neighbours using kd tree
    std::vector<k_d_tree_type::TreeNodeType::Pointer> p_tree_nodes;
    p_tree_nodes.resize(2);
    std::cout <<"\nFinding neighbours using KDTree...";
    const auto kd_start = std::chrono::high_resolution_clock::now();
    const auto number_of_neighbours = kd_tree.FindNearestNeighbours(p_tree_nodes, test_location, distance);
    const auto kd_end = std::chrono::high_resolution_clock::now();
    std::cout <<"\nFound neighbours using KDTree within " << std::chrono::duration_cast<std::chrono::microseconds>(kd_end - kd_start).count() << " us.";
    std::sort(p_tree_nodes.begin(), p_tree_nodes.begin() + number_of_neighbours, [](auto pNode1, auto pNode2) { return pNode1->GetData()->GetId() < pNode2->GetData()->GetId(); });

    // now find the neighbours using brute force
    std::cout <<"\nFinding neighbours using brute force...";
    std::vector<DataPoint<Dim>::Pointer> brute_force_search_results;
    const double squared_distance = distance * distance;
    const auto bf_start = std::chrono::high_resolution_clock::now();
    for (const auto& p_point : points) {
        const auto& pos = p_point->Coordinates();
        double current_distance_square = 0.0;
        for (std::size_t i = 0; i < Dim; ++i) {
            current_distance_square += std::pow(pos[i] - test_location[i], 2);
        }
        if (current_distance_square <= squared_distance) {
            brute_force_search_results.push_back(p_point);
        }
    }
    const auto bf_end = std::chrono::high_resolution_clock::now();
    std::cout <<"\nFound neighbours using brute force within " << std::chrono::duration_cast<std::chrono::microseconds>(bf_end - bf_start).count() << " us.";

    // now the check
    bool is_correct = number_of_neighbours == brute_force_search_results.size();
    std::sort(brute_force_search_results.begin(), brute_force_search_results.end(), [](const auto& rA, const auto& rB) { return rA->GetId() < rB->GetId(); });
    if (is_correct) {
        for (std::size_t i = 0; i < number_of_neighbours; ++i) {
            if (p_tree_nodes[i]->GetData()->GetId() != brute_force_search_results[i]->GetId()) {
                is_correct = false;
            }
        }
    }

    // std::cout <<"\nBrute force neighbours:";
    // for (const auto& p_point : brute_force_search_results) {
    //     std::cout << "\n\t" << *p_point;
    // }

    // std::cout << "\nKDTree neighbours:";
    // for (std::size_t i = 0; i < number_of_neighbours; ++i) {
    //     std::cout << "\n\t" << *(p_tree_nodes[i]->GetData());
    // }

    std::cout <<"\n ---------------- Summary ----------------";
    std::cout <<"\nTest position          = [" << test_location[0];
    for (std::size_t i = 1; i < Dim; ++i) {
        std::cout << ", " << test_location[i];
    }
    std::cout << "]";
    std::cout <<"\nKDTree neighbours      = " << number_of_neighbours;
    std::cout <<"\nBrute force neighbours = " << brute_force_search_results.size();
    std::cout <<"\nIS CORRECT             = " << (is_correct ? "YES" : "NO");
    std::cout <<"\nN                      = " << N << std::endl;

    return 0;
}