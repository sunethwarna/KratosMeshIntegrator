// System includes
#include <array>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <type_traits>
#include <chrono>

// Project includes

template<unsigned int TDim, class TDataType, class TPositionGetterFunctor>
class KDTree;

template<unsigned int TDim, class TDataType, class TPositionGetterFunctor>
class TreeNode
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using TreeNodeType = TreeNode<TDim, TDataType, TPositionGetterFunctor>;

    using Pointer = TreeNodeType const*;

    using PositionType = std::result_of_t<decltype(&TPositionGetterFunctor::GetPosition)(const TDataType&)>;

    static IndexType Counter;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TIteratorType>
    TreeNode(
        TIteratorType Begin,
        TIteratorType End,
        const IndexType Dimension)
        : mDimension(Dimension % TDim)
    {
        const IndexType number_of_items = std::distance(Begin, End);

        if (number_of_items == 1) {
            mData = *Begin;
        } else if (number_of_items == 2) {
            mData = *Begin;

            const auto& r_position = TPositionGetterFunctor::GetPosition(*(Begin + 1));
            if (r_position[mDimension] < GetDimensionalPosition()) {
                mpLeft = std::make_unique<TreeNodeType>(Begin + 1, End, mDimension + 1);
            } else {
                mpRight = std::make_unique<TreeNodeType>(Begin + 1, End, mDimension + 1);
            }
        } else {
            // first sort based on the dimension
            std::sort(Begin, End, [&](const auto &rA, const auto &rB)
                      { return TPositionGetterFunctor::GetPosition(rA)[mDimension] < TPositionGetterFunctor::GetPosition(rB)[mDimension]; });

            // now get the mid point of that dimension
            const IndexType mid_point_index = number_of_items / 2;
            mData = *(Begin + mid_point_index);

            // now create left and right tree node
            mpLeft = std::make_unique<TreeNodeType>(Begin, Begin + mid_point_index, mDimension + 1);
            mpRight = std::make_unique<TreeNodeType>(Begin + mid_point_index + 1, End, mDimension + 1);
        }
    }

    ~TreeNode() = default;

    ///@}
    ///@name Public operations
    ///@{

    TDataType GetData() const
    {
        return mData;
    }

    std::string Info(const std::string& Tabbing) const
    {
        std::stringstream info;
        info << "\n" << Tabbing << "Data: " << *mData << ", Dimension: " << mDimension;
        info << "\n" << Tabbing << "Left : " << (mpLeft ? mpLeft->Info(Tabbing + "   ") : "NONE");
        info << "\n" << Tabbing << "Right: " << (mpRight ? mpRight->Info(Tabbing + "   ") : "NONE");
        return info.str();
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const IndexType mDimension;

    TDataType mData;

    std::unique_ptr<TreeNodeType> mpLeft;

    std::unique_ptr<TreeNodeType> mpRight;

    ///@}
    ///@name Private methods
    ///@{

    PositionType GetPosition() const
    {
        return TPositionGetterFunctor::GetPosition(mData);
    }

    double GetDimensionalPosition() const
    {
        return GetPosition()[mDimension];
    }

    double GetDimensionalDistance(const PositionType& rPosition) const
    {
        return rPosition[mDimension] - GetDimensionalPosition();
    }

    double GetDistanceSquare(const PositionType& rPosition) const
    {
        const auto& position = GetPosition();
        double distance_square = 0;
        for (IndexType i = 0; i < TDim; ++i){
            distance_square += std::pow(rPosition[i] - position[i], 2);
        }
        return distance_square;
    }

    void FillTreeNodesWithinDistance(
        std::vector<TreeNodeType::Pointer>& rTreeNodes,
        IndexType& rNumberOfNeighbours,
        const PositionType& rPosition,
        const double DistanceSquare) const
    {
        // exit the recursion if enough number of neighbours found.
        if (rNumberOfNeighbours >= rTreeNodes.size()) {
            return;
        }

        Counter++;
        const double dimensional_distance = GetDimensionalDistance(rPosition);
        // std::cout << "\ndebug detail: " << this->mData->GetId() << ", dimensional distance = " << dimensional_distance << ", dimension = " << mDimension;
        if (dimensional_distance * dimensional_distance <= DistanceSquare) {
            if (GetDistanceSquare(rPosition) <= DistanceSquare) {
                rTreeNodes[rNumberOfNeighbours++] = this;
            }
            if (mpLeft) {
                mpLeft->FillTreeNodesWithinDistance(rTreeNodes, rNumberOfNeighbours, rPosition, DistanceSquare);
            }
            if (mpRight) {
                mpRight->FillTreeNodesWithinDistance(rTreeNodes, rNumberOfNeighbours, rPosition, DistanceSquare);
            }
        } else {
            if (mpLeft && dimensional_distance < 0) {
                mpLeft->FillTreeNodesWithinDistance(rTreeNodes, rNumberOfNeighbours, rPosition, DistanceSquare);
            } else if (mpRight && dimensional_distance >= 0) {
                mpRight->FillTreeNodesWithinDistance(rTreeNodes, rNumberOfNeighbours, rPosition, DistanceSquare);
            }
        }
    }

    ///@}
    ///@name Friends
    ///@{

    friend class KDTree<TDim, TDataType, TPositionGetterFunctor>;

    ///@}
};

template<unsigned int TDim, class TDataType, class TPositionGetterFunctor>
std::size_t TreeNode<TDim, TDataType, TPositionGetterFunctor>::Counter = 0;

template<unsigned int TDim, class TDataType, class TPositionGetterFunctor>
class KDTree
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using TreeNodeType = TreeNode<TDim, TDataType, TPositionGetterFunctor>;

    using PositionType = typename TreeNodeType::PositionType;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TIdIteratorType>
    KDTree(
        TIdIteratorType Begin,
        TIdIteratorType End)
    {
        std::cout <<"\nConstructing KDTree...";
        const auto start = std::chrono::high_resolution_clock::now();
        mpRootNode = std::make_unique<TreeNodeType>(Begin, End, 0);
        const auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\nConstructed KDTree witin " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " us.";
    }

    ///@}
    ///@name Public operations
    ///@{

    IndexType FindNearestNeighbours(
        std::vector<typename TreeNodeType::Pointer>& rNearesetTreeNodes,
        const PositionType& rPosition,
        const double DistanceTolerance) const
    {
        IndexType number_of_neighbours = 0;
        mpRootNode->FillTreeNodesWithinDistance(
            rNearesetTreeNodes,
            number_of_neighbours,
            rPosition,
            std::pow(DistanceTolerance, 2));
        return number_of_neighbours;
    }

    std::string Info() const
    {
        return mpRootNode->Info("");
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::unique_ptr<TreeNodeType> mpRootNode;

    ///@}
};