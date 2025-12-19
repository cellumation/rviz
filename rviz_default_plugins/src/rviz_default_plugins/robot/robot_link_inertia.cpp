#include "robot_link_inertia.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <cmath>

namespace rviz_default_plugins::robot
{
static Eigen::Vector3d calculateSqrScaledSize(const Eigen::Vector3d & principalMoments)
{
  return Eigen::Vector3d{principalMoments.y() + principalMoments.z(),
    principalMoments.z() + principalMoments.x(),
    principalMoments.x() + principalMoments.y()} - principalMoments;
}

static bool validPrincipalMoments(const Eigen::Vector3d & principalMoments, const double tolerance)
{
  if ((principalMoments.array() < 0).any() && !principalMoments.isZero(tolerance)) {
    return false;
  }
  const Eigen::Vector3d sqrScaledSize{calculateSqrScaledSize(principalMoments)};
  return sqrScaledSize.allFinite() &&
         ((sqrScaledSize.array() > 0).all() || sqrScaledSize.isZero(tolerance));
}

static Eigen::Matrix3d makeSymmetric(Eigen::Matrix3d matrix)
{
  matrix(0, 1) = matrix(1, 0) = (matrix(0, 1) + matrix(1, 0)) / 2;
  matrix(0, 2) = matrix(2, 0) = (matrix(0, 2) + matrix(2, 0)) / 2;
  matrix(1, 2) = matrix(2, 1) = (matrix(1, 2) + matrix(2, 1)) / 2;
  return matrix;
}

std::optional<InertiaEquivalentBoxResult> inertiaEquivalentBox(
  const double mass, const Eigen::Matrix3d & inertia, const double tolerance)
{
  if (!std::isfinite(mass) || !inertia.allFinite() || mass <= 0 || inertia.determinant() <= 0) {
    return std::nullopt;
  }

  const auto boxSizeFromPrincipalComponents{[&](const Eigen::Vector3d & principalMoments) {
      return (6 / mass * calculateSqrScaledSize(principalMoments.cwiseMax(0))).cwiseSqrt().eval();
    }};
  if (inertia.isDiagonal(tolerance)) {
    if (!validPrincipalMoments(inertia.diagonal(), tolerance)) {
      return std::nullopt;
    }
    return InertiaEquivalentBoxResult{boxSizeFromPrincipalComponents(inertia.diagonal()),
      Eigen::Quaterniond::Identity()};
  }

  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver{makeSymmetric(inertia)};
  if (eigen_solver.info() != Eigen::ComputationInfo::Success) {
    return std::nullopt;
  }

  if (!validPrincipalMoments(eigen_solver.eigenvalues(), tolerance)) {
    return std::nullopt;
  }

  // Make sure the eigen vectors build a right handed coordinate system.
  Eigen::Matrix3d eigen_vectors{eigen_solver.eigenvectors()};
  if (eigen_vectors.determinant() < 0) {
    eigen_vectors.col(2) = -eigen_vectors.col(2);
  }

  return InertiaEquivalentBoxResult{boxSizeFromPrincipalComponents(eigen_solver.eigenvalues()),
    Eigen::Quaterniond{eigen_vectors}};
}
}  // namespace rviz_default_plugins::robot
