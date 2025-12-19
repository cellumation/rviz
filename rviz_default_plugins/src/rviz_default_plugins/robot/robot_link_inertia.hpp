#ifndef RVIZ_DEFAULT_PLUGINS__ROBOT__ROBOT_LINK_INERTIA_HPP_
#define RVIZ_DEFAULT_PLUGINS__ROBOT__ROBOT_LINK_INERTIA_HPP_

#include <Eigen/Geometry>

#include <optional>

namespace rviz_default_plugins::robot
{
struct InertiaEquivalentBoxResult
{
  Eigen::Vector3d box_size;
  Eigen::Quaterniond rotation;
};

std::optional<InertiaEquivalentBoxResult>
inertiaEquivalentBox(double mass, const Eigen::Matrix3d & inertia, double tolerance = 1e-6);
}  // namespace rviz_default_plugins::robot

#endif  // RVIZ_DEFAULT_PLUGINS__ROBOT__ROBOT_LINK_INERTIA_HPP_
