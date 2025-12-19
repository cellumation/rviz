#include <Eigen/Geometry>

#include <gtest/gtest.h>

#include "../../../src/rviz_default_plugins/robot/robot_link_inertia.hpp"

using namespace ::testing;  // NOLINT
using namespace rviz_default_plugins::robot;  // NOLINT

static Eigen::Matrix3d calculateBoxInertia(
  const double mass, const Eigen::Vector3d & size, const Eigen::Quaterniond & rotation)
{
  assert(std::isfinite(mass) && mass >= 0.);
  assert(size.allFinite() && (size.array() >= 0.).all());
  assert(rotation.coeffs().allFinite() && rotation.normalized().isApprox(rotation));

  const Eigen::Vector3d sqrSize{size.cwiseAbs2()};
  const Eigen::Vector3d principalComponents{mass / 12. * Eigen::Vector3d{sqrSize.y() + sqrSize.z(),
      sqrSize.z() + sqrSize.x(), sqrSize.x() + sqrSize.y()}};
  const Eigen::Matrix3d rotationMatrix{rotation.toRotationMatrix()};
  return rotationMatrix * principalComponents.asDiagonal() * rotationMatrix.transpose();
}

template<typename Mat1, typename Mat2>
static testing::AssertionResult isApprox(
  const char * m1_expr, const char * m2_expr, const char * /*prec_expr*/,
  const Eigen::MatrixBase<Mat1> & m1, const Eigen::MatrixBase<Mat2> & m2,
  const typename Eigen::MatrixBase<Mat1>::Scalar prec)
{
  if (m1.isApprox(m2, prec)) {
    return testing::AssertionSuccess();
  }
  const Eigen::IOFormat fmt{Eigen::FullPrecision, 0, ", ", ";\n", "", "", "[", "]"};
  return testing::AssertionFailure() << "Expected \"" << m1_expr << "\" which is: \n" <<
         m1.format(fmt) <<
         "\nto be approximately equal to \"" << m2_expr << "\", which is: \n" <<
         m2.format(fmt) <<
         "\nwhere the norm of difference is: " << (m1 - m2).norm();
}

template<typename T>
static testing::AssertionResult isApproxBoxRotation(
  const char * q1_expr, const char * q2_expr, const char * /*prec_expr*/,
  const Eigen::Quaternion<T> & q1, const Eigen::Quaternion<T> & q2, const T prec)
{
  if (q1.isApprox(q2, prec)) {
    return testing::AssertionSuccess();
  }
  const Eigen::Vector3d axis{Eigen::AngleAxisd{q2}.axis()};
  const Eigen::Quaterniond altq2{Eigen::AngleAxisd{M_PI, axis} *q2};
  if (q1.isApprox(altq2)) {
    return testing::AssertionSuccess();
  }
  return testing::AssertionFailure() << "Expected \"" << q1_expr << "\", which is:\n" << q1 <<
         ",\nto be either approximately equal to \"" << q2_expr << "\", which is: \n" << q2 <<
         ",\nor (due to 180° box symmetry): " << altq2;
}

TEST(RobotTest, inertiaEquivalentBox)
{
  EXPECT_FALSE(inertiaEquivalentBox(-1., Eigen::Matrix3d::Zero()));
  EXPECT_FALSE(inertiaEquivalentBox(0., Eigen::Matrix3d::Zero()));
  EXPECT_FALSE(inertiaEquivalentBox(1., Eigen::Matrix3d::Zero()));
  EXPECT_FALSE(inertiaEquivalentBox(-1., Eigen::Matrix3d::Identity()));
  EXPECT_FALSE(inertiaEquivalentBox(0., Eigen::Matrix3d::Identity()));

  { // Moment of inertia matrix that doesn't satisfy triangle inequality
    const Eigen::Matrix3d inertia{(Eigen::Matrix3d{} << 2, -1, 0, -1, 2, -1, 0, -1, 2).finished()};
    EXPECT_FALSE(inertiaEquivalentBox(1., inertia));
  }
  { // Expect cube with side lengths sqrt(6) from identity inertia matrix
    const auto result{inertiaEquivalentBox(1., Eigen::Matrix3d::Identity())};
    EXPECT_TRUE(result);
    EXPECT_PRED_FORMAT3(
      isApprox, result.value().box_size, Eigen::Vector3d::Constant(std::sqrt(6)), 1e-9);
    EXPECT_PRED_FORMAT3(
      isApproxBoxRotation, result.value().rotation, Eigen::Quaterniond::Identity(), 1e-9);
  }
  { // box [1, 1, 1]
    const double mass{1.};
    const Eigen::Vector3d size(Eigen::Vector3d::Ones());
    const Eigen::Quaterniond rotation{Eigen::Quaterniond::Identity()};
    const auto result{inertiaEquivalentBox(mass, calculateBoxInertia(mass, size, rotation))};
    EXPECT_TRUE(result);
    EXPECT_TRUE(result.value().box_size.isApprox(size));
    EXPECT_TRUE(result.value().rotation.isApprox(rotation));
  }
  { // box [8, 4, 2]
    const double mass{16.};
    const Eigen::Vector3d size{8., 4., 2.};
    const Eigen::Quaterniond rotation{Eigen::Quaterniond::Identity()};
    const auto result{inertiaEquivalentBox(mass, calculateBoxInertia(mass, size, rotation))};
    EXPECT_TRUE(result);
    EXPECT_PRED_FORMAT3(isApprox, result.value().box_size, size, 1e-9);
    EXPECT_PRED_FORMAT3(isApproxBoxRotation, result.value().rotation, rotation, 1e-9);
  }
  { // box [8, 4, 1e5]
    const double mass{16.};
    const Eigen::Vector3d size{8., 4., 1e5};
    const Eigen::Quaterniond rotation{Eigen::Quaterniond::Identity()};
    const auto result{inertiaEquivalentBox(mass, calculateBoxInertia(mass, size, rotation))};
    EXPECT_TRUE(result);
    EXPECT_PRED_FORMAT3(isApprox, result.value().box_size, size, 1e-9);
    EXPECT_PRED_FORMAT3(isApproxBoxRotation, result.value().rotation, rotation, 1e-9);
  }
  { // box [4, 1e8, 7]
    const double mass{16.};
    const Eigen::Vector3d size{4., 1e8, 7.};
    const Eigen::Quaterniond rotation{Eigen::Quaterniond::Identity()};
    const auto result{inertiaEquivalentBox(mass, calculateBoxInertia(mass, size, rotation))};
    EXPECT_TRUE(result);
    EXPECT_PRED_FORMAT3(isApprox, result.value().box_size, size, 1e-9);
    EXPECT_PRED_FORMAT3(isApproxBoxRotation, result.value().rotation, rotation, 1e-9);
  }
  { // box [9, 4, 1] rotated by 90° around Z, equivalent to [4, 9, 1] after rotation
    const double mass{18.};
    const Eigen::Vector3d size{9., 4., 1.};
    const Eigen::Quaterniond rotation{Eigen::AngleAxisd{M_PI / 2, Eigen::Vector3d::UnitZ()}};
    const auto result{inertiaEquivalentBox(mass, calculateBoxInertia(mass, size, rotation))};
    EXPECT_TRUE(result);
    const Eigen::Vector3d rotatedSize{size.y(), size.x(), size.z()};
    EXPECT_PRED_FORMAT3(isApprox, result.value().box_size, rotatedSize, 1e-9);
    EXPECT_PRED_FORMAT3(
      isApproxBoxRotation, result.value().rotation, Eigen::Quaterniond::Identity(), 1e-9);
  }
  { // box [14, 9, 5] rotated by 45° around X
    const double mass{68.};
    const Eigen::Vector3d size{14., 9., 5.};
    const Eigen::Quaterniond rotation{Eigen::AngleAxisd{M_PI / 4, Eigen::Vector3d::UnitX()}};
    const auto result{inertiaEquivalentBox(mass, calculateBoxInertia(mass, size, rotation))};
    EXPECT_TRUE(result);
    EXPECT_PRED_FORMAT3(isApprox, result.value().box_size, size, 1e-9);
    EXPECT_PRED_FORMAT3(isApproxBoxRotation, result.value().rotation, rotation, 1e-9);
  }
  { // box [11, 3, 2] rotated by 22.5° around axis [0, 1, 1].normalized
    const double mass{42.};
    const Eigen::Vector3d size{11., 3., 2.};
    const Eigen::Quaterniond rotation{
      Eigen::AngleAxisd{M_PI / 8, Eigen::Vector3d{0, 1, 1}.normalized()}};
    const Eigen::Matrix3d expectedInertia{calculateBoxInertia(mass, size, rotation)};
    const auto result{inertiaEquivalentBox(mass, expectedInertia)};
    EXPECT_TRUE(result);
    EXPECT_PRED_FORMAT3(isApprox, result.value().box_size, size, 1e-9);
    // The rotation may not be the same but the resulting inertia matrix has to.
    EXPECT_PRED_FORMAT3(
      isApprox, expectedInertia,
      calculateBoxInertia(mass, result.value().box_size, result.value().rotation), 1e-9);
  }
}
