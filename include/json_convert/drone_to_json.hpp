#include <fstream>
#include <Eigen/Dense>
#include "../drone/trajectory.hpp"
#include "../drone/time_based_traj.hpp"

/*================== TRAJECTORY ====================*/

void trajectory_to_json(
const std::vector<Eigen::Vector3d>& points,
const std::vector<Cubic_Interpolation_Coeff>& sx,
const std::vector<Cubic_Interpolation_Coeff>& sy,
const std::vector<Cubic_Interpolation_Coeff>& sz,
const std::string& filename);

/*============== TIME-BASED TRAJECTORY ===============*/
void drone_motion_to_json(Time_Based_Traj time_traj, const std::string& filename);

void export_trajectories_to_json(
    const std::vector<Eigen::Vector3d>& est_traj, 
    const std::string& filename);