#pragma once

#include "trajectory.hpp"
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

const size_t PARTIAL_PATH_IN_EACH_SEGMENT = 300; /*number points of each segment  */
const size_t MAX_SPEED = 25; /*100 m/s = 360 km/h*/

/*U is interpolated param which is used to interpolated trajectory*/
class Time_Based_Traj {
    public:
    explicit Time_Based_Traj(const Trajectory& traj);

    // set velocity
    void setVelocity(double v);
    void setTimeStep(double dt);
    void setAcceleration(double a);
    void setMaxSpeed(double v_max); 

    // update in time
    void update();

    void reset();
    bool finished() const;


    // getters
    Eigen::Vector3d position() const;
    Eigen::Vector3d velocity() const;
    double getTimeStep(){return dt_;};

    double getTotalLength() const;

    private:
    // spline segment
    struct Segment {
        Cubic_Interpolation_Coeff sx;
        Cubic_Interpolation_Coeff sy;
        Cubic_Interpolation_Coeff sz;
        double u_max;
        double length;
    };

private:
    // compute arc-length using simpson approximation
    double arcLengthSimpson(
        const Segment& seg,
        double u_max,
        int N = PARTIAL_PATH_IN_EACH_SEGMENT
    ) const;

    /*Find interpolated value using bisection after having path_length*/
    double findUfromArc(
        const Segment& seg,
        double target_s
    ) const;

    /*Find segment and interpolated value.
    This function includes findUfromArc()*/
    void findSegmentAndU(
        double s,
        std::size_t& seg_id,
        double& u
    ) const;

private:
    std::vector<Segment> segments_;
    double total_length_ = 0.0;

    // trạng thái thời gian
    double s_ = 0.0;
    double speed_ = 0.0;
    double accel_ = 0.0;
    double max_speed_ = MAX_SPEED;  
    double dt_ = 1.0/30.0; /*60Hz*/
};


/*for simulation - read data from file trajectory json*/
struct FrameData {
    double t;
    Eigen::Vector3d pos;
    Eigen::Vector3d vel;   // optional, nếu chưa cần có thể bỏ
};

struct TrajectoryData {
    double dt;
    std::vector<FrameData> frames;
};

TrajectoryData load_trajectory_from_json(const std::string& filename);

std::vector<FrameData> load_frames_from_json(const std::string& filename);
