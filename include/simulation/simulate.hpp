#pragma once


#include "../../include/camera/algorithms.hpp"
#include "../../include/drone/trajectory.hpp"
#include "../../include/drone/time_based_traj.hpp"

#include "../../include/json_convert/drone_to_json.hpp"
#include "../../include/json_convert/convex_polyhedra_to_json.hpp"

#include <queue>


/*==========DRONE'S CONSTANTS==========*/
const double TIMESTEP = 1.0/30.0; /* (1/sample frequency), sample freq max is PARTIAL_PATH_IN_EACH_SEGMENT */
const int DRONE_INITIAL_SPEED = 0; /* m/s ,in case velocity = const*/
const int DRONE_INITIAL_ACC = 12;
const int DRONE_MAX_SPEED = 40;
const double DRONE_RADIUS = 1.0;

/*==========TRAJECTORY CONSTANTS==========*/
const double CYLYNDER_RADIUS = 50.0;      // meters
const int INTERPOLATION_POINTS = 10;       // q0 + 8 random + qn

const Eigen::Vector3d STARTING_POINT(700, 120,350);
// const Eigen::Vector3d ENDING_POINT(200,20,-150); /*ending point*/
const Eigen::Vector3d ENDING_POINT(-250,60,-250);

const double ACC_DEVIATION = 5.0; /*from 0.1 to 2.0 m/s^2*/ /*q_x = q_y= q_z = ACC_DEVIATION * ACC_DEVIATION*/

/*==========FILTER"S CONSTANTS==========*/
const int REACQUIRE_THRESHOLD_FRAMES = 30;
const int TRACKING_THRESHOLD_FRAMES = 5;
const double ALG_PROCESSING_TIMESTEP = TIMESTEP;
const double ROT_CAM_LEADING_TIME = TIMESTEP * 6;
const int PF_PARTICLES_NUMBER = 1000;

const int LATENCY = 3; /* 3 frames late*/
const int LEAD_FRAMES_NUMBER = 6; /*filter estimates 6 frames forward*/

/*For plot reset frames info*/
struct LostObjFrame {
    int frame;
    Eigen::Vector3d pos;
    std::vector<Pixel_Obj> list;
};

struct LostTimePerFrameInfo{
    int reacquire_count = 0;
    int reacquire_cam_numbers;
} ;
struct LostFrameInfo {
    int frame;
    int lost_count = 0;
    std::vector<LostTimePerFrameInfo> lost_time_info;
    bool isSuccess; /*False -> lost obj | true -> reacquire succesfully*/
};
/*========================================================================*/

/*==========RESULT REPRESENTATION==========*/
struct Rmse_based_frame {
    int frame;
    double rmse;
    int rmse_count = 1;
};

/*==========SIMULATION==========*/

enum TrackingState
{
    RESET_STATE,
    TRACKING_STATE,
    REACQUIRE_STATE
};

enum ReacquireMode {
    SEARCH_MODE,
    FOCUS_MODE
};

struct SimulationResult {
    double pf_penalty_count = 0.0;
    double number_recover_from_reacquire = 0.0;
    double successful_redetection_percentage = 0.0;
    std::vector<Rmse_based_frame> rmse_list_tracking;
    std::vector<Rmse_based_frame> rmse_list_det_succ;
    std::vector<Rmse_based_frame> rmse_list_det_fail;
    
};

class Simulation {
    private:
        std::vector<Camera> cameras;
        SimulationResult result;

        /*Tracking variable*/
        int current_frame = 0;
        int total_frames;

        ReacquireMode reacquire_mode = SEARCH_MODE;

        // lưu hypothesis đang được focus
        Eigen::Vector3d focus_target;

        /*debug parameter*/
        std::vector<LostObjFrame> lost_frame_list; 

    public:
        Simulation();
        std::vector<LostFrameInfo> lost_frame_list_info; /*Information of each obj list*/
        SimulationResult getResult(){return result;}
        void do_simulate();
        void show_result();
        std::vector<Rmse_based_frame> temp_det_rmse_list;
};

void save_rmse_to_json(const std::vector<Rmse_based_frame>& rmse_history,
                       const std::string& filename, const int latency);