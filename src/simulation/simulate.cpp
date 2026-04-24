#include "../../include/simulation/simulate.hpp"

/*===========================================================================*/
/*===========================================================================*/

/*restore angles after scanning*/
void point_all_cameras_to_object(std::vector<Camera>& cameras, const Eigen::Vector3d& obj_world) {
    for (auto& cam : cameras) {
        Eigen::Vector3d cam_pos = cam.getPosition();

        Eigen::Vector3d direction = obj_world - cam_pos;

        double target_az_rad = std::atan2(direction.x(), direction.z());

        double horizontal_dist = std::sqrt(direction.x() * direction.x() + direction.z() * direction.z());
        
        double target_el_rad = std::atan2(direction.y(), horizontal_dist);
        cam.set_camera_angles(target_az_rad, target_el_rad);
    }
}

void export_pf_json(
    std::vector<Polyhedron>& fov_polys,
    std::vector<Particle> particles,
    const Point& obj_pos,
    const Point& est_pos,
    const std::string& filename
) {
    json J;

    // ===== polyhedra =====
    J["polyhedra"] = json::object();
    int i = 1;
    for(Polyhedron pol : fov_polys){
        std::string cam_id = std::to_string(i++); 
        polyhedron_to_json(pol, J["polyhedra"][cam_id]);
    }

    J["pos"] = { 
        CGAL::to_double(obj_pos.x()), 
        CGAL::to_double(obj_pos.y()), 
        CGAL::to_double(obj_pos.z()) 
    };

    J["est_pos"] = { 
        CGAL::to_double(est_pos.x()), 
        CGAL::to_double(est_pos.y()), 
        CGAL::to_double(est_pos.z()) 
    };
    
    J["particles"] = json::array();

    for (const auto& p : particles) {
        json p_json;

        p_json["pos"] = {
            p.x(0),  // x
            p.x(3),  // y
            p.x(6)   // z
        };

        p_json["weight"] = p.weight;

        J["particles"].push_back(p_json);
    }

    std::ofstream(filename) << J.dump(2);
}

void save_rmse_to_json(const std::vector<Rmse_based_frame>& rmse_history,
                       const std::string& filename, const int latency)
{
    json j;
    j["rmse_per_frame"] = json::array();

    for (size_t i = 0; i < rmse_history.size(); ++i)
    {
        json frame_entry;
        frame_entry["frame"] = rmse_history[i].frame;         
        frame_entry["rmse"]  = rmse_history[i].rmse;

        j["rmse_per_frame"].push_back(frame_entry);
    }

    std::ofstream file(filename);
    if (file.is_open())
    {
        file << j.dump(4); 
        file.close();
    }
}

Simulation::Simulation(){
     try {
        /*================================================================*/
        /*======================== SETUP CAMERA ==========================*/
        /*================================================================*/

        Camera Cam1 = Camera::createFromFile("./config/camera1_config.txt");
        Camera Cam2 = Camera::createFromFile("./config/camera2_config.txt");
        Camera Cam3 = Camera::createFromFile("./config/camera3_config.txt");
        Camera Cam4 = Camera::createFromFile("./config/camera4_config.txt");
        Camera Cam5 = Camera::createFromFile("./config/camera5_config.txt");
        Camera Cam6 = Camera::createFromFile("./config/camera6_config.txt");
        Camera Cam7 = Camera::createFromFile("./config/camera7_config.txt");
        Camera Cam8 = Camera::createFromFile("./config/camera8_config.txt");
        std::cout << "Khởi tạo camera thành công!" << std::endl;
        Cam1.set_camera_angles_initialized(66.37 , 15.83);
        Cam2.set_camera_angles_initialized(64.69 , 16.82);
        Cam3.set_camera_angles_initialized(68.6 , 17.29);
        Cam4.set_camera_angles_initialized(70.07 , 16.21);
        Cam5.set_camera_angles_initialized(64.51 , 14.57);
        Cam6.set_camera_angles_initialized(58.91 , 17.30);
        Cam7.set_camera_angles_initialized(71.31 , 19.02);
        Cam8.set_camera_angles_initialized(75.03 , 15.54);

        cameras.push_back(Cam1);
        cameras.push_back(Cam2);
        cameras.push_back(Cam3);
        cameras.push_back(Cam4);
        // cameras.push_back(Cam5);
        // cameras.push_back(Cam6);
        // cameras.push_back(Cam7);
        // cameras.push_back(Cam8);
        /*================================================================*/
        /*===================== BUILD TRAJECTORY =========================*/
        /*================================================================*/
        /*build for the first time, then comment this section*/

         Trajectory traj(
            STARTING_POINT,
            ENDING_POINT,
            CYLYNDER_RADIUS,
            INTERPOLATION_POINTS
        );

        traj.build();

        Time_Based_Traj time_based_trajectory(traj);
        time_based_trajectory.setTimeStep(TIMESTEP);
        time_based_trajectory.setVelocity(DRONE_INITIAL_SPEED);
        time_based_trajectory.setAcceleration(DRONE_INITIAL_ACC);
        time_based_trajectory.setMaxSpeed(DRONE_MAX_SPEED);

        trajectory_to_json(
            traj.getInterpolationPoints(),
            traj.getSx(),
            traj.getSy(),
            traj.getSz(),
            "trajectory.json"
        ); 

        drone_motion_to_json(time_based_trajectory, "trajectory_motion.json");

    } catch (const std::exception& e) {
        std::cerr << "Lỗi: " << e.what() << std::endl;
    }
}

void Simulation::do_simulate(){

    /*================================================================*/
    /*==================== GET DRONE"S POSITION ======================*/
    /*================================================================*/

    std::cout << "start do_simulate ()\n";
    std::vector<FrameData> traj = load_frames_from_json("trajectory_motion.json");
    total_frames  = traj.size();
    /*variable for debug */
    std::vector<Eigen::Vector3d> result_traj; /* for drawing*/

    /*we build fov of a camera changes. In this simulation, that feature is set*/
    for(auto& cam : cameras){
        cam.initial_fov();
    }

    ParticleFilter pf(PF_PARTICLES_NUMBER); /*gen 500 particles*/
    UnscentedKalmanFilter ukf(ACC_DEVIATION, MEASUREMENT_DEVIATION);
    bool is_initialized = false;  /* for initializing ukf state and pf state*/
    int one_cam_seen_obj_frame_count = 0; /*Scan if only one camera sees object continously*/

    bool is_intersection_initial_plot = false; /*After compute the intersection of fovs in the 1st loop, this value becomes true*/

    /*For delayed traj frame*/
    std::deque<Eigen::Vector3d> traj_buffer;

    /*state control*/
    TrackingState state = TRACKING_STATE;

    Eigen::Vector3d current_position(0,0,0); /* (current frame -- LATENCY) position*/

    /*For debugging*/
    bool is_json_pf = false;
    int count_pf = 0;
    int reacquire_frame_count = 0; /*For remaining stable status after reacquire successfully*/
    bool is_require_state_initialized = false;
    int focus_time_per_frame_count = 0;
    LostFrameInfo lost_frame_info; /*For debug lost frame */

    point_all_cameras_to_object(cameras, traj[0].pos);

    while(current_frame < total_frames){

        /*Delay system*/
        // traj_buffer.push_back(traj[current_frame].pos);
        // if(traj_buffer.size() > LATENCY){
        //     current_position = traj_buffer.front();
        //     traj_buffer.pop_front();
        // } else {
        //     current_frame++;
        //     continue;
        // }
        
        current_position = traj[current_frame].pos;

        /*================================================================*/
        /*======================== PREPROCESSING =========================*/
        /*================================================================*/
        std::vector<Camera*> obj_visible_cameras; /*list of cameras which can see obj*/
        std::vector<Camera*> obj_hidden_cameras; /* list of camera which can not see obj*/
        std::vector<Camera*> pf_cameras;

        for (auto& cam : cameras) {
            
            /*compute min and max distance from each camera to trajectory (for showing results)*/
            double distance = (traj[current_frame].pos - cam.getPosition()).norm();
            if(distance < cam.min_distance) cam.min_distance = distance;
            if(distance > cam.max_distance) cam.max_distance = distance;
            
            /*project obj on screen*/
            cam.calculate_pixel_projection(current_position, DRONE_RADIUS);

            pf_cameras.push_back(&cam);
            if (cam.obj_is_visible) {     
                obj_visible_cameras.push_back(&cam);
                cam.obj_visible_frame_count++;

            } else {
                obj_hidden_cameras.push_back(&cam);
            } 
        }

        int cameras_seen_obj = obj_visible_cameras.size();

        if(cameras_seen_obj <= 1) {
            one_cam_seen_obj_frame_count++;
        } else
        {
            one_cam_seen_obj_frame_count = 0;
        }

        /*==============================================================*/
        /*====================ALGORITHM SIMULATION======================*/
        /*==============================================================*/
        
        switch(state)
        {
            /*=============================================================*/
            /* ======== SCANNING CASE (ALL CAMERAS LOST OBJ'S SIGHT)=======*/
            /*=============================================================*/
            case RESET_STATE:
            {
            
                /*Print lost frame list on console for debugging*/
                std::cout << "[DEBUG] Lost target at frame " << current_frame << std::endl;
                // std::cout << "Real position: (" << traj[current_frame].pos[0] << ", " 
                //             << traj[current_frame].pos[1] << ", "
                //             << traj[current_frame].pos[2] << ")\n";
                // LostObjFrame lost_frame;
                // lost_frame.frame = current_frame;
                // lost_frame.pos = traj[current_frame].pos;

                // for(auto& cam : cameras) {
                //     cam.calculate_pixel_projection(traj[current_frame].pos, DRONE_RADIUS);
                //     // std::cout << "cam " << cam.getID() << " pos: (" << cam.get_camera_angles().azimuth * 180 / M_PI
                //     //             << ", " << cam.get_camera_angles().elevation * 180 / M_PI<< ") | (u,v) = (" 
                //     //     << cam.pix_info.u << ", " << cam.pix_info.v <<")\n";
                //     Pixel_Obj pix = cam.pix_info;
                //     lost_frame.list.push_back(pix);
                // }
                // lost_frame_list.push_back(lost_frame);

                /*RECOVERY AFTER SCANNING*/
                result.pf_penalty_count++; /*scanning time*/
                is_initialized = false;
                point_all_cameras_to_object(cameras, traj[current_frame].pos);
                state = TRACKING_STATE;
                break;
            }

            /*============================================================ */
            /*===========================UKF HERE========================= */
            /*============================================================ */
            case TRACKING_STATE:
            {
                if (!is_initialized) {
                    Polyhedron init_poly = calculate_total_intersection(obj_visible_cameras);
                    if (!init_poly.is_empty()) {
                        Point start_pt = find_volume_centroid(init_poly);
                        Eigen::VectorXd init_state = Eigen::VectorXd::Zero(9);
                        init_state(0) = CGAL::to_double(start_pt.x()); 
                        init_state(3) = CGAL::to_double(start_pt.y()); 
                        init_state(6) = CGAL::to_double(start_pt.z());
                        ukf.setState(init_state);
                        ukf.setCovariance(Eigen::MatrixXd::Identity(9, 9) * 1000.0);
                    }
                    is_initialized = true;
                }

                ukf.predict(ALG_PROCESSING_TIMESTEP); /*answer "if there is no n */                                    
                if(cameras_seen_obj > 0){
                    ukf.update(obj_visible_cameras); /*best estimation right now & answer "what measurements think about state we've predicted yet"*/
                }

                GaussState lead;
                lead.mean = ukf.getState();
                lead.covariance = ukf.getCovariance();
                // Early frame → lead 1 frame
                // if (is_recently_initialized) {
                    // lead = ukf.getLeadState(TIMESTEP, ACC_DEVIATION);
                // }
                // // Late frame → lead latency camera
                // else {
                //     lead = ukf.getLeadState(TIMESTEP*LEAD_FRAMES_NUMBER, ACC_DEVIATION);
                // }

                Polyhedron target_intersection = calculate_total_intersection(obj_visible_cameras);

                // optimize_all_cameras_jointly_weighted(obj_visible_cameras,lead, 0.001); /*Rotate cameras*/
                double px_lead = lead.mean(0);
                double py_lead = lead.mean(3);
                double pz_lead = lead.mean(6);
                Eigen::Vector3d lead_pos(px_lead,py_lead,pz_lead);
                rotate_cameras(obj_visible_cameras, lead_pos);
                rotate_cameras(obj_hidden_cameras,lead_pos);

                // save real trajectory (no leading) for plotting
                Eigen::VectorXd x_cur = ukf.getState();
                double px = x_cur(0);
                double py = x_cur(3);
                double pz = x_cur(6);
                Eigen::Vector3d estimate_pos = Eigen::Vector3d(px, py, pz); 
                result_traj.push_back(estimate_pos);

                /*Propagating estimated postition forward so as to compute rmse*/
                // GaussState rmse_state = ukf.getLeadState(LATENCY*TIMESTEP, ACC_DEVIATION);
                GaussState rmse_state;
                rmse_state.mean = ukf.getState();
                rmse_state.covariance = ukf.getCovariance();

                double rmse_x = rmse_state.mean(0);
                double rmse_y = rmse_state.mean(3);
                double rmse_z = rmse_state.mean(6);
                Eigen::Vector3d rmse_pos(rmse_x,rmse_y,rmse_z);
                Eigen::Vector3d error = rmse_pos - current_position;
                Rmse_based_frame rmse_frame;
                rmse_frame.frame = current_frame;rmse_frame.rmse = error.norm();
                result.rmse_list_tracking.push_back(rmse_frame);
                
                /*export json for plotting FOV of cameras information at the beginning*/
                if(!is_intersection_initial_plot){
                    const Point starting_p(STARTING_POINT(0),STARTING_POINT(1),STARTING_POINT(2));
                    std::vector<Polyhedron> fov_polys;
                    for(auto &cam: cameras){
                        fov_polys.push_back(cam.get_FOV()); 
                    }
                    export_intersection_json(fov_polys, target_intersection, starting_p, "intersect.json");
                    is_intersection_initial_plot = true;
                }

                /*print camera tracking state on console*/
                std::cout << "[" << current_frame << "] UKF mean: (" << px << ", " << py << ", " << pz << ")" ;
                std::cout << " | Cams: " << obj_visible_cameras.size() <<  std::endl;
                std::cout << "Real position: (" << traj[current_frame].pos[0] << ", " 
                            << traj[current_frame].pos[1] << ", "
                            << traj[current_frame ].pos[2] << ") | " 
                            << "rmse is " << rmse_frame.rmse << "\n";

                for(auto cam_print : obj_visible_cameras){
                    Camera* cam_tes = cam_print; 
                    CameraAngles current = cam_tes->get_camera_angles();
                    CameraAngles target = cam_tes->get_total_object_to_cam_angles();
                    std::cout 
                            << "[Visible] Cam "<< cam_tes->getID() <<" Pos: (" << current.azimuth * 180 / M_PI << ", " << current.elevation * 180 / M_PI << ") | "
                            << "Target: (" << target.azimuth * 180 / M_PI << ", " << target.elevation * 180 / M_PI << ") | "
                            << "Diff: (" << (target.azimuth - current.azimuth) * 180 / M_PI << ") | (u,v) = (" 
                            << cam_tes->pix_info.u << ", " << cam_tes->pix_info.v << ")\n";
                }

                // for(auto cam_print : obj_hidden_cameras){
                //     Camera* cam_tes = cam_print; 
                //     CameraAngles current = cam_tes->get_camera_angles();
                //     CameraAngles target = cam_tes->get_total_object_to_cam_angles();
                //     std::cout 
                //             << "[Invisible] Cam "<< cam_tes->getID() <<" Pos: (" << current.azimuth * 180 / M_PI<< ", " << current.elevation * 180 / M_PI<< ") | "
                //             << "Target: (" << target.azimuth * 180 / M_PI << ", " << target.elevation * 180 / M_PI << ") | "
                //             << "Diff: (" <<( target.azimuth - current.azimuth) * 180 / M_PI<< ") | (u,v) = (" 
                //             << cam_tes->pix_info.u << ", " << cam_tes->pix_info.v << ")\n";
                // }

                if(one_cam_seen_obj_frame_count >= TRACKING_THRESHOLD_FRAMES) {
                    is_initialized = false;
                    state = REACQUIRE_STATE;
                    lost_frame_info.lost_time_info.clear();
                } 
                break;
            }


            /*============================================================ */
            /*===========================PF HERE========================= */
            /*============================================================ */
            case REACQUIRE_STATE:
            {
                reacquire_frame_count++;

                // std::cout << "REQUIRE STATE at frame [ " << current_frame << " ]\n";
                Eigen::MatrixXd P = ukf.getCovariance();

                P(0,0) *= 10;  P(3,3) *= 10; P(6,6) *= 10;
                P(1,1) *= 10;  P(4,4) *= 10; P(7,7) *= 10;

                if(!is_initialized){
                    pf.init_particles_from_ukf(ukf.getState(), P);
                    is_initialized = true;

                    /*=================================================================*/
                    /* THIS FOR DEBUG */
                    lost_frame_info.frame = current_frame;
                    /*=================================================================*/
                }

                /*predict and update based on gaussian state pf*/
                
                pf.predict(ALG_PROCESSING_TIMESTEP, ACC_DEVIATION);
                pf.update_weights_focus(obj_visible_cameras, obj_hidden_cameras);
                
                /*========= NEW KMEANS CLUSTERING ================*/
                Eigen::Vector3d focus_target;
                
                // ===== SWITCH MODE =====
                if (!obj_visible_cameras.empty())
                {
                    reacquire_mode = FOCUS_MODE;
                }
                else {
                    reacquire_mode = SEARCH_MODE;
                }
                
                GaussState rmse_state;
                rmse_state.mean = pf.get_estimated_state();

                double rmse_x = rmse_state.mean(0);
                double rmse_y = rmse_state.mean(3);
                double rmse_z = rmse_state.mean(6);
                Eigen::Vector3d rmse_pos(rmse_x,rmse_y,rmse_z);
                Eigen::Vector3d error = rmse_pos - current_position;
                Rmse_based_frame rmse_frame;
                rmse_frame.frame = current_frame;rmse_frame.rmse = error.norm();

                // ===== SEARCH MODE =====
                if (reacquire_mode == SEARCH_MODE)
                {
                    is_require_state_initialized = false;
                    /*============= weighted k means ==============*/
                    auto hypotheses = weighted_kmeans_return_centers(pf.get_particles(), CLUSTERING_NUMBER, 5, 1e-4);
                    for(int i = 0; i < obj_hidden_cameras.size(); i++){
                        int k = i % CLUSTERING_NUMBER;
                        Eigen::Vector3d target(
                            hypotheses[k](0),
                            hypotheses[k](3),
                            hypotheses[k](6)
                        );
                        auto* cam = obj_hidden_cameras[i];
                        cam->rotate_to_direction_vector(target - cam->getPosition());
                    }
                    // ===== GMM =====
                    // auto gmm = fit_gmm_em(pf.get_particles(), CLUSTERING_NUMBER, 5);

                    // // sort hypotheses based on weight
                    // std::vector<int> order(gmm.size());
                    // std::iota(order.begin(), order.end(), 0); /*create a list of hypotheses*/

                    // std::sort(order.begin(), order.end(),
                    //         [&](int a, int b){
                    //             return gmm[a].weight > gmm[b].weight;
                    //         });

                    // // assign mỗi camera 1 hypothesis khác nhau
                    // for (int i = 0; i < obj_hidden_cameras.size(); ++i)
                    // {
                    //     int k = order[i % gmm.size()];

                    //     Eigen::Vector3d target(
                    //         gmm[k].mean(0),
                    //         gmm[k].mean(3),
                    //         gmm[k].mean(6));

                    //     auto* cam = obj_hidden_cameras[i];
                    //     cam->rotate_to_direction_vector(target - cam->getPosition());
                    // }
                    temp_det_rmse_list.push_back(rmse_frame);
                }

                // ===== FOCUS MODE =====
                else if (reacquire_mode == FOCUS_MODE)
                {
                    /*================== FOR DEBUG =============================*/
                    if(!is_require_state_initialized){
                        LostTimePerFrameInfo lost_time_per_frame_info;
                        focus_time_per_frame_count++;
                        lost_time_per_frame_info.reacquire_count = focus_time_per_frame_count;
                        lost_time_per_frame_info.reacquire_cam_numbers = cameras_seen_obj;
                        lost_frame_info.lost_time_info.push_back(lost_time_per_frame_info);
                        is_require_state_initialized = true;
                        std::cout << "lost_time_per_frame_info is "<< lost_time_per_frame_info.reacquire_count << "\n";
                    }
                    /*==============================================================*/
                    pf.resample_if_needed();
                    Eigen::VectorXd state_pf = pf.get_estimated_state();

                    focus_target = Eigen::Vector3d(
                        state_pf(0),
                        state_pf(3),
                        state_pf(6));
                    
                    // tất cả camera quay về target này
                    rotate_cameras(obj_hidden_cameras, focus_target);
                    rotate_cameras(obj_visible_cameras, focus_target);

                    temp_det_rmse_list.push_back(rmse_frame);
                }
                
                /*============= FOR PLOTTING ===============*/
                Eigen::VectorXd x_cur = pf.get_estimated_state();
                double px = x_cur(0);
                double py = x_cur(3);
                double pz = x_cur(6);
                Eigen::Vector3d estimate_pos = Eigen::Vector3d(px, py, pz); 
                result_traj.push_back(estimate_pos);

                // std::cout << "Real position: (" << traj[current_frame].pos[0] << ", " 
                //             << traj[current_frame].pos[1] << ", "
                //             << traj[current_frame ].pos[2] << ") | " 
                //             << "rmse is " << rmse_frame.rmse << "\n";
                //  for(auto cam_print : obj_hidden_cameras){
        
                //     Camera* cam_tes = cam_print; 
                //     CameraAngles current = cam_tes->get_camera_angles();
                //     CameraAngles target = cam_tes->get_total_object_to_cam_angles();
                //     std::cout 
                //             << "Cam "<< cam_tes->getID() <<" Pos: (" << current.azimuth << ", " << current.elevation << ") | "
                //             << "Target: (" << target.azimuth * 180 / M_PI << ", " << target.elevation * 180 / M_PI << ") | "
                //             << "Diff: (" << (target.azimuth - current.azimuth) * 180 / M_PI << ") | (u,v) = (" 
                //             << cam_tes->pix_info.u << ", " << cam_tes->pix_info.v << ")\n";
                // }
                
                 /*export json for plotting FOV of cameras information at the beginning*/
                if(!is_json_pf){
                    const Point obj_p(current_position(0), current_position(1), current_position(2));
                    const Point est_p(x_cur(0),x_cur(3), x_cur(6));

                    std::vector<Polyhedron> fov_polys;
                    for(auto &cam: cameras){
                        fov_polys.push_back(cam.get_FOV()); 
                    }
                    std::string file_name = "pf_json/" + std::to_string(count_pf) + "_pf.json";
                    export_pf_json(fov_polys, pf.get_particles(), obj_p, est_p, file_name);
                    count_pf++;
                }

                /*Change state to tracking state*/
                if (obj_visible_cameras.size() >= 4)
                {   
                    reacquire_frame_count = 0;
                    state = TRACKING_STATE;
                    result.number_recover_from_reacquire++;
                    Eigen::MatrixXd pf_covariance = pf.compute_ellipsoid_cov();
                    std::cout << "it can run till here\n";
                    ukf.setState(x_cur);
                    ukf.setCovariance(pf_covariance);
                    is_json_pf = true;
                    lost_frame_info.isSuccess = true;
                    lost_frame_list_info.push_back(lost_frame_info);
                    focus_time_per_frame_count = 0;

                    for(auto& rmse_f : temp_det_rmse_list)
                        result.rmse_list_det_succ.push_back(rmse_f);
                    temp_det_rmse_list.clear();
                } 

                /*Change state to reset state*/
                if(reacquire_frame_count > REACQUIRE_THRESHOLD_FRAMES)
                {
                    state = RESET_STATE; 
                    is_json_pf = true; 
                    reacquire_frame_count = 0;
                    lost_frame_info.isSuccess = false;
                    lost_frame_list_info.push_back(lost_frame_info);
                    focus_time_per_frame_count = 0;

                    for(auto& rmse_f : temp_det_rmse_list)
                        result.rmse_list_det_fail.push_back(rmse_f);
                    temp_det_rmse_list.clear();
                }

                break;
            }
        }
            
        current_frame++;   
    }
    if (result.number_recover_from_reacquire == 0 && result.pf_penalty_count == 0){
        result.successful_redetection_percentage = 1;
    } else {
        result.successful_redetection_percentage = result.number_recover_from_reacquire / (result.number_recover_from_reacquire + result.pf_penalty_count);
    }
    export_trajectories_to_json(result_traj, "two_traj.json");
   
}

void Simulation::show_result(){
    std::cout << "=========================\n";
    std::cout << "current frame is " << current_frame << " in total " << total_frames <<" frames.\n";
    std::cout << "=========================\n";

    std::cout << "The missing frame in each camera is: \n";
    for(int i =0 ; i < cameras.size(); i++){
        int miss_frames = total_frames - cameras[i].obj_visible_frame_count;
        double miss_frames_percentage =  1.0 * miss_frames / (double) total_frames; 
        std::cout << "\t camera " << i+1 << " has missed " 
                << miss_frames << " frames."
                << "(" << miss_frames_percentage * 100.0 << "%)."<<'\n';
        std::cout << "\t min distance is " << cameras[i].min_distance 
                << " and max distance is " <<  cameras[i].max_distance <<'\n';
        std::cout << "====================================\n";
    } 
    
    std::cout << "\ntotal scanning time is " << result.pf_penalty_count << "\n";
    std::cout << "\nnumber recover from REACQUIRE STATE " << result.number_recover_from_reacquire << "\n";
    std::cout << "\n percentage recover from REACQUIRE STATE " << result.successful_redetection_percentage << "\n";

    std::cout << "==================================\n";
    for(auto curFrame : lost_frame_list){
        std::cout << "[Frame : " << curFrame.frame << "]  ";
        std::cout << "( " << curFrame.pos[0] << " "<< curFrame.pos[1] << " " << curFrame.pos[2] << " )\n";
        for (auto lost_obj: curFrame.list)
        {
            std::cout << "pixel info u is " << lost_obj.u << " and v is " << lost_obj.v << "\n";
        }
    }
    std::cout << "==================================\n";
    int reacquire_fail_with_no_focus_state_count = 0;
    int reacquire_fail_with_focus_state_count = 0;

    for(auto curFrame : lost_frame_list_info){
        std::string reacuire_state = curFrame.isSuccess ? "Success" : "False";
        std::cout << "[Frame : " << curFrame.frame << "] | " 
                    << "[Reacquie : " << reacuire_state << "] | \n";
        std::cout << "number of lost time in a frame is " << curFrame.lost_time_info.size() << "\n";  
        for (auto lost_obj: curFrame.lost_time_info)
        { 
            std::cout << "reacuire time " << lost_obj.reacquire_count 
            << " and number of camera can see is " << lost_obj.reacquire_cam_numbers << "\n";
        }

        if(!curFrame.isSuccess){
            if(curFrame.lost_time_info.size() > 0)
                reacquire_fail_with_focus_state_count++;
            else 
                reacquire_fail_with_no_focus_state_count++;
        }
    }
    std::cout << "\nreacquire_fail_with_focus_state_count is " << reacquire_fail_with_focus_state_count 
            << " |\n reacquire_fail_with_no_focus_state_count is " << reacquire_fail_with_no_focus_state_count
            << std::endl;
}