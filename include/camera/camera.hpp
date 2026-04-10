#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <map>

#include "convex_polyhedron.hpp"

struct Pixel_Obj {  /* This struct for simulation only*/
    double u,v;
    int size_px;
};

struct CameraAngles {
    double azimuth; 
    double elevation;
};

class Camera {
    private:
        const int id;
        const double sensor_w, sensor_h;
        const double left_azimuth_max, right_azimuth_max;
        const double up_elevation_max, down_elevation_max;
        const int res_w_max, res_h_max;
        const int frames_per_sec;
        const double recog_range_max;
        const double focal_length_min, focal_length_max;
        const double catch_obj_scale;

    /*camera changing properties*/
        double focal_length; /*we fixed this in simulation and = focal_length_max */
        int resolution_w; /*pixels, we fixed this in simulation and its value equals res_w_max*/
        double range; /*(meters), we fixed this in simulation and its value equals recog_range_max*/

        /*angle formed by camera principal axis and initial setup axis*/
        CameraAngles cam_angles;
        
        /*angles formed by line between camera screen and object and camera principal axis 
        + cam_angles */
        CameraAngles total_obj_to_cam_angles;
        
        Polyhedron field_of_view;

        Polyhedron local_fov_poly; 
        // DKHierarchy cached_local_dk;
        
        /*Positions can change or not*/
        Point position; 
        
    public:
        bool obj_is_visible;
        Pixel_Obj pix_info;
        int obj_visible_frame_count = 0;
        double min_distance=1000;
        double max_distance=0;
        const double MAX_ROTATION_ANGLE = 4; /*deg*/
        
        /*Constructor*/
        Camera(
            int id,
            double s_w_max, 
            double s_h_max,
            double l_az_max, 
            double r_az_max,
            double l_el_max, 
            double r_el_max,
            int p_w_max, 
            int p_h_max,
            int frames_p_s,
            double r_range_max,
            double f_length_min, 
            double f_length_max,
            double scale,
            double position_x,
            double position_y,
            double position_z
        );

        static Camera createFromFile(const std::string& filename);

        /*Getters*/
        double getSensorW(){return sensor_w;}
        double getSensorH(){return sensor_h;} 
        double getLeftAzimuthMax(){return left_azimuth_max;} 
        double getRightAzimuthMax(){return right_azimuth_max;} 
        double getUpElevatorMax(){return up_elevation_max;}
        double getDownElevatorMax(){return down_elevation_max;} 
        int getResWmax(){return res_w_max;}
        int getResHmax(){return res_h_max;} 
        int getFramesPerSec(){return frames_per_sec;}
        double getRecogRangeMax(){return recog_range_max;}
        double getFocalLengthMin(){return focal_length_min;} 
        double getFocalLengthMax(){return focal_length_max;}
        double getCatchObjScale(){return catch_obj_scale;}

        double getFocalLength(){return focal_length;}
        double getResolutionWidth(){return resolution_w;}
        int getID(){return id;}

        Polyhedron get_FOV() const;
        
        Eigen::Vector3d get_camera_direction();

        void set_camera_angles_initialized(double az, double el){
            cam_angles.azimuth = az * M_PI / 180;
            cam_angles.elevation = el * M_PI / 180;
        }

        CameraAngles get_camera_angles() const{return cam_angles;}
        void set_camera_angles(double az, double el){
            cam_angles.azimuth = az;
            cam_angles.elevation = el;
        }

        void calculate_pixel_projection(
            const Eigen::Vector3d& obj_world,
            double r
        );

        CameraAngles get_total_object_to_cam_angles();

        Eigen::Vector3d getPosition() const;

        CameraAngles get_half_fov_angles() const;

        // DKHierarchy get_dk_hierarchy() const;
        void initial_fov();

        CameraAngles calculate_angles_to(Point target_p);
        void rotate_to_direction_vector(Eigen::Vector3d dir_vect);
        
};


Polyhedron calculate_total_intersection(const std::vector<Camera*>& cam_list);

/*This function is used for converting 
CGAL::Point3 type to Eigen::Vector3d type*/
Eigen::Vector3d toEigen(const Point& p);

