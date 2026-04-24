#include "../../include/camera/camera.hpp"

/*This function is used for converting 
CGAL::Point3 type to Eigen::Vector3d type*/
Eigen::Vector3d toEigen(const Point& p)
{
    return Eigen::Vector3d(
        CGAL::to_double(p.x()),
        CGAL::to_double(p.y()),
        CGAL::to_double(p.z())
    );
}

inline double clamp(double v, double vmin, double vmax) {
    return std::max(vmin, std::min(vmax, v));
}

Polyhedron calculate_total_intersection(const std::vector<Camera*>& cam_list) {
    if (cam_list.empty()) return Polyhedron();
    std::vector<Polyhedron> polys;
    for(const Camera* cam: cam_list){
        polys.push_back(cam->get_FOV());
    }
    
    Polyhedron result = compute_intersection_LP(polys);

    return result;
}

Camera::Camera(
    int cam_id,
    double s_w, 
    double s_h,
    double l_az_max, 
    double r_az_max,
    double u_el_max, 
    double d_el_max,
    int res_w_max, 
    int res_h_max,
    int frames_p_s,
    double r_range_max,
    double f_length_min, 
    double f_length_max,
    double scale,
    double position_x,
    double position_y,
    double position_z
):
    id(cam_id),
    sensor_w(s_w), 
    sensor_h(s_h),
    left_azimuth_max(l_az_max), 
    right_azimuth_max(r_az_max),
    up_elevation_max(u_el_max), 
    down_elevation_max(d_el_max),
    res_w_max(res_w_max), 
    res_h_max(res_h_max),
    frames_per_sec(frames_p_s),
    recog_range_max(r_range_max),
    focal_length_min(f_length_min), 
    focal_length_max(f_length_max),
    catch_obj_scale(scale)
{
    focal_length = f_length_max;
    resolution_w = res_w_max;
    range = recog_range_max;
    position = Point(position_x,position_y,position_z);

    cam_angles.azimuth = 0.0;
    cam_angles.elevation = 0.0;
} 

Camera Camera::createFromFile(const std::string& filename){
    std::ifstream file(filename);
    std::map<std::string, double> conf;
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Không thể mở file config: " + filename);
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                conf[key] = std::stod(value);
            }
        }
    }

    return Camera(
        conf["cam_id"],
        conf["sensor_w"], conf["sensor_h"],
        conf["left_azimuth_max"], conf["right_azimuth_max"],
        conf["up_elevation_max"], conf["down_elevation_max"],
        static_cast<int>(conf["res_w_max"]), static_cast<int>(conf["res_h_max"]),
        static_cast<int>(conf["frames_per_sec"]),
        conf["recog_range_max"],
        conf["focal_length_min"], conf["focal_length_max"],
        conf["catch_obj_scale"],
        conf["pos_x"], conf["pos_y"], conf["pos_z"]
    );
}

void Camera::initial_fov() {
    double half_w = (sensor_w / (2.0 * focal_length)) * range;
    double half_h = (sensor_h / (2.0 * focal_length)) * range;

    std::vector<Point> pts = {
        Point(0, 0, 0),                               
        Point(half_w,  half_h, range),                 
        Point(-half_w,  half_h, range),                
        Point(-half_w, -half_h, range),                
        Point(half_w, -half_h, range)                  
    };

    local_fov_poly.clear();
    CGAL::convex_hull_3(pts.begin(), pts.end(), local_fov_poly);
}

CameraAngles Camera::get_total_object_to_cam_angles() {
 
    Eigen::Vector3d forward = get_camera_direction();
    Eigen::Vector3d world_up(0, 1, 0);
    Eigen::Vector3d right = world_up.cross(forward);
    if (right.squaredNorm() < 1e-9) right = Eigen::Vector3d::UnitX();
    right.normalize();
    Eigen::Vector3d up = forward.cross(right);
    up.normalize();

    double pixel_pitch = sensor_w / (double)res_w_max;
    double dx_m = (pix_info.u - res_w_max / 2.0) * pixel_pitch;
    double dy_m = (res_h_max / 2.0 - pix_info.v) * pixel_pitch; 

    /*direction vector in camera coordinate*/
    Eigen::Vector3d dir_cam(dx_m, dy_m, focal_length);
    dir_cam.normalize();

    /*direction vector in world coordinate*/
    Eigen::Vector3d dir_world = dir_cam.x() * right + 
                               dir_cam.y() * up + 
                               dir_cam.z() * forward;

    double total_az = std::atan2(dir_world.x(), dir_world.z()) ;
    
    double total_el = std::asin(dir_world.y() / dir_world.norm());

    total_obj_to_cam_angles.azimuth = total_az;
    total_obj_to_cam_angles.elevation = total_el;

    return total_obj_to_cam_angles;
}

Eigen::Vector3d Camera::getPosition() const{
    return toEigen(position);
}

Polyhedron Camera::get_FOV() const {

    double az_rad = cam_angles.azimuth ;
    double el_rad = cam_angles.elevation;

    double vx = std::sin(az_rad) * std::cos(el_rad); 
    double vy = std::sin(el_rad);               
    double vz = std::cos(az_rad) * std::cos(el_rad); 

    Eigen::Vector3d fwd(vx,vy,vz);
    Eigen::Vector3d world_up(0, 1, 0);
    Eigen::Vector3d rgt = world_up.cross(fwd);
    Eigen::Vector3d up = fwd.cross(rgt);
    fwd.normalize();
    rgt.normalize();
    up.normalize();

    CGAL::Aff_transformation_3<Kernel> M(
        rgt.x(), up.x(), fwd.x(), position.x(),
        rgt.y(), up.y(), fwd.y(), position.y(),
        rgt.z(), up.z(), fwd.z(), position.z()
    );

    Polyhedron transformed_fov = local_fov_poly;
    std::transform(transformed_fov.points_begin(), transformed_fov.points_end(),
                   transformed_fov.points_begin(), M);

    return transformed_fov;
}

CameraAngles Camera::calculate_angles_to(Point target) {
    Eigen::Vector3d dir = toEigen(target) - toEigen(position);
    CameraAngles ang;
    ang.azimuth = std::atan2(dir.x(), dir.z())  ;
    double horizontal_dist = std::sqrt(dir.x()*dir.x() + dir.z()*dir.z());
    ang.elevation = std::atan2(dir.y(), horizontal_dist);
    return ang;
}

void Camera::rotate_to_direction_vector(Eigen::Vector3d dir_vect)
{
    // ===== target angles =====
    double target_az = std::atan2(dir_vect.x(), dir_vect.z());
    double horizontal_dist = std::sqrt(dir_vect.x()*dir_vect.x() + dir_vect.z()*dir_vect.z());
    double target_el = std::atan2(dir_vect.y(), horizontal_dist);

    // ===== current angles =====
    double current_az = cam_angles.azimuth;
    double current_el = cam_angles.elevation;

    // ===== max delta (4 độ → rad) =====
    const double MAX_DELTA = 4.0 * M_PI / 180.0;

    // ===== delta =====
    double delta_az = target_az - current_az;
    double delta_el = target_el - current_el;

    // ===== normalize góc về [-pi, pi] (rất quan trọng) =====
    auto normalize_angle = [](double a) {
        while (a > M_PI) a -= 2 * M_PI;
        while (a < -M_PI) a += 2 * M_PI;
        return a;
    };

    delta_az = normalize_angle(delta_az);
    delta_el = normalize_angle(delta_el);

    // ===== clamp =====
    delta_az = std::clamp(delta_az, -MAX_DELTA, MAX_DELTA);
    delta_el = std::clamp(delta_el, -MAX_DELTA, MAX_DELTA);

    // ===== update =====
    cam_angles.azimuth += delta_az;
    cam_angles.elevation += delta_el;
}

Eigen::Vector3d Camera::get_camera_direction(){
    double az_rad = cam_angles.azimuth ;
    double el_rad = cam_angles.elevation;

    double vx = std::sin(az_rad) * std::cos(el_rad); 
    double vy = std::sin(el_rad);               
    double vz = std::cos(az_rad) * std::cos(el_rad); 

    Eigen::Vector3d forward(vx, vy, vz);
    return forward;
}

/*=================================================================*/
/*================ For rotating cameras ====================*/

void rotate_cameras(
    std::vector<Camera*> &cameras, 
    const Eigen::Vector3d target_pos)
{
    for(Camera* camera: cameras){
        Eigen::Vector3d dir = target_pos - camera->getPosition();
        camera->rotate_to_direction_vector(dir);
    }
}

/*========================================================*/
/*=================  SIMULATION ONLY =====================*/
/*========================================================*/

/*input : object coordinate , object radius
  output: object's pixel info*/
void Camera::calculate_pixel_projection(
    const Eigen::Vector3d& obj_world, double r
) {

    Eigen::Vector3d forward = get_camera_direction();
    Eigen::Vector3d world_up(0, 1, 0);
    
    Eigen::Vector3d right = world_up.cross(forward);
    if (right.squaredNorm() < 1e-9) right = Eigen::Vector3d::UnitX(); 
    Eigen::Vector3d up = forward.cross(right);

    forward.normalize();
    right.normalize();
    up.normalize();

    Eigen::Vector3d cam_pos = toEigen(position);
    Eigen::Vector3d rel_world = obj_world - cam_pos;

    double x_c = rel_world.dot(right);  
    double y_c = rel_world.dot(up);      
    double z_c = rel_world.dot(forward); 

    if (z_c <= 1e-3) { 
        obj_is_visible = false;
        pix_info.u = -1;
        pix_info.v = -1;
        return;
    }

    double x_m = focal_length * (x_c / z_c);
    double y_m = focal_length * (y_c / z_c);

    double pixel_pitch = sensor_w / (double)res_w_max; 
    
    double dist_to_center = rel_world.norm(); 
    double diameter_m = focal_length * (2.0 * r / dist_to_center);
    pix_info.size_px = static_cast<int>(std::round(diameter_m / pixel_pitch));
    
    pix_info.u = (x_m / pixel_pitch) + (res_w_max / 2.0);
    pix_info.v = (res_h_max / 2.0) - (y_m / pixel_pitch);

    // 8. FOV Check
    if (pix_info.u >= 0 && pix_info.u < res_w_max && 
        pix_info.v >= 0 && pix_info.v < res_h_max) {
        obj_is_visible = true;
    } else {
        obj_is_visible = false;
    }
}

