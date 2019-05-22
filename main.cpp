#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  float ref_vel = 0;  //starting velocity
  h.onMessage([&ref_vel,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;

    int target_lane; // target lane for trajectory generaton
    float max_vel = 22.2;  //50 MPH converted to mps

    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

            car_speed = car_speed*0.44704;  //speed converted to mps unit

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

            int prev_size = previous_path_x.size();
  			       if (prev_size > 0)
  			          {
  				              car_s = end_path_s;
  			          }
            //determine current lane
            int current_lane =  div(car_d,4).quot;

            // declare flags indicating weather space surrounding the vehilce is blocked
            bool front_blocked = false;
            bool left_blocked = false;
            bool right_blocked = false;
            float close_thr_dist = max_vel*0.02*prev_size;

           // process each car from the simulator one by one in the for loop
           for (int i =0; i<sensor_fusion.size(); i++)
               {    float other_car_d = sensor_fusion[i][6];
                    int other_car_lane =  div(other_car_d,4).quot;

                    float vx = sensor_fusion[i][3];
                    float vy = sensor_fusion[i][4];
                    float other_car_speed = sqrt(vx*vx + vy*vy);

                    //predict position of this car
                    float other_car_s_end = sensor_fusion[i][5];
                    other_car_s_end += (other_car_speed*0.02*prev_size);
                    // how close is it to our car
                    float dist_btwn_cars = other_car_s_end - end_path_s;
                    // check if this car is blocking lanes surrounding our car
                     if ((other_car_lane == current_lane) && (dist_btwn_cars < close_thr_dist) && (other_car_s_end > end_path_s))
                        {front_blocked = true;}

                     if (current_lane == 0)
                         {left_blocked = true;}
                     else if ((other_car_lane == current_lane - 1) && (dist_btwn_cars < (close_thr_dist + 10)) && (dist_btwn_cars > -10))
                        {left_blocked = true;}

                     if (current_lane == 2)
                         {right_blocked = true;}
                     else if ((other_car_lane == current_lane + 1) && (dist_btwn_cars < (close_thr_dist + 10)) && (dist_btwn_cars > -10))
                         {right_blocked = true;}
                  }

          //decide on the maneouver, which in turns decides the target lane and reference velocity
          if (front_blocked)   //check if front space in the current lane is blocked
             { if (!left_blocked)  //if left lane is open, change to left lane
                  {if(ref_vel > max_vel-2)
                       {ref_vel -= 0.19;}
                   target_lane = current_lane - 1;}
               else if (!right_blocked) // if right lane is open, change to right lane
                  { if(ref_vel > max_vel-2)
                       {ref_vel -= 0.19;}
                   target_lane = current_lane + 1;}
               else    // if none of the lane open, slow down
                  {ref_vel = fmax(0,(ref_vel - 0.19));
                   target_lane = current_lane;}
             }
          else {   //front space is clear, so drive at max allowable speed
              ref_vel = fmin(max_vel, (ref_vel + 0.19));
              target_lane = current_lane;
               }

            /* Attemp to use quintic polynomial trajectory but abanded
            double s_t_f = 1;  // time span for trajectory generation
            double s_f = car_s + car_speed*straj_time_final;
            double s_dot_f = 22.35;
            double s_doubledot_f = 100;

            double s_f = car_d;
            double d_dot_f = 0;
            double d_doubledot_f = 0;

            Eigen::VectorXd start_d,start_s, end_d, end_s;
            start_s = {car_s, car_speed*cos(car_yaw),s_doubledot_f};
            end_s = {s_f, s_dot_f, s_doubledot_f};

            start_d = {car_d, car_speed*sin(car_yaw),d_doubledot_f};
            end_d = {d_f, d_dot_f, d_doubledot_f};   */


            //Create two vectors to save spline points in (x,y) coordinates
      			vector<double> ptsx;
      			vector<double> ptsy;

            //declare two vectors to send (x,y) points to simulator
            vector<double> next_x_vals;
            vector<double> next_y_vals;

      			//reference states from the previous trajectory
      			double ref_x, ref_y, ref_x_prev, ref_y_prev, ref_yaw;

      			//Looks at the spare points of previous trajectory
      			if (prev_size < 2)
            {
              // if car is starting (previous trajectory is empty, back-calculate car position
      				ref_x_prev = car_x - cos(car_yaw);
      				ref_y_prev = car_y - sin(car_yaw);

              ref_x = car_x;
        			ref_y = car_y;

              ref_yaw = deg2rad(car_yaw);
      			}
      			else
            { //read the last two previous trajectory points
      				ref_x = previous_path_x[prev_size-1];
      				ref_y = previous_path_y[prev_size-1];

      				ref_x_prev = previous_path_x[prev_size-2];
      				ref_y_prev = previous_path_y[prev_size-2];

      				ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);
      			}
            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);

            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);

      			//In Frenet add 3 evely spaced points at 30m distances
      			vector<double> next_wp0 = getXY(car_s+30, (2 + 4*target_lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
      			vector<double> next_wp1 = getXY(car_s+60, (2 + 4*target_lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
      			vector<double> next_wp2 = getXY(car_s+90, (2 + 4*target_lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);

            // add three pointst to spline array
      			ptsx.push_back(next_wp0[0]);
      			ptsx.push_back(next_wp1[0]);
      			ptsx.push_back(next_wp2[0]);

      			ptsy.push_back(next_wp0[1]);
      			ptsy.push_back(next_wp1[1]);
      			ptsy.push_back(next_wp2[1]);

      			//reseting the points to reference values before crating spline
      			for (int i = 0; i< ptsx.size();i++){
      				double shift_x = ptsx[i] - ref_x;
      				double shift_y = ptsy[i] - ref_y;
      				ptsx[i] = (shift_x * cos(0-ref_yaw)-shift_y * sin(0-ref_yaw));
      				ptsy[i] = (shift_x * sin(0-ref_yaw)+shift_y * cos(0-ref_yaw));
      			}

      			//create a spline
      			tk::spline s;

      			//set (x,y) points to the spline
      			s.set_points(ptsx,ptsy);

            //Now interpolate the spline to correct number of points for simulator

      			//start the trajectory with unsed previous points
      			for (int i = 0; i<previous_path_x.size();i++){
      				next_x_vals.push_back(previous_path_x[i]);
      				next_y_vals.push_back(previous_path_y[i]);
      			}

      			//fill the rest of trajectory with interpolated spline points
             //choose some target distance based on number of points in trajectory
             // 1 second long trajectory is chosen i.e. 50 points
      			double target_x = 25.0;  // 25 meter distance longer than car can travel in 1 s.
      			double target_y = s(target_x); //evaluating spline at the end of target
      			double target_dist = sqrt(target_x*target_x + target_y*target_y);

      			double x_add_on = 0;
      			//fill up the rest of trajectory
      			for (int i = 1; i<= 50-previous_path_x.size();i++){
      				double N = (target_dist/(0.02*ref_vel));  // number of divisions/slices of target distance
      				double x_point = x_add_on +  target_x/N;  // add 1 slice to X
      				double y_point = s(x_point); // caluclate correspoinding y

      				x_add_on = x_point;

      				double x_ref = x_point;
      				double y_ref = y_point;

      				//rotate the spline back to simulator(map) coordinates
      				x_point = (x_ref * cos(ref_yaw)-y_ref * sin(ref_yaw));
      				y_point = (x_ref * sin(ref_yaw)+y_ref * cos(ref_yaw));
              //shift the spline back to simulator(map) coordinates
      				x_point += ref_x;
      				y_point += ref_y;

      				next_x_vals.push_back(x_point);
      				next_y_vals.push_back(y_point);
      			}

            /*
            double dist_inc = 0.44;
            for (int i = 0; i < 50; i++)
            { double next_s = car_s + (i+1)*dist_inc;
              double next_d = car_d;
              vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
              next_x_vals.push_back(xy[0]);
              next_y_vals.push_back(xy[1]);
            }    */

            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
