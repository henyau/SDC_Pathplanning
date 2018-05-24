#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
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
double ms2mph(double x) {return x*2.23694;}
double mph2ms(double x) {return x*0.44704;}

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

//returns the length of a subsection of a spline, if x_i< start control point or if x_f>end cp,
//spline library uses linear extrapolation.
double computeSplineLength(tk::spline s, double x_i, double x_f)
{//to accurately compute length should use a Gaussian quadrature
 //but first try a simple chorded approx to see 
    double accum = 0.0;
    double x_len = x_f-x_i;
    int intervals = 25;
    double x_inter = x_len/intervals;

    for(int i = 0; i<intervals; i++)
    {
        accum += distance(x_i, s(x_i), x_i+x_inter,s(x_i+x_inter)); 
        x_i += x_inter;    
    }


    return accum;
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
  //simulator refresh time
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
  double ref_speed = 0.0;             
  int lane = 1;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &ref_speed, &lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
       
            double del_t = 0.02;//time interval where "move" is executed
            double lane_width = 4.0;

            // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
    	    vector<double> next_y_vals;

          	//define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            //the vehicle simulator lacks dynamics and appears to just teleport to prescribed locations along path
            const double max_speed = 49.0; //mph,  conversion from mi/hr to m/s is done later 

            //create a simple state machine to slow down or change lanes if needed

            /*
             use zones to decide what(if any) action should be taken for any specific vehicle spotted
            |a|b|c|
            |d|0|e|
            |f| |g|*/

            bool a_s,b_s,c_s,d_s,e_s,f_s,g_s;// these states are relative to the car, so if in left lane then adf are true, meaning blocked
            a_s = b_s = c_s = d_s = e_s = f_s = g_s = false;
            if(lane == 2)
                c_s = e_s = g_s = true;
            else if(lane == 0)
                a_s = d_s = f_s = true;

            double gap = 30;
            //keep track of how fast the vehicles ahead are travelling
            double a_speed = max_speed;
            double b_speed = max_speed;
            double c_speed = max_speed;

            double b_dist = 1000.0;
            double a_dist = 1000.0;
            double c_dist = 1000.0;
            for (int i = 0; i < sensor_fusion.size(); i++) 
            {//id, xpos,ypos, xvel, yvel, s,d
             //in MAP coordinates
                double sense_d = sensor_fusion[i][6];//sensor_v[6];
                int sense_lane;
                if ( sense_d>=0 && sense_d<4)
                    sense_lane = 0;
                else if (sense_d>=4 && sense_d<8)
                    sense_lane = 1;
                else if (sense_d>=8 && sense_d<=12)
                    sense_lane = 2;
                else//ignore oncoming traffic
                    continue;

                double sense_s = sensor_fusion[i][5];//sensor_v[5];

                //adjusted buffers to be more conservative
                if(sense_lane-lane == 0) // same lane               
                {
                    if(sense_s > car_s-2.0 && (sense_s-car_s)<gap)
                    {    
                        b_s = true;
                        b_speed = distance(0,0,sensor_fusion[i][3],sensor_fusion[i][4]);
                        b_dist = fmin(sense_s-car_s,b_dist);
                    }
                }
                else if((sense_lane-lane) == 1)//right lane
                {
                    if(sense_s > car_s+2.0 && sense_s<car_s+2.0*gap)
                    {
                        c_s = true;
                        c_speed = distance(0,0,sensor_fusion[i][3],sensor_fusion[i][4]);
                        c_dist = fmin(sense_s-car_s, c_dist);
                    }
                    else if(sense_s > car_s-4.0 && sense_s <= car_s+4.0)
                        e_s = true;
                    else if(sense_s <= car_s-4.0 && sense_s >= car_s-15.0 )                       
                        g_s = true;
                }
                else if((sense_lane-lane) == -1)//left lane
                {
                    if(sense_s > car_s+2.0 && sense_s<car_s+2.0*gap)
                    {
                        a_s = true;
                        a_speed = distance(0,0,sensor_fusion[i][3],sensor_fusion[i][4]);
                        a_dist = fmin(sense_s-car_s, a_dist);
                   }
                    else if(sense_s > car_s-4.0 && sense_s <=car_s +4.0)
                        d_s = true;
                    else if(sense_s <= car_s-4.0 && sense_s >= car_s-15.0 )
                        f_s = true;                
                    
                }     
            }

            //check ahead to see if we need to slow down or change lanes
            if(b_s)
            {
                if(!d_s && !f_s && !a_s)//prefer to pass on left if possible
                    lane--;
                else if(!e_s && !g_s && !c_s)
                    lane++;
                //a and c are still blocked but as least they are faster than b
                else if(!d_s && !f_s && a_speed>b_speed && a_dist-b_dist>0.0)
                    lane--;
                else if(!e_s && !g_s && c_speed>b_speed && c_dist-b_dist>0.0)
                    lane++;
                //b is still faster, but slow down to avoid collision
                else if(ref_speed>ms2mph(b_speed))
                    ref_speed -= fmin((ref_speed-ms2mph(b_speed)),0.1);
                if (b_dist<6.0) // if too close slow down a little to get a gap
                    ref_speed -= 0.05;
            
            }
            else
            {
                if(ref_speed<max_speed)
                    ref_speed+=0.25;
            }

            //list of control points for spline
		    vector<double> ptsx;
    		vector<double> ptsy;

	    	double ref_x = car_x;
    		double ref_y = car_y;
	    	double ref_yaw = deg2rad(car_yaw);

		    //two conditions, either starting out or we've been going for awhile
            int prev_size = previous_path_x.size();
            //not enough points, linearly extrapolate a previous position
	    	if(prev_size <2)
    		{
                ptsx.push_back(car_x-cos(car_yaw));
                ptsy.push_back(car_y-sin(car_yaw));
                
                ptsx.push_back(car_x);
                ptsy.push_back(car_y);              
           }
            //use the previous end points to continue 
            else
            {
                ref_x = previous_path_x[prev_size-1];
                ref_y = previous_path_y[prev_size-1];
  
                double ref_x_prev = previous_path_x[prev_size-2];
                double ref_y_prev = previous_path_y[prev_size-2];
                ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev); //store end yaw
                
                
                ptsx.push_back(previous_path_x[prev_size-2]);
                ptsy.push_back(previous_path_y[prev_size-2]);

                ptsx.push_back(previous_path_x[prev_size-1]);
                ptsy.push_back(previous_path_y[prev_size-1]);              
        
            }

            //use Frenet coordinates to create next control points of spline
            //use regular s intervals
            
            double spline_cp_int = 50;
            if(b_s)//likely following car, going slow, need to shorten intervals
                spline_cp_int = 40;


            for(int i = 1; i<=3; i++)//order of spline
            {
                vector<double> cp = getXY(car_s+spline_cp_int*i, (0.5*lane_width+lane_width*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
                ptsx.push_back(cp[0]);
                ptsy.push_back(cp[1]);
            }

            //want control points to be have zero origin and oriented to last path position
            for(int i = 0; i<ptsx.size(); i++)
            {
                double shift_x = ptsx[i]-ref_x;
                double shift_y = ptsy[i]-ref_y;

                ptsx[i] = shift_x*cos(-ref_yaw)-shift_y*sin(-ref_yaw);
                ptsy[i] = shift_x*sin(-ref_yaw)+shift_y*cos(-ref_yaw);
            }
            //create interpolating spline
            tk::spline traj_spline;
            traj_spline.set_points(ptsx, ptsy);

            //copy previous path, then append new points
            for(int i = 0; i<prev_size; i++)
            {
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
            }
          
            //the spline just gives points in space without consideration to speed
            //we need to choose points such that the reference speed is obtained

            //first compute the length to some setpoint ahead
            double target_x = 30.0;
            double target_dist = computeSplineLength(traj_spline, 0, target_x);

            //each point in the path read and executed ever 0.02sec
            //use 50 points
            double x_accum = 0.0; //starting at x = 0, increment by 
            
            //ref_speed = 49.0; // give a little buffer for 50 mph

            //50 is hardcoded somewhere in simulator, breaks an assert if >50
            for(int i = 1; i<=50-previous_path_x.size(); i++)
            {
                double x_point = x_accum + del_t*ref_speed*0.44704; //converts mi/hr to m/s
                double y_point = traj_spline(x_point);

                x_accum = x_point;

                double x_ref = x_point;
                double y_ref = y_point;

                x_point = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw));
                y_point = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw));

                x_point += ref_x;//spline is 0 origin, shift to last point in path
                y_point += ref_y;

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }
            
           	
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
