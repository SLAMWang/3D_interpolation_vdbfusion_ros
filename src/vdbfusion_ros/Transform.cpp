// MIT License
//
// # Copyright (c) 2022 Saurabh Gupta, Ignacio Vizzo, Cyrill Stachniss, University of Bonn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "Transform.hpp"

#include <geometry_msgs/Transform.h>
#include <geometry_msgs/TransformStamped.h>
#include <ros/ros.h>
#include <tf2_ros/transform_listener.h>

#include <deque>

#include "sophus/se3.hpp"

using geometry_msgs::Transform;
using geometry_msgs::TransformStamped;

void invertTransform(TransformStamped& tf) {
    tf.transform = SE3ToTransform(TransformToSE3(tf.transform).inverse());
}

Transform interpolate(const Transform& tf_old, const Transform& tf_new, const double alpha) {
    auto T_old = TransformToSE3(tf_old);
    auto T_new = TransformToSE3(tf_new);
    return SE3ToTransform(T_old * Sophus::SE3d::exp(alpha * ((T_old.inverse() * T_new).log())));
}

vdbfusion::Transform::Transform(ros::NodeHandle& nh) : buffer_(ros::Duration(50, 0)), tf_(buffer_) {
    ROS_INFO("Transform init");
    nh.getParam("/use_tf_transforms", use_tf2_);
    if (use_tf2_) {
        nh.getParam("/parent_frame", parent_frame_);
        nh.getParam("/child_frame", child_frame_);
    } else {
        std::string tf_topic;
        nh.getParam("/tf_topic", tf_topic);
        const int queue_size = 500;
        tf_sub_ = nh.subscribe(tf_topic, queue_size, &vdbfusion::Transform::tfCallback, this);

        float tx, ty, tz, x, y, z, w;
        nh.getParam("/tx", tx);
        nh.getParam("/ty", ty);
        nh.getParam("/tz", tz);
        nh.getParam("/x", x);
        nh.getParam("/y", y);
        nh.getParam("/z", z);
        nh.getParam("/w", w);
        static_tf_.transform.translation.x = tx;
        static_tf_.transform.translation.y = ty;
        static_tf_.transform.translation.z = tz;
        static_tf_.transform.rotation.x = x;
        static_tf_.transform.rotation.y = y;
        static_tf_.transform.rotation.z = z;
        static_tf_.transform.rotation.w = w;

        bool invert_static_tf;
        nh.getParam("/invert_static_tf", invert_static_tf);
        if (invert_static_tf) {
            invertTransform(static_tf_);
        }
    }
}

void vdbfusion::Transform::tfCallback(const TransformStamped& transform_msg) {
    tf_queue_.push_back(transform_msg);
}

bool vdbfusion::Transform::lookUpTransform(const ros::Time& timestamp,
                                           const ros::Duration& tolerance,
                                           TransformStamped& transform) {
    if (use_tf2_) {
        return lookUpTransformTF2(parent_frame_, child_frame_, timestamp, tolerance, transform);
    } else {
        return lookUpTransformQ(timestamp, tolerance, transform);
    }
}

bool vdbfusion::Transform::lookUpTransformTF2(const std::string& parent_frame,
                                              const std::string& child_frame,
                                              const ros::Time& timestamp,
                                              const ros::Duration& tolerance,
                                              TransformStamped& transform) {
    if (buffer_.canTransform(parent_frame_, child_frame_, timestamp, tolerance)) {
        transform = buffer_.lookupTransform(parent_frame_, child_frame_, timestamp, tolerance);
        return true;
    }
    return false;
}

bool vdbfusion::Transform::lookUpTransformQ(const ros::Time& timestamp,
                                            const ros::Duration& tolerance,
                                            TransformStamped& transform) {
    if (tf_queue_.empty()) {
        ROS_WARN_STREAM_THROTTLE(30, "No match found for transform timestamp: "
                                         << timestamp << " as transform queue is empty.");
        return false;
    }
    bool match_found = false;
    std::deque<geometry_msgs::TransformStamped>::iterator it = tf_queue_.begin();
    for (; it != tf_queue_.end(); ++it) {
        if (it->header.stamp > timestamp) {
            if ((it->header.stamp - timestamp).toNSec() < tolerance.toNSec()) {
                match_found = true;
            }
            break;
        }

        if ((timestamp - it->header.stamp).toNSec() < tolerance.toNSec()) {
            match_found = true;
            break;
        }
    }

    if (match_found) {
        transform = *it;
    } else {
        if (it == tf_queue_.begin() || it == tf_queue_.end()) {
            ROS_WARN_STREAM_THROTTLE(30, "No match found for transform timestamp: " << timestamp);
            return false;
        }
        TransformStamped tf_newest = *it;
        int64_t newest_timestamp_ns = it->header.stamp.toNSec();
        --it;
        TransformStamped tf_oldest = *it;
        int64_t oldest_timestamp_ns = it->header.stamp.toNSec();

        double alpha = 0.0;
        if (newest_timestamp_ns != oldest_timestamp_ns) {
            alpha = static_cast<double>(timestamp.toNSec() - oldest_timestamp_ns) /
                    static_cast<double>(newest_timestamp_ns - oldest_timestamp_ns);
        }
        transform.transform = interpolate(tf_oldest.transform, tf_newest.transform, alpha);
    }

    transform.transform = SE3ToTransform(TransformToSE3(transform.transform) *
                                         (TransformToSE3(static_tf_.transform).inverse()));

    tf_queue_.erase(tf_queue_.begin(), it);
    return true;
}
