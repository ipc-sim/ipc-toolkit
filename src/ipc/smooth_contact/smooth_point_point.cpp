// #include "smooth_point_point.hpp"
// #include <queue>
// #include <ipc/utils/logger.hpp>
// #include <fstream>

// namespace ipc {

// /// @brief 
// /// @param v 
// /// @param ray normalized
// /// @param neighbors counter-clockwise
// /// @return 
// bool is_outside_object(
//     const Eigen::Ref<const RowVector3<double>>& v,
//     const Eigen::Ref<const RowVector3<double>>& ray,
//     const Eigen::Matrix<double, -1, 3> &neighbors)
// {
//     constexpr double tolerance = 1e-14;
//     assert(neighbors.rows() > 2);
    
//     enum class PrimitiveType {edge, face};
//     struct primitive {
//         int id;
//         PrimitiveType type;
//         double dist;

//         primitive(int id_, PrimitiveType type_, double dist_) : id(id_), type(type_), dist(dist_)
//         { }
//     };

//     primitive q(-1, PrimitiveType::face, std::numeric_limits<double>::max()), p(-1, PrimitiveType::face, std::numeric_limits<double>::max());
//     {
//         RowVector3<double> t_prev = (neighbors.row(neighbors.rows()-1) - v).normalized();
//         for (int f = 0; f < neighbors.rows(); f++)
//         {
//             RowVector3<double> t = (neighbors.row(f) - v).normalized();
            
//             // edges
//             {
//                 double tv = ray.dot(t);
//                 double nv = (ray - tv * t).norm();

//                 primitive cur(f, PrimitiveType::edge, nv);
//                 if (tv >= 0)
//                 {
//                     if (cur.dist < q.dist)
//                         q = cur;
//                 }
//                 else
//                 {
//                     if (cur.dist < p.dist)
//                         p = cur;
//                 }
//             }

//             RowVector3<double> normal = t_prev.cross(t).normalized();
//             {
//                 double a = ray.dot(normal.cross(t_prev).normalized());
//                 double b = ray.dot(normal.cross(t).normalized());
//                 double nv = abs(ray.dot(normal));
                
//                 primitive cur(f, PrimitiveType::face, nv);
//                 if (a >= 0 && b <= 0)
//                 {
//                     if (cur.dist < q.dist)
//                         q = cur;
//                 }
//                 else if (a <= 0 && b >= 0)
//                 {
//                     if (cur.dist < p.dist)
//                         p = cur;
//                 }
//             }
//             t_prev = t;
//         }
//     }

//     primitive top = (q.id >= 0) ? q : p;

//     const int f = top.id;
//     if (top.type == PrimitiveType::face)
//     {
//         RowVector3<double> t = (neighbors.row(f) - v).normalized();
//         RowVector3<double> t_prev = (neighbors.row((f + neighbors.rows() - 1) % neighbors.rows()) - v).normalized();
//         RowVector3<double> normal = t_prev.cross(t).normalized();
//         double val = ray.dot(normal);
//         if (abs(val) < tolerance)
//         {
//             logger().warn("normal computation may be incorrect! id {}, err {}", f, val);
//             // logger().warn("{}, {}, {}", v, ray, neighbors);
//             std::ofstream file("wrong_face" + std::to_string(abs(v(0) - neighbors(0,0))) + ".txt");
//             file << "t " << t << " tprev " << t_prev<< "\n";
//             file << "face " << f << " normal " << normal << " ray " << ray << "\n";
//             for (int i = 0; i < v.size(); i++)
//                 file << v(i) << ", ";
//             file << "\n";
//             for (int i = 0; i < ray.size(); i++)
//                 file << ray(i) << ", ";
//             file << "\n";
//             for (int i = 0; i < neighbors.rows(); i++)
//             {
//                 for (int j = 0; j < neighbors.cols(); j++)
//                     file << neighbors(i, j) << ", ";
//                 file << "\n";
//             }
//             file.close();
//             // exit(0);
//         }
//         return val > 0;
//     }
//     else // if (top.type == PrimitiveType::edge)
//     {
//         RowVector3<double> t0 = (neighbors.row((f + neighbors.rows() - 1) % neighbors.rows()) - v).normalized();
//         RowVector3<double> t1 = (neighbors.row(f) - v).normalized();
//         RowVector3<double> t2 = (neighbors.row((f + 1) % neighbors.rows()) - v).normalized();
//         RowVector3<double> n1 = t0.cross(t1).normalized(), n2 = t1.cross(t2).normalized();
//         double val = std::max(ray.dot(n1), ray.dot(n2));
//         if (abs(val) < tolerance)
//         {
//             logger().warn("normal computation may be incorrect! id {}, err {}", f, val);
//             // logger().warn("{}, {}, {}", v, ray, neighbors);
//             std::ofstream file("wrong_edge" + std::to_string(abs(v(0) - neighbors(0,0))) + ".txt");
//             for (int i = 0; i < v.size(); i++)
//                 file << v(i) << ", ";
//             file << "\n";
//             for (int i = 0; i < ray.size(); i++)
//                 file << ray(i) << ", ";
//             file << "\n";
//             for (int i = 0; i < neighbors.rows(); i++)
//             {
//                 for (int j = 0; j < neighbors.cols(); j++)
//                     file << neighbors(i, j) << ", ";
//                 file << "\n";
//             }
//             file.close();
//             // exit(0);
//         }
//         return val > 0;
//     }
// }

// bool smooth_point3_term_type(
//     const Eigen::Ref<const RowVector3<double>>& v,
//     const Eigen::Ref<const RowVector3<double>>& direc,
//     const Eigen::Matrix<double, -1, 3> &neighbors,
//     const double &alpha)
// {
//     RowVector3<double> t;
//     assert(neighbors.rows() > 2);

//     bool tangent_term = true;
//     for (int a = 0; a < neighbors.rows(); a++)
//     {
//         t = neighbors.row(a) - v;
//         tangent_term = tangent_term && (direc.dot(t) / t.norm() / alpha > -1);
//     }

//     if (!tangent_term)
//         return false;
    
//     if (!is_outside_object(v, -direc, neighbors))
//         return false;
    
//     return true;
// }

// }