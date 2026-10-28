#include "read_gltf.hpp"

#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/math/math.hpp> // nearest_rotation
#include <ipc/utils/logger.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <igl/edges.h>
#include <tiny_gltf.h>

#include <cmath>
#include <numeric>
#include <unordered_set>

namespace ipc::rigid {

namespace {

    /// @brief The node's transform relative to its parent.
    ///
    /// glTF stores this either as a column-major 4×4 matrix or as separate
    /// translation, rotation, and scale (TRS), applied as T·R·S.
    Eigen::Affine3d local_transform(const tinygltf::Node& node)
    {
        Eigen::Affine3d transform = Eigen::Affine3d::Identity();
        if (node.matrix.size() == 16) {
            transform.matrix() =
                Eigen::Map<const Eigen::Matrix4d>(node.matrix.data());
            return transform;
        }
        if (node.translation.size() == 3) {
            transform.translate(
                Eigen::Vector3d(
                    node.translation[0], node.translation[1],
                    node.translation[2]));
        }
        if (node.rotation.size() == 4) {
            // NOTE: glTF uses (x, y, z, w)
            transform.rotate(
                Eigen::Quaterniond(
                    node.rotation[3], node.rotation[0], node.rotation[1],
                    node.rotation[2]));
        }
        if (node.scale.size() == 3) {
            transform.scale(
                Eigen::Vector3d(node.scale[0], node.scale[1], node.scale[2]));
        }
        return transform;
    }

    /// @brief Append one primitive's elements to (V, F, E).
    ///
    /// Triangles go to F and lines to E (points add vertices only). The
    /// primitive's indices are offset by the vertices already in V, so a
    /// multi-primitive mesh becomes one indexed mesh.
    ///
    /// @return False (with a warning) if the primitive cannot be read as a
    ///         dense float triangle, line, or point list; V, F, and E are then
    ///         left unchanged.
    bool extract_mesh_data(
        const tinygltf::Model& model,
        const tinygltf::Primitive& primitive,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        Eigen::MatrixXi& E)
    {
        // tinygltf reports an unspecified mode as -1 (glTF's default is
        // TRIANGLES).
        int simplex_size;
        switch (primitive.mode) {
        case -1:
        case TINYGLTF_MODE_TRIANGLES:
            simplex_size = 3;
            break;
        case TINYGLTF_MODE_LINE:
            simplex_size = 2;
            break;
        case TINYGLTF_MODE_POINTS:
            simplex_size = 1;
            break;
        default:
            logger().warn(
                "[GLTF] skipping primitive with mode {} (only triangles, "
                "lines, and points are supported)",
                primitive.mode);
            return false;
        }

        // --- 1. Vertices (POSITION) ---
        const auto position = primitive.attributes.find("POSITION");
        if (position == primitive.attributes.end()) {
            logger().warn("[GLTF] skipping primitive without POSITION");
            return false;
        }
        const tinygltf::Accessor& vAccessor = model.accessors[position->second];
        if (vAccessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT
            || vAccessor.type != TINYGLTF_TYPE_VEC3 || vAccessor.bufferView < 0
            || vAccessor.sparse.isSparse) {
            logger().warn(
                "[GLTF] skipping primitive whose POSITION is not a dense float "
                "VEC3 accessor");
            return false;
        }
        const tinygltf::BufferView& vBufferView =
            model.bufferViews[vAccessor.bufferView];
        const tinygltf::Buffer& vBuffer = model.buffers[vBufferView.buffer];

        // The stride is usually 3 floats, but interleaved data is larger.
        const int vByteStride = vAccessor.ByteStride(vBufferView);
        if (vByteStride <= 0) {
            logger().warn("[GLTF] skipping primitive with an invalid stride");
            return false;
        }
        const unsigned char* vData =
            &vBuffer.data[vBufferView.byteOffset + vAccessor.byteOffset];

        // --- 2. Elements (INDICES) ---
        // An unindexed primitive lists each element's vertices in order.
        const size_t num_indices = primitive.indices < 0
            ? vAccessor.count
            : model.accessors[primitive.indices].count;
        if (num_indices % simplex_size != 0) {
            logger().warn(
                "[GLTF] skipping primitive with {} indices (not a list of "
                "{}-vertex elements)",
                num_indices, simplex_size);
            return false;
        }

        std::vector<int> indices(num_indices);
        if (primitive.indices < 0) {
            std::iota(indices.begin(), indices.end(), 0);
        } else {
            const tinygltf::Accessor& iAccessor =
                model.accessors[primitive.indices];
            if (iAccessor.bufferView < 0) {
                logger().warn(
                    "[GLTF] skipping primitive whose indices have no buffer");
                return false;
            }
            const tinygltf::BufferView& iBufferView =
                model.bufferViews[iAccessor.bufferView];
            const unsigned char* iData =
                &model.buffers[iBufferView.buffer]
                     .data[iBufferView.byteOffset + iAccessor.byteOffset];
            for (size_t i = 0; i < num_indices; ++i) {
                switch (iAccessor.componentType) {
                case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
                    indices[i] =
                        int(reinterpret_cast<const uint32_t*>(iData)[i]);
                    break;
                case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
                    indices[i] =
                        int(reinterpret_cast<const uint16_t*>(iData)[i]);
                    break;
                case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
                    indices[i] = int(iData[i]);
                    break;
                default:
                    logger().warn(
                        "[GLTF] skipping primitive with index component type "
                        "{}",
                        iAccessor.componentType);
                    return false;
                }
            }
        }

        // --- 3. Append ---
        const Eigen::Index vOffset = V.rows();
        V.conservativeResize(vOffset + Eigen::Index(vAccessor.count), 3);
        for (size_t i = 0; i < vAccessor.count; ++i) {
            const auto* p =
                reinterpret_cast<const float*>(vData + i * size_t(vByteStride));
            V.row(vOffset + Eigen::Index(i)) << double(p[0]), double(p[1]),
                double(p[2]);
        }

        if (simplex_size > 1) {
            Eigen::MatrixXi& elements = simplex_size == 3 ? F : E;
            const Eigen::Index offset = elements.rows();
            const size_t num_elements = num_indices / simplex_size;
            elements.conservativeResize(
                offset + Eigen::Index(num_elements), simplex_size);
            for (size_t i = 0; i < num_elements; ++i) {
                for (int j = 0; j < simplex_size; ++j) {
                    elements(offset + Eigen::Index(i), j) =
                        int(vOffset) + indices[simplex_size * i + j];
                }
            }
        }
        return true;
    }

    void process_node(
        const tinygltf::Model& model,
        const tinygltf::Node& node,
        const Eigen::Affine3d& parent_transform,
        std::vector<Eigen::MatrixXd>& rest_positions,
        std::vector<Eigen::MatrixXi>& edges,
        std::vector<Eigen::MatrixXi>& faces,
        std::vector<double>& densities,
        std::vector<Pose>& initial_poses,
        std::vector<std::string>& names)
    {
        // A node's transform is relative to its parent, so we accumulate the
        // world transform down the hierarchy.
        const Eigen::Affine3d transform =
            parent_transform * local_transform(node);
        const Eigen::Matrix3d L = transform.linear();

        // A node scaled to zero on every axis (directly or by an ancestor) has
        // no geometry left, so we skip the body entirely (but still recurse
        // to its children).
        const bool zero_scale = L.cwiseAbs().maxCoeff() < 1e-12;

        if (node.mesh > -1 && zero_scale) {
            logger().info("Skipping zero-scale body: {}", node.name);
        } else if (node.mesh > -1) {
            const tinygltf::Mesh& mesh = model.meshes[node.mesh];
            logger().info("Found Body: {} (Mesh: {})", node.name, mesh.name);

            Eigen::MatrixXd V(0, 3);
            Eigen::MatrixXi F(0, 3), line_edges(0, 2);
            for (const auto& primitive : mesh.primitives) {
                extract_mesh_data(model, primitive, V, F, line_edges);
            }

            if (V.rows() == 0) {
                logger().warn(
                    "[GLTF] skipping body {} with no readable primitives",
                    node.name);
            } else {
                // We split the world linear map into a rotation and a
                // stretch, L = R·S (polar decomposition). The stretch (scale
                // and any shear from non-uniform parent scales) is baked into
                // the rest positions, and the rotation goes into the pose.
                const Eigen::Matrix3d R = nearest_rotation(L);
                const Eigen::Matrix3d S = R.transpose() * L;
                if (!S.isIdentity(1e-12)) {
                    V = V * S.transpose();
                }
                if (L.determinant() < 0) {
                    // A mirroring transform flips the triangles inside out,
                    // so we restore their winding.
                    F.col(1).swap(F.col(2));
                }

                Pose pose = Pose::identity(3);
                pose.position = transform.translation();
                pose.rotation = rotation_matrix_to_vector(R);
                logger().info(
                    "  - Position: {}, Rotation: {}", pose.position.transpose(),
                    pose.rotation.transpose());

                // The body's edges are its triangles' edges plus any line
                // primitives (codimensional edges).
                Eigen::MatrixXi E(0, 2);
                if (F.rows() > 0) { // igl::edges does not accept empty faces
                    igl::edges(F, E);
                }
                if (line_edges.rows() > 0) {
                    E.conservativeResize(E.rows() + line_edges.rows(), 2);
                    E.bottomRows(line_edges.rows()) = line_edges;
                }

                rest_positions.push_back(std::move(V));
                faces.push_back(std::move(F));
                edges.push_back(std::move(E));
                initial_poses.push_back(std::move(pose));

                // We only read geometry here. Physical properties (density,
                // body type, velocities, ...) are set by the caller after
                // loading, so every body starts with a default density.
                densities.emplace_back(1e3);

                // We prefer the node name and fall back to the mesh name. An
                // empty name lets RigidBodies generate a unique "Body{i}".
                names.emplace_back(!node.name.empty() ? node.name : mesh.name);
            }
        }

        // Recursively process children (important for nested objects/stacks)
        for (int childIdx : node.children) {
            process_node(
                model, model.nodes[childIdx], transform, rest_positions, edges,
                faces, densities, initial_poses, names);
        }
    }

    /// @brief The root nodes to load: the default scene's, or every node that
    /// is no other node's child when the file has no usable scene.
    std::vector<int> root_nodes(const tinygltf::Model& model)
    {
        const int scene = model.defaultScene > -1 ? model.defaultScene : 0;
        if (scene < int(model.scenes.size())) {
            logger().info(
                "Loading scene: {} with {} root nodes.",
                model.scenes[scene].name, model.scenes[scene].nodes.size());
            return model.scenes[scene].nodes;
        }

        logger().warn("[GLTF] no scene to load; using every root node");
        std::unordered_set<int> children;
        for (const tinygltf::Node& node : model.nodes) {
            children.insert(node.children.begin(), node.children.end());
        }
        std::vector<int> roots;
        for (int i = 0; i < int(model.nodes.size()); ++i) {
            if (children.count(i) == 0) {
                roots.push_back(i);
            }
        }
        return roots;
    }

} // namespace

std::pair<std::shared_ptr<RigidBodies>, std::vector<Pose>> read_gltf(
    const std::string& filename,
    const bool convert_planes,
    const std::unordered_map<std::string, VectorMax6b>& is_dof_fixed)
{
    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;

    // .glb is the binary format, .gltf is the ASCII/JSON format
    bool ret;
    if (filename.size() >= 4
        && filename.compare(filename.size() - 4, 4, ".glb") == 0) {
        ret = loader.LoadBinaryFromFile(&model, &err, &warn, filename);
    } else {
        ret = loader.LoadASCIIFromFile(&model, &err, &warn, filename);
    }

    if (!warn.empty()) {
        logger().warn("[GLTF] {}", warn);
    }

    if (!err.empty()) {
        logger().error("[GLTF] {}", err);
    }

    if (!ret) {
        logger().error("Failed to parse glTF file: {}", filename);
        return { nullptr, {} };
    }

    std::vector<Eigen::MatrixXd> rest_positions;
    std::vector<Eigen::MatrixXi> edges;
    std::vector<Eigen::MatrixXi> faces;
    std::vector<double> densities;
    std::vector<Pose> initial_poses;
    std::vector<std::string> names;

    for (int nodeIdx : root_nodes(model)) {
        process_node(
            model, model.nodes[nodeIdx], Eigen::Affine3d::Identity(),
            rest_positions, edges, faces, densities, initial_poses, names);
    }

    // We make the names unique here, over every mesh node in file order,
    // rather than leaving it to build_from_meshes, which only sees the bodies
    // left after plane pruning. That way a body's name depends only on the
    // file, so callers can predict it and use it as a key below.
    names = RigidBodies::make_unique_names(names);

    // Look up each body's DOF mask by name; an empty mask leaves every DOF
    // free. We warn about unmatched names since they are most likely typos.
    std::vector<VectorMax6b> body_is_dof_fixed(names.size());
    for (const auto& [name, mask] : is_dof_fixed) {
        bool found = false;
        for (size_t i = 0; i < names.size(); ++i) {
            if (names[i] == name) {
                body_is_dof_fixed[i] = mask;
                found = true;
            }
        }
        if (!found) {
            logger().warn(
                "is_dof_fixed: no body named '{}' in {}", name, filename);
        }
    }

    // build_from_meshes prunes plane bodies from initial_poses, names, and the
    // DOF masks in place, so they all stay aligned with the returned bodies.
    std::shared_ptr<RigidBodies> bodies = RigidBodies::build_from_meshes(
        rest_positions, edges, faces, densities, initial_poses, convert_planes,
        body_is_dof_fixed, names);
    return { bodies, initial_poses };
}

} // namespace ipc::rigid