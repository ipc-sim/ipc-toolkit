#include "read_gltf.hpp"

#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/utils/logger.hpp>

#include <Eigen/Core>
#include <igl/edges.h>
#include <tiny_gltf.h>

namespace ipc::rigid {

namespace {

    void extract_mesh_data(
        const tinygltf::Model& model,
        const tinygltf::Primitive& primitive,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F)
    {
        // --- 1. Extract Vertices (POSITION) ---
        // Look for the "POSITION" attribute in the mesh primitive
        const tinygltf::Accessor& vAccessor =
            model.accessors[primitive.attributes.at("POSITION")];
        const tinygltf::BufferView& vBufferView =
            model.bufferViews[vAccessor.bufferView];
        const tinygltf::Buffer& vBuffer = model.buffers[vBufferView.buffer];

        // Resize Eigen matrix: rows = num vertices, cols = 3 (x, y, z)
        const size_t vOffset = V.rows();
        V.conservativeResize(V.rows() + vAccessor.count, 3);

        // Pointer to the start of the position data
        const float* vData = reinterpret_cast<const float*>(
            &vBuffer.data[vBufferView.byteOffset + vAccessor.byteOffset]);

        // Calculate stride (usually 3 * sizeof(float), but check for
        // interleaved data)
        int vStride = vAccessor.ByteStride(vBufferView) / sizeof(float);

        for (size_t i = 0; i < vAccessor.count; ++i) {
            V(vOffset + i, 0) = static_cast<double>(vData[i * vStride + 0]);
            V(vOffset + i, 1) = static_cast<double>(vData[i * vStride + 1]);
            V(vOffset + i, 2) = static_cast<double>(vData[i * vStride + 2]);
        }

        // --- 2. Extract Faces (INDICES) ---
        const tinygltf::Accessor& iAccessor =
            model.accessors[primitive.indices];
        const tinygltf::BufferView& iBufferView =
            model.bufferViews[iAccessor.bufferView];
        const tinygltf::Buffer& iBuffer = model.buffers[iBufferView.buffer];

        // glTF supports 3 vertices per triangle. MatrixXi: rows = num faces,
        // cols = 3
        const size_t fOffset = F.rows();
        F.conservativeResize(F.rows() + iAccessor.count / 3, 3);

        const unsigned char* iDataPtr =
            &iBuffer.data[iBufferView.byteOffset + iAccessor.byteOffset];

        // Faces can be stored as scalar types: unsigned short, unsigned int, or
        // unsigned byte
        for (size_t i = 0; i < iAccessor.count / 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                size_t idx = i * 3 + j;
                if (iAccessor.componentType
                    == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
                    F(fOffset + i, j) = static_cast<int>(
                        reinterpret_cast<const uint32_t*>(iDataPtr)[idx]);
                } else if (
                    iAccessor.componentType
                    == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                    F(fOffset + i, j) = static_cast<int>(
                        reinterpret_cast<const uint16_t*>(iDataPtr)[idx]);
                } else if (
                    iAccessor.componentType
                    == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE) {
                    F(fOffset + i, j) = static_cast<int>(
                        reinterpret_cast<const uint8_t*>(iDataPtr)[idx]);
                }
            }
        }
    }

    void process_node(
        const tinygltf::Model& model,
        const tinygltf::Node& node,
        std::vector<Eigen::MatrixXd>& rest_positions,
        std::vector<Eigen::MatrixXi>& edges,
        std::vector<Eigen::MatrixXi>& faces,
        std::vector<double>& densities,
        std::vector<Pose>& initial_poses)
    {
        // 1. Check if this node is a visual 'body' (has a mesh)
        if (node.mesh > -1) {
            const tinygltf::Mesh& mesh = model.meshes[node.mesh];
            logger().info("Found Body: {} (Mesh: {})", node.name, mesh.name);
            rest_positions.emplace_back();
            faces.emplace_back();
            for (const auto& primitive : mesh.primitives) {
                extract_mesh_data(
                    model, primitive, rest_positions.back(), faces.back());
            }
            edges.emplace_back();
            igl::edges(faces.back(), edges.back());

            // 2. Extract Transform (Position, Rotation, Scale)
            // glTF uses column-major matrices or TRS (Translation, Rotation,
            // Scale)
            initial_poses.push_back(Pose::Identity(3)); // Default pose
            if (node.translation.size() == 3) {
                logger().info(
                    "  - Position: [{}, {}, {}]", node.translation[0],
                    node.translation[1], node.translation[2]);
                initial_poses.back().position << node.translation[0],
                    node.translation[1], node.translation[2];
            }

            if (node.scale.size() == 3) {
                logger().info(
                    "  - Scale: [{}, {}, {}]", node.scale[0], node.scale[1],
                    node.scale[2]);

                rest_positions.back() = rest_positions.back().array().rowwise()
                    * Eigen::Vector3d(
                          node.scale[0], node.scale[1], node.scale[2])
                          .transpose()
                          .array();
            }

            if (node.rotation.size() == 4) {
                logger().info(
                    "  - Rotation (quaternion): [{}, {}, {}, {}]",
                    node.rotation[0], node.rotation[1], node.rotation[2],
                    node.rotation[3]);
                // NOTE: glTF uses (x, y, z, w)
                Eigen::AngleAxisd aa(
                    Eigen::Quaterniond(
                        node.rotation[3], node.rotation[0], node.rotation[1],
                        node.rotation[2]));
                initial_poses.back().rotation = aa.angle() * aa.axis();
            }

            // 3. TODO: Add Physics Extension Loading Here
            // This is where you will check node.extensions for
            // "OMI_physics_body"
            densities.emplace_back(1e3); // Default density
        }

        // 4. Recursively process children (important for nested objects/stacks)
        for (int childIdx : node.children) {
            process_node(
                model, model.nodes[childIdx], rest_positions, edges, faces,
                densities, initial_poses);
        }
    }

} // namespace

std::pair<std::shared_ptr<RigidBodies>, std::vector<Pose>>
read_gltf(const std::string& filename, const bool convert_planes)
{
    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;

    // .glb is the binary format, .gltf is the ASCII/JSON format
    bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, filename);

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

    // Process the default scene
    int sceneToLoad = model.defaultScene > -1 ? model.defaultScene : 0;
    const tinygltf::Scene& scene = model.scenes[sceneToLoad];

    logger().info(
        "Loading scene: {} with {} root nodes.", scene.name,
        scene.nodes.size());

    std::vector<Eigen::MatrixXd> rest_positions;
    std::vector<Eigen::MatrixXi> edges;
    std::vector<Eigen::MatrixXi> faces;
    std::vector<double> densities;
    std::vector<Pose> initial_poses;

    for (int nodeIdx : scene.nodes) {
        process_node(
            model, model.nodes[nodeIdx], rest_positions, edges, faces,
            densities, initial_poses);
    }

    return { RigidBodies::build_from_meshes(
                 rest_positions, edges, faces, densities, initial_poses,
                 convert_planes),
             initial_poses };
}

} // namespace ipc::rigid