#include "write_gltf.hpp"

#include <tiny_gltf.h>

#include <numeric>

namespace ipc::rigid {

bool write_gltf(
    const std::string& filename,
    const RigidBodies& bodies,
    const std::list<std::vector<Pose>>& poses,
    double timestep,
    bool embed_buffers,
    bool write_binary,
    bool prettyPrint)
{
    typedef float Float;
    int float_component_type;
    if (std::is_same<Float, float>::value) {
        float_component_type = TINYGLTF_COMPONENT_TYPE_FLOAT;
    } else {
        // assert(std::is_same<Float, double>::value);
        float_component_type = TINYGLTF_COMPONENT_TYPE_DOUBLE;
    }

    using namespace tinygltf;

    assert(bodies.dim() == 3);
    size_t num_bodies = bodies.num_bodies();
    assert(poses.size() > 0);
    size_t num_steps = poses.size();

    Model model;
    model.defaultScene = 0;

    model.asset.version = "2.0";
    model.asset.generator = "RigidIPC";

    Scene& scene = model.scenes.emplace_back();
    scene.name = "RigidIPCSimulation";
    scene.nodes.resize(num_bodies);
    std::iota(scene.nodes.begin(), scene.nodes.end(), 0);

    model.nodes.resize(num_bodies);
    model.animations.resize(1);
    Animation& animation = model.animations[0];
    animation.name = "Simulation";
    animation.channels.resize(2 * num_bodies);
    animation.samplers.resize(2 * num_bodies);
    model.meshes.resize(num_bodies);
    // Four accessors per bidy: vertices, faces, translations, and rotations
    model.accessors.resize(4 * num_bodies + 1);
    {
        Accessor& accessor = model.accessors[2 * num_bodies];
        accessor.name = "Times";
        accessor.bufferView = 2 * num_bodies;
        accessor.componentType = float_component_type;
        accessor.count = num_steps;
        accessor.minValues.push_back(timestep);
        accessor.maxValues.push_back(num_steps * timestep);
        accessor.type = TINYGLTF_TYPE_SCALAR;
    }

    for (int i = 0; i < num_bodies; i++) {
        Node& node = model.nodes[i];
        node.mesh = i;
        std::string body_name = "Body" + std::to_string(i);
        node.name = body_name;
        node.translation = { { poses.front()[i].position.x(),
                               poses.front()[i].position.y(),
                               poses.front()[i].position.z() } };
        // NOTE: Use identity for the first frame by applying the initial
        // rotation to the vertices and removing the rotation from further
        // rotations.
        // Eigen::Quaternion<double> q = poses[0][i].construct_quaternion();
        Eigen::Quaternion<double> q = Eigen::Quaternion<double>::Identity();
        node.rotation = { { q.x(), q.y(), q.z(), q.w() } };

        animation.channels[2 * i + 0].sampler = 2 * i;
        animation.channels[2 * i + 0].target_node = i;
        animation.channels[2 * i + 0].target_path = "translation";
        animation.channels[2 * i + 1].sampler = 2 * i + 1;
        animation.channels[2 * i + 1].target_node = i;
        animation.channels[2 * i + 1].target_path = "rotation";

        animation.samplers[2 * i + 0].input = 2 * num_bodies;
        animation.samplers[2 * i + 0].output = 2 * num_bodies + 2 * i + 1;
        animation.samplers[2 * i + 0].interpolation = "LINEAR";
        animation.samplers[2 * i + 1].input = 2 * num_bodies;
        animation.samplers[2 * i + 1].output = 2 * num_bodies + 2 * i + 2;
        animation.samplers[2 * i + 1].interpolation = "LINEAR";

        Mesh& mesh = model.meshes[i];
        mesh.name = body_name;
        Primitive& primitive = mesh.primitives.emplace_back();
        primitive.attributes["POSITION"] = 2 * i;
        primitive.indices = 2 * i + 1;
        primitive.mode = TINYGLTF_MODE_TRIANGLES;

        Accessor* accessor = &model.accessors[2 * i];
        accessor->name = body_name + "Vertices";
        accessor->bufferView = 2 * i;
        accessor->componentType = float_component_type;
        accessor->count = bodies.num_body_vertices(i);
        // accessor->max = ...;
        // accessor->min = ...;
        accessor->type = TINYGLTF_TYPE_VEC3;

        accessor = &model.accessors[2 * i + 1];
        accessor->name = body_name + "Faces";
        accessor->bufferView = 2 * i + 1;
        accessor->componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
        accessor->count = 3 * bodies.num_body_faces(i);
        accessor->type = TINYGLTF_TYPE_SCALAR;

        accessor = &model.accessors[2 * num_bodies + 2 * i + 1];
        accessor->name = body_name + "Translations";
        accessor->bufferView = 2 * num_bodies + 2 * i + 1;
        accessor->componentType = float_component_type;
        accessor->count = num_steps;
        accessor->type = TINYGLTF_TYPE_VEC3;

        accessor = &model.accessors[2 * num_bodies + 2 * i + 2];
        accessor->name = body_name + "Rotations";
        accessor->bufferView = 2 * num_bodies + 2 * i + 2;
        accessor->componentType = float_component_type;
        accessor->count = num_steps;
        accessor->type = TINYGLTF_TYPE_VEC4;
    }

    ///////////////////////////////////////////////////////////////////////////

    model.bufferViews.resize(4 * num_bodies + 1);
    size_t byte_offset = 0;
    for (int i = 0; i < num_bodies; i++) {
        std::string body_name = "Body" + std::to_string(i);

        BufferView* buffer_view = &model.bufferViews[2 * i];
        buffer_view->name = body_name + "Vertices";
        buffer_view->buffer = 0;
        buffer_view->byteLength =
            sizeof(Float) * bodies.body_rest_positions(i).size();
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;

        buffer_view = &model.bufferViews[2 * i + 1];
        buffer_view->name = body_name + "Faces";
        buffer_view->buffer = 0;
        buffer_view->byteLength =
            sizeof(unsigned int) * bodies.body_faces(i).size();
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;
    }
    {
        BufferView& buffer_view = model.bufferViews[2 * num_bodies];
        buffer_view.name = "Times";
        buffer_view.buffer = 0;
        buffer_view.byteLength = num_steps * sizeof(Float);
        buffer_view.byteOffset = byte_offset;
        byte_offset += buffer_view.byteLength;
    }
    for (int i = 0; i < num_bodies; i++) {
        std::string body_name = "Body" + std::to_string(i);

        BufferView* buffer_view =
            &model.bufferViews[2 * num_bodies + 2 * i + 1];
        buffer_view->name = body_name + "Translations";
        buffer_view->buffer = 0;
        buffer_view->byteLength = 3 * sizeof(Float) * num_steps;
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;

        buffer_view = &model.bufferViews[2 * num_bodies + 2 * i + 2];
        buffer_view->name = body_name + "Rotations";
        buffer_view->buffer = 0;
        buffer_view->byteLength = 4 * sizeof(Float) * num_steps;
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;
    }

    ///////////////////////////////////////////////////////////////////////////

    std::vector<unsigned char> byte_data(byte_offset);
    size_t byte_i = 0;
    for (int i = 0; i < num_bodies; i++) {
        // NOTE: Apply the initial rotation to the vertices
        Eigen::MatrixXd V =
            bodies.body_rest_positions(i) * bodies[i].R0().transpose();
        for (int r = 0; r < V.rows(); r++) {
            for (int c = 0; c < V.cols(); c++) {
                Float v = V(r, c);
                std::memcpy(&byte_data[byte_i], &v, sizeof(Float));
                byte_i += sizeof(Float);
            }
        }

        Eigen::MatrixXi F = bodies.body_faces(i);
        for (int r = 0; r < F.rows(); r++) {
            for (int c = 0; c < F.cols(); c++) {
                unsigned int fij = F(r, c);
                std::memcpy(&byte_data[byte_i], &fij, sizeof(unsigned int));
                byte_i += sizeof(unsigned int);
            }
        }
    }
    for (int i = 0; i < num_steps; i++) {
        Float t = i * timestep;
        std::memcpy(&byte_data[byte_i], &t, sizeof(Float));
        byte_i += sizeof(Float);
    }
    for (int i = 0; i < num_bodies; i++) {
        for (const auto& poses_j : poses) {
            Eigen::Vector3d p = poses_j[i].position;
            for (int d = 0; d < p.size(); d++) {
                Float pd = p[d];
                std::memcpy(&byte_data[byte_i], &pd, sizeof(Float));
                byte_i += sizeof(Float);
            }
        }

        const Eigen::Quaternion<double> q0(Eigen::Matrix3d(bodies[i].R0()));
        for (const auto& poses_j : poses) {
            // NOTE: Apply the inverse of the initial rotation to the
            // quaternion to get the rotation relative to the input orientation.
            Eigen::Quaternion<double> quat =
                poses_j[i].quaternion() * q0.inverse();
            Eigen::Vector4d q(quat.x(), quat.y(), quat.z(), quat.w());
            for (int d = 0; d < q.size(); d++) {
                Float qd = q[d];
                std::memcpy(&byte_data[byte_i], &qd, sizeof(Float));
                byte_i += sizeof(Float);
            }
        }
    }
    assert(byte_i == byte_data.size());

    model.buffers.emplace_back().data = byte_data;
    return TinyGLTF().WriteGltfSceneToFile(
        &model, filename,
        /*embedImages=*/true, embed_buffers, prettyPrint, write_binary);
}

} // namespace ipc::rigid
