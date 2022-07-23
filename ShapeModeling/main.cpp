#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

/*** insert any libigl headers here ***/
#include <igl/facet_components.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/boundary_loop.h>
#include <igl/loop.h>
#include <igl/edge_topology.h>
#include <igl/barycenter.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include "gui_utils.h"

#define PI 3.1415926f

enum MouseMode { NONE, FACE_SELECT, VERTEX_SELECT, TRANSLATE, ROTATE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;
// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd colors_per_face;

std::set<int> selected_faces;
int selected_v = -1;
void update_display(igl::opengl::glfw::Viewer& viewer);

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1')
    {
        //add your code for computing vertex to face relations here
        //store in VF,VFi
        VF.clear();
        VFi.clear();

        size_t vertNum = V.rows();
        VF.resize(vertNum);
        VFi.resize(vertNum);

        for (size_t fi = 0; fi < F.rows(); ++fi)
        {
            for (size_t i = 0; i < F.cols(); ++i)
            {
                VF[F(fi, i)].push_back(fi);
                VFi[F(fi, i)].push_back(i);
            }
        }

    }
    if (key == '2')
    {
        //add your code for computing vertex to vertex relations here
        //store in VV
        int maxCoeff = 0;

        for (size_t fi = 0; fi < F.rows(); ++fi)
        {
            for (size_t i = 0; i < F.cols(); ++i)
            {

                maxCoeff = std::max(F(fi, i), maxCoeff);
            }
        }

        VV.resize(maxCoeff + 1);
        for (size_t fi = 0; fi < F.rows(); ++fi)
        {
            for (size_t i = 0; i < F.cols(); ++i)
            {
                int s = F(fi, i);
                int d = F(fi, (i + 1) % F.cols());
                VV[s].push_back(d);
                VV[d].push_back(s);
            }
        }

        for (size_t Index = 0; Index < VV.size(); ++Index)
        {
            std::sort(VV[Index].begin(), VV[Index].end());
            VV[Index].erase(std::unique(VV[Index].begin(), VV[Index].end()), VV[Index].end());
        }
    }




    if (key == '3')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_face.setZero(F.rows(), 3);
        //add your code for computing per-face connected components here
        //store the component labels in cid
        igl::facet_components(F, cid);
        //compute colors for the faces based on components
        //store the colors in component_colors_per_face
        igl::jet(cid,true,colors_per_face);
        //set the viewer colors
        viewer.data().set_colors(colors_per_face);
    }

    if (key == '4') {
        Eigen::MatrixXd Vout = V;
        Eigen::MatrixXi Fout = F;
        // Add your code for sqrt(3) subdivision here.
        Eigen::MatrixXd BC;
        igl::barycenter(V, F, BC);
        Vout.setZero(V.rows()+F.rows(),V.cols());
        Vout << V, BC;
        Eigen::MatrixXi EV, FE, EF;
        igl::edge_topology(V, F, EV, FE, EF);
        Fout.setZero(F.rows()*3, F.cols());

        int c = 0;
        for (int i = 0; i < EV.rows(); i++)
        {
            //if is border
            if (EF(i, 0) == -1 || EF(i, 1) == -1) {
                if (EF(i, 0) == -1) {
                    int v0 = EV(i, 0);
                    int v1 = EV(i, 1);
                    int v = EF(i, 1) + V.rows();

                    Fout.row(c) = Eigen::Vector3i(v1, v0, v);
                }
                else {
                    int v0 = EV(i, 0);
                    int v1 = EV(i, 1);
                    int v = EF(i, 0) + V.rows();

                    Fout.row(c) = Eigen::Vector3i(v0, v1, v);
                }
                c++;
            }
            else {
                int face1 = EF(i, 0) + V.rows();
                int face2 = EF(i, 1) + V.rows();
                int v0 = EV(i, 0);
                int v1 = EV(i, 1);
                Fout.row(c) = Eigen::Vector3i(face2, face1, v0);
                Fout.row(c + 1) = Eigen::Vector3i(face1, face2, v1);
                c += 2;
            }
        }

        igl::adjacency_list(F, VV);
        Eigen::MatrixXd sumVec;
        sumVec.setZero(1, 3);
        for (int i = 0; i < V.rows(); i++) {
            int n = VV[i].size();
            sumVec.setZero();
            for (int j = 0; j < n; j++) {
                sumVec += V.row(VV[i][j]); //VV[i][j] is index of j-th vector adjacent to i-th vector
            }
            double a = (4 - 2 * cos((2 * PI) / n)) / 9;
            Vout.row(i) = (1 - a) * V.row(i) + a / n * sumVec;
        }
        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        update_display(viewer);
    }

    return true;
}

std::set<int> get_v_from_faces_idx(const Eigen::MatrixXi& F, std::set<int>& face_idx) {
    std::set<int> v_set;
    for (auto f : face_idx) {
        v_set.insert(F(f, 0)); v_set.insert(F(f, 1)); v_set.insert(F(f, 2));
    }
    return v_set;
}

void extrude(igl::opengl::glfw::Viewer& viewer) {
    Eigen::MatrixXd Vout = V;
    Eigen::MatrixXi Fout = F;

    // Get selected faces
    Eigen::MatrixXi sF(selected_faces.size(), 3); int idx = 0;
    for (auto it = selected_faces.begin(); it != selected_faces.end(); it++) { sF.row(idx++) = F.row(*it); }

    // Assert selected faces are connected
    Eigen::VectorXi comp; igl::facet_components(sF, comp);
    if (comp.maxCoeff() != 0) { std::cout << "Error: Not a single connected component, #face_comp =  " << comp << std::endl; return; }

    // 1) Get the boundary vertices surrounding the selected faces
    std::vector<int> bnd_loop;
    igl::boundary_loop(sF, bnd_loop);

    // 2) Duplicate boundary vertices
    Vout.resize(V.rows() + bnd_loop.size(), 3);
    for (int i = 0; i < V.rows(); i++) Vout.row(i) = V.row(i); // set vertices as old vertices
    for (int i = 0; i < bnd_loop.size(); i++) { Vout.row(V.rows() + i) = V.row(bnd_loop[i]); } // create new vertices as duplicates of the faces boundary

    // 3) Compute direction T: The average of selected face normals
    Eigen::RowVector3d T; T.setZero(); Eigen::MatrixXd FN;
    igl::per_face_normals(V, F, FN);
    for (auto it = selected_faces.begin(); it != selected_faces.end(); it++) { T += FN.row(*it); }
    T.normalize(); T *= 0.25 * (V.row(bnd_loop[1]) - V.row(bnd_loop[0])).norm();

    // 4) Offset old vertices by T
    std::set<int> inner_v = get_v_from_faces_idx(F, selected_faces);
    for (auto v : inner_v) {
        Vout.row(v) += T;
    }

    // 5) Update Fout 
    Fout.resize(F.rows() + 2 * bnd_loop.size(), 3); // 2 new faces per new edge (= per new vertex)
    for (int i = 0; i < F.rows(); i++) Fout.row(i) = F.row(i); // set first 'F.rows()' faces as the old faces

    // Add your code for updating Fout here

    // 5.1) Get the set of faces containing the old boundary vertices (hint: call igl::vertex_triangle_adjacency on the old 'F')
    std::set<int> Fset;
    std::vector<std::vector<int>> VFold, VFiold;
    igl::vertex_triangle_adjacency(V, F, VFold, VFiold);
    for (int i = 0; i < sF.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            std::vector<int> adjFaces = VFold[sF(i, j)];
            for (int k = 0; k < adjFaces.size(); k++)
            {
                Fset.insert(adjFaces[k]);
            }
        }
    }
    // 5.2) Get the "outer" set of faces containing the boundary vertices 
    //      (hint: call std::set_difference to compute the difference between the previously computed set of faces, and the selected faces)
    std::set<int> Fouter;
    std::set_difference(Fset.begin(), Fset.end(), selected_faces.begin(), selected_faces.end(), std::inserter(Fouter, Fouter.end()));
    // 5.3) Edit old outer faces indices, replacing the old vertices with the indices of the duplicated boundary vertices
    for (auto it = Fouter.begin(); it != Fouter.end(); it++)
    {
        int vertIndex = -1;
        for (int j = 0; j < Fout.cols(); j++)
        {
            vertIndex = Fout(*it, j);
            for (int k = 0; k < bnd_loop.size(); k++)
            {
                if (vertIndex == bnd_loop[k]) {
                    Fout(*it, j) = V.rows() + k;
                    break;
                }
            }
        }
    }
    // 5.4) Add new faces, 2 per edge
    int f_idx = F.rows();
    for (int i = 0; i < bnd_loop.size(); i++) {
        int v1, v2, v3, v4;
        // set v1,v2,v3,v4 correctly
        v1 = bnd_loop[(i + 1) % bnd_loop.size()];
        v2 = bnd_loop[i];
        v3 = V.rows() + i;
        v4 = V.rows() + (i + 1) % bnd_loop.size();

        Fout.row(f_idx++) << v1, v2, v3;
        Fout.row(f_idx++) << v3, v4, v1;
    }

    // 6) Check that the new mesh is a manifold (call is_edge_manifold, is_vertex_manifold on Vout,Fout)
    Eigen::MatrixXi B;
    if (!igl::is_edge_manifold(Fout) || !igl::is_vertex_manifold(Fout, B)) {
        std::cout << "Error: non-manifold!"<<std::endl;
    }
    
    // 7) Update V,F
    V = Vout; // uncomment for your code to take effect
    F = Fout; // uncomment for your code to take effect

    // Update gui and move to edit-translate mode
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(), 3); // number of faces has changed
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    for (auto f : selected_faces) { colors_per_face.row(f) << 1, 0, 0; }
    viewer.data().set_colors(colors_per_face);
    mouse_mode = TRANSLATE;
}

void clear_selection(igl::opengl::glfw::Viewer& viewer) {
    selected_faces.clear();
    selected_v = -1;
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(), 3);
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(colors_per_face);
}

void export_mesh() {
    std::string f = igl::file_dialog_save();
    igl::writeOFF(f, V, F);
}


bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
    down_mouse_x = viewer.current_mouse_x;
    down_mouse_y = viewer.current_mouse_y;

    if (mouse_mode == FACE_SELECT) {
        int f = pick_face(viewer, down_mouse_x, down_mouse_y, V, F);
        if (f != -1) {
            selected_faces.insert(f);
            selected_v = -1;
            // update face colors
            //colors_per_face.setConstant(colors_per_face.rows(), colors_per_face.cols(), 0);
            colors_per_face.row(f) << 1, 0, 0;
            viewer.data().set_colors(colors_per_face);
        }

    }
    else if (mouse_mode == VERTEX_SELECT) {
        int v = pick_vertex(viewer, down_mouse_x, down_mouse_y, V, F);
        if (v != -1) {
            selected_v = v;
            selected_faces.clear();
            update_display(viewer);
            viewer.data().set_points(V.row(selected_v), Eigen::RowVector3d(1, 0, 0));
        }
    }
    else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) {
        if (!selected_faces.empty()) {
            int f = pick_face(viewer, down_mouse_x, down_mouse_y, V, F);
            if (std::find(selected_faces.begin(), selected_faces.end(), f) != selected_faces.end()) {
                doit = true;
            }
        }
        else if (selected_v != -1) {
            int v = pick_vertex(viewer, down_mouse_x, down_mouse_y, V, F);
            if (v == selected_v) {
                viewer.data().set_points(V.row(selected_v), Eigen::RowVector3d(1, 0, 0));
                doit = true;
            }
        }
    }
    return false;
}

Eigen::RowVector3d get_face_avg(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::set<int>& selected_faces) {

    Eigen::RowVector3d avg; avg << 0, 0, 0;
    std::set<int> v_set = get_v_from_faces_idx(F, selected_faces);
    for (auto v : v_set) {
        avg += V.row(v);
    }
    avg /= v_set.size();

    return avg;
}

bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) {
    if (!doit)
        return false;
    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) {
        if (!selected_faces.empty()) {
            Eigen::RowVector3d face_avg_pt = get_face_avg(V, F, selected_faces);
            std::set<int> v_idx = get_v_from_faces_idx(F, selected_faces);
            if (mouse_mode == TRANSLATE) {
                Eigen::Vector3f translation = computeTranslation(viewer,
                    mouse_x,
                    down_mouse_x,
                    mouse_y,
                    down_mouse_y,
                    face_avg_pt);

                for (auto v_i : v_idx) { V.row(v_i) += translation.cast<double>(); }
            }
            else { // ROTATE
                Eigen::Vector4f rotation = computeRotation(viewer,
                    mouse_x,
                    down_mouse_x,
                    mouse_y,
                    down_mouse_y,
                    face_avg_pt);
                for (auto v_i : v_idx) {
                    Eigen::RowVector3f goalPosition = V.row(v_i).cast<float>();
                    goalPosition -= face_avg_pt.cast<float>();
                    igl::rotate_by_quat(goalPosition.data(), rotation.data(), goalPosition.data());
                    goalPosition += face_avg_pt.cast<float>();
                    V.row(v_i) = goalPosition.cast<double>();
                }
            }
            viewer.data().set_mesh(V, F);
            down_mouse_x = mouse_x;
            down_mouse_y = mouse_y;
            return true;
        }
        else if ((selected_v != -1) && (mouse_mode == TRANSLATE)) {
            Eigen::Vector3f translation = computeTranslation(viewer,
                mouse_x,
                down_mouse_x,
                mouse_y,
                down_mouse_y,
                V.row(selected_v));
            V.row(selected_v) += translation.cast<double>();
            viewer.data().set_mesh(V, F);
            viewer.data().set_points(V.row(selected_v), Eigen::RowVector3d(1, 0, 0));
            down_mouse_x = mouse_x;
            down_mouse_y = mouse_y;
            return true;
        }
    }
    return false;
}

bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
    doit = false;
    return false;
}

void update_display(igl::opengl::glfw::Viewer& viewer) {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(), 3);
    viewer.data().set_colors(colors_per_face);
}

int main(int argc, char* argv[]) {
    // Show the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_down = callback_key_down;

    if (argc == 2)
    {
        // Read mesh
        igl::readOFF(argv[1], V, F);

    }
    else
    {
        // Read mesh
        igl::readOFF("../data/honda.off", V, F);
    }

    viewer.data().set_mesh(V, F);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(V, F);
    
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // We can also use a std::vector<std::string> defined dynamically
            ImGui::Combo("Mouse Mode", (int*)&mouse_mode, "None\0Faces Select\0Vertex Select\0Translate\0Rotate");

            // Add a button
            if (ImGui::Button("Extrude", ImVec2(-1, 0)))
            {
                extrude(viewer);
            }
            // Add a button
            if (ImGui::Button("Clear selection", ImVec2(-1, 0)))
            {
                clear_selection(viewer);
            }
            // Add a button
            if (ImGui::Button("Export mesh", ImVec2(-1, 0)))
            {
                export_mesh();
            }
        }
    };

    viewer.callback_mouse_down = callback_mouse_down;
    viewer.callback_mouse_up = callback_mouse_up;
    viewer.callback_mouse_move = callback_mouse_move;
    update_display(viewer);
    viewer.launch();
}