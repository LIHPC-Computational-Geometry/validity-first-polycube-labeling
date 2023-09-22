#include <geogram/mesh/mesh.h>          // for MeshElementsFlags::MESH_FACETS
#include <geogram_gfx/mesh/mesh_gfx.h>  // for set_picking_mode()
#include <geogram_gfx/GLUP/GLUP.h>      // for glupUnProject()
#include <geogram/basic/attributes.h>   // for Attribute

#include "LabelingViewerApp.h"

#include <fmt/core.h>

#include <deque>

#include "containers.h" // for VECTOR_CONTAINS()

// values for paint_mode_
#define PAINT_MODE_DISABLED     0
#define PAINT_MODE_PENCIL       1
#define PAINT_MODE_BUCKET_FILL  2
#define PAINT_MODE_PIPETTE      3

using namespace GEO;

const char* label_selection_[] = {
    "0 = +X",
    "1 = -X",
    "2 = +Y",
    "3 = -Y",
    "4 = +Z",
    "5 = -Z"
};

class LabelingPainterApp : public LabelingViewerApp {
public:

    LabelingPainterApp() : LabelingViewerApp("labeling_painter") {
        paint_mode_ = PAINT_MODE_DISABLED;
        selected_label_ = 0;
        picked_facet_id_as_pixel_[0] = 0xff;
        picked_facet_id_as_pixel_[1] = 0xff;
        picked_facet_id_as_pixel_[2] = 0xff;
        picked_facet_id_as_pixel_[3] = 0xff;
        picked_facet_id_ = index_t(-1);
    }

private:

    void state_transition(State new_state) override {
        // if user change state -> disable paint mode
        paint_mode_ = PAINT_MODE_DISABLED;
        LabelingViewerApp::state_transition(new_state);
    }

    void labeling_visu_mode_transition(int new_mode) override {
        // if user change visu mode -> disable paint mode
        paint_mode_ = PAINT_MODE_DISABLED;
        update_static_labeling_graph(allow_boundaries_between_opposite_labels_); // recompute charts, boundaries and corners from the labeling facet attribute
        LabelingViewerApp::labeling_visu_mode_transition(new_mode);
    }

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(state_ == LabelingViewerApp::State::labeling) {
			ImGui::SeparatorText("Paint mode");
            switch (paint_mode_) {
                case PAINT_MODE_DISABLED:
                    ImGui::TextUnformatted(fmt::format("Current : {} disabled (camera pan)",icon_UTF8("hand-pointer")).c_str());
                    break;
                case PAINT_MODE_PENCIL:
                    ImGui::TextUnformatted(fmt::format("Current : {} pencil",icon_UTF8("pencil-alt")).c_str());
                    break;
                case PAINT_MODE_BUCKET_FILL:
                    ImGui::TextUnformatted(fmt::format("Current : {} bucket fill",icon_UTF8("fill-drip")).c_str());
                    break;
                case PAINT_MODE_PIPETTE:
                    ImGui::TextUnformatted(fmt::format("Current : {} pipette",icon_UTF8("eye-dropper")).c_str());
                    break;
                default:
                    geo_assert_not_reached;
            }
            if(ImGui::Button(icon_UTF8("hand-pointer"))) {
                paint_mode_ = PAINT_MODE_DISABLED;
            }
            ImGui::SameLine();
            if(ImGui::Button(icon_UTF8("pencil-alt"))) {
                if(labeling_visu_mode_ != VIEW_RAW_LABELING) {
                    labeling_visu_mode_transition(VIEW_RAW_LABELING); // avoid this function if the labeling visu mode is already "raw labeling", because there is the re-computation of the labeling graph inside
                }
                paint_mode_ = PAINT_MODE_PENCIL;
            }
            ImGui::SameLine();
            if(ImGui::Button(icon_UTF8("fill-drip"))) {
                if(labeling_visu_mode_ != VIEW_RAW_LABELING) {
                    labeling_visu_mode_transition(VIEW_RAW_LABELING); // avoid this function if the labeling visu mode is already "raw labeling", because there is the re-computation of the labeling graph inside
                }
                paint_mode_ = PAINT_MODE_BUCKET_FILL;
            }
            ImGui::SameLine();
            if(ImGui::Button(icon_UTF8("eye-dropper"))) {
                if(labeling_visu_mode_ != VIEW_RAW_LABELING) {
                    labeling_visu_mode_transition(VIEW_RAW_LABELING); // avoid this function if the labeling visu mode is already "raw labeling", because there is the re-computation of the labeling graph inside
                }
                paint_mode_ = PAINT_MODE_PIPETTE;
            }
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 4);
            ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) selected_label_)); // change the text color
            if (ImGui::BeginCombo("Label", label_selection_[selected_label_]))
            {
                // based on ext/geogram/src/lib/geogram_gfx/third_party/imgui/imgui_demo.cpp "combo 1" code snippet
                for (index_t n = 0; n < 6; n++) // for each label in [0,5]
                {
                    const bool is_selected = (selected_label_ == n);
                    ImGui::PushStyleColor(ImGuiCol_Text, labeling_colors_.color_as_ImVec4( (std::size_t) n)); // change the text color
                    if (ImGui::Selectable(label_selection_[n], is_selected))
                        selected_label_ = n;
                    ImGui::PopStyleColor();

                    if (is_selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            ImGui::PopStyleColor();
		}
	}

    void cursor_pos_callback( double x, double y, int source ) override {
        // based on ext/geogram/src/lib/geogram_gfx/gui/simple_application.cpp cursor_pos_callback()
        geo_argused(source);
        cursor_pos_ = vec2(x,y); // Added. Store x and y. Unit = pixels from top-left corner
        mouse_xy_ = vec2(
            double(x) / double(get_width()),
            double(y) / double(get_height())
        );
        mouse_xy_ *= 2.0;
        mouse_xy_ -= vec2(1.0, 1.0);
        if(three_D_ && mouse_op_ == MOUSE_ROTATE) {
            switch(mouse_target_) {
            case MOUSE_NOTARGET:
                break;
            case MOUSE_OBJECT:
                object_rotation_.drag(mouse_xy_);
                if(fixed_clip_) {
                clip_rotation_.drag(mouse_xy_);			
                }
                break;
            case MOUSE_LIGHT:
                light_rotation_.drag(mouse_xy_);		    
                break;
            case MOUSE_CLIP:
                clip_rotation_.drag(mouse_xy_);		    
                break;
            }
        } else if(mouse_op_ == MOUSE_ZOOM && mouse_target_ == MOUSE_OBJECT) {
            double R = mouse_xy_.y - mouse_down_xy_.y;
            double fact = (1.0 + R);
            fact = std::min(fact, 2.0);
            fact = std::max(fact, 0.1);
            zoom_ = zoom_down_ * fact;
        } else if( mouse_op_ == MOUSE_TRANSLATE) {
            double dx = mouse_xy_.x - mouse_down_xy_.x;	    	    
            double dy = mouse_xy_.y - mouse_down_xy_.y;
            if(paint_mode_ != PAINT_MODE_DISABLED) { // Added
                if((dx > 0.05) || (dy > 0.05)) { // If significative mouse movement
                    fmt::println(Logger::out("painting"),"Translation disabled in paint mode"); Logger::out("painting").flush();
                }
                return;
            }
            if(mouse_target_ == MOUSE_OBJECT) {
            object_translation_.x += 2.0 * dx / zoom_;
            object_translation_.y -= 2.0 * dy / zoom_;
            if(fixed_clip_) {
                clip_translation_.x += 2.0 * dx / zoom_;
                clip_translation_.y -= 2.0 * dy / zoom_;
            }
            } else if(mouse_target_ == MOUSE_CLIP) {
            clip_translation_.x += 2.0 * dx / zoom_;
            clip_translation_.y -= 2.0 * dy / zoom_;
            }
            mouse_down_xy_ = mouse_xy_;
        }
    }

    void mouse_button_callback(int button, int action, int mods, int source) override {
        SimpleApplication::mouse_button_callback(button, action, mods, source);
        if((paint_mode_ != PAINT_MODE_DISABLED) && (action==EVENT_ACTION_DOWN) && (button == 0)) { // if left click while in paint mode

            // see https://github.com/BrunoLevy/geogram/discussions/88
            // based on https://github.com/BrunoLevy/GraphiteThree/blob/main/src/lib/OGF/renderer/context/rendering_context.cpp get_picked_point()

            picked_facet_id_ = pick(MESH_FACETS);

            if(picked_facet_id_ == index_t(-1)) { // if no facet under the cursor
                return;
            }

            if(picked_facet_id_ >= mesh_.facets.nb()) { // if not in the range of facet indices
                fmt::println(Logger::err("painting"),"picked_facet_id_ = {} is greater of equal than number of facets = {}",picked_facet_id_,mesh_.facets.nb()); Logger::err("painting").flush();
                return;
            }

            Attribute<index_t> label(mesh_.facets.attributes(),LABELING_ATTRIBUTE_NAME); // fetch labeling

            if(paint_mode_ == PAINT_MODE_PIPETTE) {
                selected_label_ = label[picked_facet_id_];
                return;
            }

            index_t previous_label = label[picked_facet_id_];
            label[picked_facet_id_] = selected_label_; // change label of picked facet
            fmt::println(Logger::out("painting"),"Label of facet {} changed from {} to {}",picked_facet_id_,previous_label,selected_label_); Logger::out("painting").flush();

            if(paint_mode_ != PAINT_MODE_BUCKET_FILL) {
                return; // if pencil mode, stop here
            }

            if(previous_label == selected_label_) {
                return; // useless bucket fill
            }

            // bucket fill mode

            std::deque<index_t> facets_to_paint;
            index_t current_facet = picked_facet_id_;
            index_t adjacent_facet = index_t(-1);
            unsigned int count_facets_changed = 1;

            FOR(le,3) {
                adjacent_facet = mesh_.facets.adjacent(current_facet,le);
                if((label[adjacent_facet] == previous_label) && !VECTOR_CONTAINS(facets_to_paint,adjacent_facet)) {
                    facets_to_paint.push_back(adjacent_facet);
                }
            }

            while(!facets_to_paint.empty()) {
                current_facet = facets_to_paint.front();
                facets_to_paint.pop_front();

                label[current_facet] = selected_label_; // change label of picked facet
                count_facets_changed++;

                FOR(le,3) {
                    adjacent_facet = mesh_.facets.adjacent(current_facet,le);
                    if((label[adjacent_facet] == previous_label) && !VECTOR_CONTAINS(facets_to_paint,adjacent_facet)) {
                        facets_to_paint.push_back(adjacent_facet);
                    }
                }
            }

            fmt::println(Logger::out("painting"),"Label of {} facets changed with bucket fill tool",count_facets_changed); Logger::out("painting").flush();
        }
    }

private:

    int paint_mode_;                                // see PAINT_MODE_* macros
    index_t selected_label_;                        // label to paint
    Memory::byte picked_facet_id_as_pixel_[4];      // intermediate representation of the picked facet id
    index_t picked_facet_id_;                       // last picked facet
};

int main(int argc, char** argv) {
    LabelingPainterApp app;
	app.start(argc,argv);
    return 0;
}