#include "LabelingViewerApp.h"

// values for paint_mode_
#define PAINT_MODE_DISABLED     0
#define PAINT_MODE_PENCIL       1
#define PAINT_MODE_BUCKET_FILL  2

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
        current_paint_mode_ = PAINT_MODE_DISABLED;
        previous_paint_mode_ = PAINT_MODE_DISABLED;
        selected_label_ = 0;
    }

private:

    void draw_scene() override {
        // if user change state/visu mode -> disable paint mode
        if( (previous_state_ != current_state_) || (previous_labeling_visu_mode_ != current_labeling_visu_mode_) )  {
            current_paint_mode_ = PAINT_MODE_DISABLED;
        }
        // if user activate paint mode -> set raw labeling visu mode
        if( (previous_paint_mode_ != current_paint_mode_) && (current_paint_mode_ != PAINT_MODE_DISABLED) ) {
            current_labeling_visu_mode_ = VIEW_RAW_LABELING;
        }
        LabelingViewerApp::draw_scene();
    }

	// add buttons for labeling operators on the "object properties" panel
    void draw_object_properties() override {
		LabelingViewerApp::draw_object_properties();
		if(current_state_ == LabelingViewerApp::State::labeling) {
			ImGui::Separator();

            ImGui::Text("Paint mode");
            ImGui::RadioButton("Disabled",&current_paint_mode_,PAINT_MODE_DISABLED);
            ImGui::RadioButton("Pencil",&current_paint_mode_,PAINT_MODE_PENCIL);
            ImGui::RadioButton("Bucket fill",&current_paint_mode_,PAINT_MODE_BUCKET_FILL);

            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 4);
            if (ImGui::BeginCombo("Label", label_selection_[selected_label_]))
            {
                // based on ext/geogram/src/lib/geogram_gfx/third_party/imgui/imgui_demo.cpp "combo 1" code snippet
                for (int n = 0; n < 6; n++) // for each label in [0,5]
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
		}
	}

private:

    int previous_paint_mode_, current_paint_mode_;
    int selected_label_;
};

int main(int argc, char** argv) {
    LabelingPainterApp app;
	app.start(argc,argv);
    return 0;
}