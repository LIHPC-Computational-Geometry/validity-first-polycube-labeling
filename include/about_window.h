#pragma once

#include <fmt/core.h>

void draw_about_window(std::string name, std::string geogram_version) {
    // modification of ext/geogram/src/lib/geogram_gfx/gui/simple_application.cpp draw_about()
    ImGui::Separator();
    if(ImGui::BeginMenu(icon_UTF8("info") + " About...")) {
        ImGui::Text("%s : a GEOGRAM application",name.c_str());
        ImGui::Text("https://github.com/BrunoLevy/geogram");
        ImGui::Text("GEOGRAM version:%s",geogram_version.c_str());
         
        ImGui::Separator();
        
        ImGui::Text("Projects and organizations behind this application:");
        ImGui::Text("\n");

        ImGui::Text("SALOME - The open source platform for numerical simulation");
        ImGui::Text("https://www.salome-platform.org/");
        
        ImGui::Text("\n");

        ImGui::Text("CEA - French Alternative Energies and Atomic Energy Commission");
        ImGui::Text("https://www.cea.fr/english/Pages/Welcome.aspx");
       
        ImGui::Text("\n");

        ImGui::Text("LIHPC - High Performance Computing Laboratory for Calculation and Simulation");
        ImGui::Text("https://www-hpc.cea.fr/en/LIHPC.html");

        ImGui::EndMenu();
    }
}