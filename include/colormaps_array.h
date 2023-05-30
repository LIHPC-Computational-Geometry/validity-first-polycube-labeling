#pragma once

// define an ordered list of Geogram colormaps, to allow indexing with macros

#define COLORMAP_FRENCH			0
#define COLORMAP_BLACK_WHITE	1
#define COLORMAP_VIRIDIS		2
#define COLORMAP_RAINBOW		3
#define COLORMAP_CEI_60757		4
#define COLORMAP_INFERNO		5
#define COLORMAP_MAGMA			6
#define COLORMAP_PARULA			7
#define COLORMAP_PLASMA			8
#define COLORMAP_BLUE_RED		9

#define LAST_XPM_COLORMAP (COLORMAP_BLUE_RED)

#define COLORMAP_LABELING       10

#define TO_GL_TEXTURE_INDEX(colormap_index) (colormap_index+2) // it seems apps based on SimpleApplication have 2 GL textures defined before the colormaps. One of them is the Geogram logo ? see SimpleApplication::GL_initialize()

extern const char colormap_name[LAST_XPM_COLORMAP+1][15];

extern const char **colormap_xpm[LAST_XPM_COLORMAP+1];