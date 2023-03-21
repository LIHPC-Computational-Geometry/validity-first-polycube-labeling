#include <geogram_gfx/gui/colormaps/french.xpm>
#include <geogram_gfx/gui/colormaps/black_white.xpm>
#include <geogram_gfx/gui/colormaps/viridis.xpm>
#include <geogram_gfx/gui/colormaps/rainbow.xpm>
#include <geogram_gfx/gui/colormaps/cei_60757.xpm>
#include <geogram_gfx/gui/colormaps/inferno.xpm>
#include <geogram_gfx/gui/colormaps/magma.xpm>
#include <geogram_gfx/gui/colormaps/parula.xpm>
#include <geogram_gfx/gui/colormaps/plasma.xpm>
#include <geogram_gfx/gui/colormaps/blue_red.xpm>

#include "colormaps_array.h"

const char colormap_name[LAST_COLORMAP+1][15] = {
	"french",
	"black_white",
	"viridis",
	"rainbow",
	"cei_60757",
	"inferno",
	"magma",
	"parula",
	"plasma",
	"blue_red",
};

const char **colormap_xpm[LAST_COLORMAP+1] = {
	french_xpm,
	black_white_xpm,
	viridis_xpm,
	rainbow_xpm,
	cei_60757_xpm,
	inferno_xpm,
	magma_xpm,
	parula_xpm,
	plasma_xpm,
	blue_red_xpm,
};