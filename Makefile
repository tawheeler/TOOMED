########################################################################
##############################  Makefile  ##############################
########################################################################

EXE = TOOMED
IMGUI_DIR = src/imgui
SOURCES = src/toomed.cpp src/assets_exporter.cpp src/delaunay_mesh.cpp src/game_map.cpp src/geometry_utils.cpp src/math_utils.cpp
SOURCES += $(IMGUI_DIR)/imgui.cpp $(IMGUI_DIR)/imgui_demo.cpp $(IMGUI_DIR)/imgui_draw.cpp $(IMGUI_DIR)/imgui_tables.cpp $(IMGUI_DIR)/imgui_widgets.cpp
SOURCES += $(IMGUI_DIR)/backends/imgui_impl_sdl2.cpp $(IMGUI_DIR)/backends/imgui_impl_sdlrenderer2.cpp
OBJS = $(addsuffix .o, $(basename $(notdir $(SOURCES))))
OBJS = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(basename $(notdir $(SOURCES)))))
OBJS = obj/toomed.o obj/assets_exporter.o obj/delaunay_mesh.o obj/game_map.o obj/geometry_utils.o obj/math_utils.o
OBJS += obj/imgui/imgui.o obj/imgui/imgui_demo.o obj/imgui/imgui_draw.o obj/imgui/imgui_tables.o obj/imgui/imgui_widgets.o
OBJS += obj/imgui/backends/imgui_impl_sdl2.o obj/imgui/backends/imgui_impl_sdlrenderer2.o

UNAME_S := $(shell uname -s)

CXXFLAGS = -std=c++17 -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends
CXXFLAGS += -g -Wall -Wformat
LIBS = /usr/local/lib/libSDL2.a -lpthread

##---------------------------------------------------------------------
## BUILD FLAGS PER PLATFORM
##---------------------------------------------------------------------

# $(info SOURCES: $(SOURCES))
# $(info OBJS: $(OBJS))

ifeq ($(UNAME_S), Linux) #LINUX
	ECHO_MESSAGE = "Linux"
	LIBS += -lGL -ldl `sdl2-config --libs`

	CXXFLAGS += `sdl2-config --cflags`
	CFLAGS = $(CXXFLAGS)
endif

ifeq ($(UNAME_S), Darwin) #APPLE
	ECHO_MESSAGE = "Mac OS X"
	LIBS += -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo `sdl2-config --libs`
	LIBS += -L/usr/local/lib -L/opt/local/lib

	CXXFLAGS += `sdl2-config --cflags`
	CXXFLAGS += -I/usr/local/include -I/opt/local/include
	CFLAGS = $(CXXFLAGS)
endif

ifeq ($(OS), Windows_NT)
	ECHO_MESSAGE = "MinGW"
	LIBS += -lgdi32 -lopengl32 -limm32 `pkg-config --static --libs sdl2`

	CXXFLAGS += `pkg-config --cflags sdl2`
	CFLAGS = $(CXXFLAGS)
endif

# $(info LIBS: $(LIBS))
# $(info CXXFLAGS: $(CXXFLAGS))

##---------------------------------------------------------------------
## BUILD RULES
##---------------------------------------------------------------------

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	
all: $(EXE)
	@echo Build complete for $(ECHO_MESSAGE)

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(EXE) $(OBJS)
