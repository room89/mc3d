BUILD_DIR ?= build
BUILD_TYPE ?= Release
BUILD_TESTING ?= ON
CMAKE ?= cmake
CTEST ?= ctest
CONFIG ?= configs/example.json
OUTPUT_DIR ?= output
LOAD_MODE ?=
CONTROL_FILE ?= $(OUTPUT_DIR).control
PYTHON ?= python3
ARGS ?=

.PHONY: all help configure build test run run-example run-smoke load clean

all: build

help:
	@printf '%s\n' \
	  'make build                         Configure and build MC3dSolver' \
	  'make test                          Build and run the CTest suite' \
	  'make run CONFIG=path OUTPUT_DIR=new_dir [LOAD_MODE=balanced] [ARGS="..."]' \
	  'make load OUTPUT_DIR=dir LOAD_MODE=quiet  Change a running calculation load' \
	  'LOAD_MODE: full, balanced, quiet, pause; CONTROL_FILE defaults to OUTPUT_DIR.control' \
	  'make run-example                   Run configs/example.json' \
	  'make run-smoke                     Run a small built-in sample case' \
	  '                                    Simulation files go to $(OUTPUT_DIR)/' \
	  'make clean                         Remove the build directory'

configure:
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DBUILD_TESTING=$(BUILD_TESTING)

build: configure
	$(CMAKE) --build "$(BUILD_DIR)" --parallel

test: BUILD_TESTING = ON
test: configure
	$(CMAKE) --build "$(BUILD_DIR)" --target mc3d_tests MC3dSolver --parallel
	$(CTEST) --test-dir "$(BUILD_DIR)" --output-on-failure

load:
	$(PYTHON) tools/control_simulation.py . "$(LOAD_MODE)" --file "$(CONTROL_FILE)" --create-parent

run: build
ifneq ($(strip $(LOAD_MODE)),)
	$(MAKE) load LOAD_MODE="$(LOAD_MODE)" CONTROL_FILE="$(CONTROL_FILE)"
endif
	"$(BUILD_DIR)/MC3dSolver" --config "$(CONFIG)" --output-dir "$(OUTPUT_DIR)" $(if $(strip $(LOAD_MODE)),--runtime-control "$(CONTROL_FILE)") $(ARGS)

run-example: CONFIG = configs/example.json
run-example: run

run-smoke: build
ifneq ($(strip $(LOAD_MODE)),)
	$(MAKE) load LOAD_MODE="$(LOAD_MODE)" CONTROL_FILE="$(CONTROL_FILE)"
endif
	"$(BUILD_DIR)/MC3dSolver" --output-dir "$(OUTPUT_DIR)" --lx 1 --ly 1 --lz 1 --ncx 2 --ncy 2 --ncz 2 --np 1 --end-time 0.01 --geometry-type none $(if $(strip $(LOAD_MODE)),--runtime-control "$(CONTROL_FILE)") $(ARGS)

clean:
	$(CMAKE) -E rm -rf $(BUILD_DIR)

ifneq ($(strip $(LOAD_MODE)),)
ifeq ($(filter $(LOAD_MODE),full balanced quiet pause),)
$(error LOAD_MODE must be full, balanced, quiet or pause)
endif
ifneq ($(findstring --runtime-control,$(ARGS)),)
$(error Use LOAD_MODE or ARGS=--runtime-control, not both)
endif
endif

.PHONY: check-resolution
check-resolution: build
	$(CMAKE) -E make_directory "$(OUTPUT_DIR)"
	"$(BUILD_DIR)/mc3d_resolution_comparison" "$(OUTPUT_DIR)/resolution_comparison.csv"

PHYSICS_MODE ?= quick
.PHONY: check-physics
check-physics: build
	$(PYTHON) tools/check_physics.py --binary "$(BUILD_DIR)/MC3dSolver" --mode "$(PHYSICS_MODE)" --output-dir "$(OUTPUT_DIR)" --load-mode "$(if $(LOAD_MODE),$(LOAD_MODE),balanced)"
