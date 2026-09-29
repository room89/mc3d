BUILD_DIR ?= build
BUILD_TYPE ?= Release
BUILD_TESTING ?= ON
CMAKE ?= cmake
CTEST ?= ctest
CONFIG ?= configs/example.json
ARGS ?=

.PHONY: all help configure build test run run-example run-smoke clean

all: build

help:
	@printf '%s\n' \
	  'make build                         Configure and build MC3dSolver' \
	  'make test                          Build and run the CTest suite' \
	  'make run CONFIG=path [ARGS="..."] Run with a configuration file' \
	  'make run-example                   Run configs/example.json' \
	  'make run-smoke                     Run a small built-in sample case' \
	  'make clean                         Remove the build directory'

configure:
	$(CMAKE) -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DBUILD_TESTING=$(BUILD_TESTING)

build: configure
	$(CMAKE) --build $(BUILD_DIR) --parallel

test: BUILD_TESTING = ON
test: configure
	$(CMAKE) --build $(BUILD_DIR) --target mc3d_tests --parallel
	$(CTEST) --test-dir $(BUILD_DIR) --output-on-failure

run: build
	$(BUILD_DIR)/MC3dSolver --config $(CONFIG) $(ARGS)

run-example: CONFIG = configs/example.json
run-example: run

run-smoke: build
	$(BUILD_DIR)/MC3dSolver --lx 1 --ly 1 --lz 1 --ncx 2 --ncy 2 --ncz 2 --np 1 --end-time 0.01 --geometry-type none

clean:
	$(CMAKE) -E rm -rf $(BUILD_DIR)
