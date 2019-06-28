#!/bin/sh

walk_dir () {
    for pathname in "$1"/*; do
        if [ -d "$pathname" ]; then
            walk_dir "$pathname"
        elif [ -e "$pathname" ]; then
            case "$pathname" in
                *.hpp|*.cpp)
                    printf '%s\n' "$pathname"
                    clang-format -i "$pathname"
            esac
        fi
    done
}

SRC_DIR=./src

walk_dir "$SRC_DIR"
