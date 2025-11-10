# MC3D

MC3D — численная библиотека для моделирования свободномолекулярных течений на основе
метода Монте-Карло.

## Сборка

```bash
cmake -S /home/lostpointer/mc3d -B /home/lostpointer/mc3d/build
cmake --build /home/lostpointer/mc3d/build
```

## Запуск

Выполнимый файл появится в `build/mc3d`:

```bash
/home/lostpointer/mc3d/build/mc3d
```

Программа интерактивно запрашивает параметры расчёта (Kn, Cu, размеры сетки и т. д.)
и сохраняет результаты в файл `data.dat`.

## Требования

- CMake 3.12+
- Компилятор с поддержкой C++23

## Форматирование

Для приведения исходников к Google C++ Style используйте:

```bash
find /home/lostpointer/mc3d -path /home/lostpointer/mc3d/build -prune -o \( -name '*.cpp' -o -name '*.h' \) -print0 | \
  xargs -0 clang-format -style=Google -i
```

