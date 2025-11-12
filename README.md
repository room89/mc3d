Monte-Carlo 3D Flow Solver
===========================

MC3D — экспериментальный решатель разреженных потоков, основанный на
стохастическом моделировании движения и столкновений частиц в трёхмерной
области. Код исторически писался для исследовательских задач и сейчас находится
в процессе модернизации (перевод на современный стиль C++, чистка API и
инфраструктуры сборки).

## Возможности
- геометрия твёрдого тела, задаваемая наборами полигонов (генерация клина,
  пирамиды, куба либо импорт STL);
- сетка ячеек (`CellCluster`) с контролем граничных условий
  (зеркальные, периодические, свободные, «гипер-свободные»);
- моделирование движения частиц, столкновения, сбор статистики по ячейкам;
- запись данных по плотности, скоростям, времени шага.

## Сборка
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Исполняемый файл `MC3dSolver` находится в директории `build/`.

## Запуск
```bash
./build/MC3dSolver [опции]
```

Все параметры расчёта теперь можно задать либо из командной строки, либо через
конфигурационный файл. Если параметры не указаны, используются значения по
умолчанию (аналогичные старому сценарию в `main.cpp` с гипер-свободными
границами).

### Примеры CLI
- задать размеры области и число ячеек:
  ```bash
  ./build/MC3dSolver --lx 2.0 --ly 1.0 --lz 0.5 --ncx 40 --ncy 20 --ncz 10
  ```
- выбрать тип граничных условий:
  ```bash
  ./build/MC3dSolver --boundary-x-neg mirror --boundary-x-pos mirror \
      --boundary-y-neg periodic --boundary-y-pos periodic \
      --boundary-z-neg free --boundary-z-pos free
  ```
- выбрать геометрию тела:
  ```bash
  ./build/MC3dSolver --geometry-type wedge --geometry-wedge-alpha 0.349
  ```

### Конфигурационный файл
Любую опцию можно перенести в JSON-файл:
```json
{
  "lx": 1.5,
  "ly": 1.0,
  "lz": 1.0,
  "ncx": 24,
  "ncy": 16,
  "ncz": 8,
  "np": 80,
  "kn": 0.03,
  "cu": 0.5,
  "end_time": 2.0,
  "boundary_x_neg": "mirror",
  "boundary_x_pos": "mirror",
  "geometry_type": "cube",
  "geometry_cube_width": 0.4,
  "geometry_cube_length": 0.4,
  "geometry_cube_height": 0.4
}
```

Запуск с файлом:
```bash
./build/MC3dSolver --config sample_config.json
```
Опции из командной строки имеют приоритет над значениями из файла и могут
дополнять несколько файлов `--config`. См. также готовый пример
`configs/example.json`.

Помимо JSON поддерживается упрощённый YAML (плоские ключи вида `ключ: значение`)
и старый формат `ключ = значение` — все варианты можно комбинировать.

В ходе выполнения в корне проекта появятся файлы `data*.dat`, `speed.dat`,
`time.dat`, `*.net`.

## Структура исходников
- `src/cell.*`, `src/cell_cluster.*` — ядро solver'а (ячейки, кластеры, шаги по
  времени);
- `src/boundary*`, `src/free_boundary*` — граничные условия;
- `src/geometry.*`, `src/poligon.*` — построение и обработка геометрии;
- `src/particle.*` — модель частицы и логика столкновений;
- `src/utils/` — вспомогательный код (логгер, таймер, потоковый пул, генератор
  случайных чисел);
- `src/warning.*`, `src/exception.*`, `src/exit_code.h` — обработка исключительных
  ситуаций;
- `tools/clang-format.sh` — массовое форматирование.

## Визуализация результатов

Solver записывает результаты в текстовые файлы с разделителем `;`:
`data.dat`, `data<t>.dat`, `speed.dat`, `time.dat`. Для их оперативного просмотра
добавлены вспомогательные скрипты на Python 3 (`numpy`, `pandas`, `matplotlib`):

- `tools/plot_snapshot.py` — строит срез по выбранной плоскости (`xy`, `xz`, `yz`)
  и отображает выбранное поле (`ro`, `T`, `vx`, `vy`, `vz`, `E`), при желании
  накладывая векторную карту скоростей. Для управления разрешением можно задать
  `--dpi` (плотность пикселей) и/или `--figsize WIDTH HEIGHT` (размер холста в
  дюймах):
  ```bash
  python3 tools/plot_snapshot.py --input build/data.dat --plane xy \
      --z-value 0.0 --thickness 0.05 --field T --show-cbar \
      --dpi 300 --figsize 10 8 \
      --show-velocity --quiver-step 3 --quiver-color white
  ```
  Аргументы `--show-velocity`, `--quiver-step`, `--quiver-scale` и
  `--quiver-color` управляют наложением стрелок скорости на выбранной плоскости.
- `tools/generate_gif.py` — собирает последовательность файлов `data*.dat`
  в анимированный GIF. Поддерживает те же настройки среза и отображения, что и
  `plot_snapshot.py`, включая `--dpi` и `--figsize` для контроля разрешения кадров.
  Пример:
  ```bash
  python3 tools/generate_gif.py --pattern "build/data*.dat" --plane xy \
      --z-value 0.0 --thickness 0.05 --field T --show-velocity \
      --quiver-step 3 --fps 8 --dpi 200 --figsize 12 9 --out plots/run.gif
  ```
  Дополнительно можно задать `--frame-step`, `--max-frames`, `--vmin/--vmax` для
  контроля длительности и цветовой шкалы.
- `tools/monitor_simulation.py` — следит за последним файлом `data<t>.dat` и
  обновляет график каждые несколько секунд (удобно запускать параллельно расчёту):
  ```bash
  python3 tools/monitor_simulation.py --directory build --plane xy --coord 0.0
  ```

Необходимые пакеты:
```bash
python3 -m pip install numpy pandas matplotlib pillow
```
