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
./build/MC3dSolver
```

По умолчанию пример в `main.cpp` создаёт клиновидное тело, настраивает кластер
ячееек, задаёт граничные условия и запускает расчёт на конечный промежуток
времени. В ходе выполнения в корне проекта появятся файлы `data*.dat`,
`speed.dat`, `time.dat`, `*.net`.

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
  накладывая векторную карту скоростей:
  ```bash
  python3 tools/plot_snapshot.py --input build/data.dat --plane xy \
      --z-value 0.0 --thickness 0.05 --field T --show-cbar \
      --show-velocity --quiver-step 3 --quiver-color white
  ```
  Аргументы `--show-velocity`, `--quiver-step`, `--quiver-scale` и
  `--quiver-color` управляют наложением стрелок скорости на выбранной плоскости.
- `tools/generate_gif.py` — собирает последовательность файлов `data*.dat`
  в анимированный GIF. Поддерживает те же настройки среза и отображения, что и
  `plot_snapshot.py`. Пример:
  ```bash
  python3 tools/generate_gif.py --pattern "build/data*.dat" --plane xy \
      --z-value 0.0 --thickness 0.05 --field T --show-velocity \
      --quiver-step 3 --fps 8 --out plots/run.gif
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
