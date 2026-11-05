# Готовые проверочные запуски

Из корня репозитория (`make run` сам обновит сборку):

```bash
make run CONFIG=configs/checks/hs_low_population.json OUTPUT_DIR=output/check-hs LOAD_MODE=balanced
make run CONFIG=configs/checks/vhs_hs_limit.json OUTPUT_DIR=output/check-vhs-hs LOAD_MODE=balanced
make run CONFIG=configs/checks/vhs_argon_equilibrium.json OUTPUT_DIR=output/check-argon LOAD_MODE=balanced
make run CONFIG=configs/checks/none_cube_heating.json OUTPUT_DIR=output/check-wall LOAD_MODE=balanced
make run CONFIG=configs/checks/cut_cube.json OUTPUT_DIR=output/check-cut-cube LOAD_MODE=balanced
make run CONFIG=configs/checks/cut_wedge.json OUTPUT_DIR=output/check-cut-wedge LOAD_MODE=balanced
```

Каждый каталог должен быть новым или пустым. Все случаи используют seed 123
и один поток. `LOAD_MODE=balanced` добавляет отдых между шагами; для этих
однопоточных проверок число потоков остаётся равным одному. Ожидается `status: completed` в `manifest.json`, файлы
`initial.dat` и `final.dat` и конечные численные значения.

| Конфигурация | Что проверяет |
|---|---|
| `checks/hs_low_population.json` | HS с двумя начальными частицами на ячейку, периодическая область |
| `checks/vhs_hs_limit.json` | Предел VHS при ω=0,5; физические файлы должны совпасть с HS побайтово |
| `checks/vhs_argon_equilibrium.json` | Короткий равновесный VHS-запуск с профилем аргона и 20 частицами на ячейку |
| `checks/cut_cube.json` | Смещённый куб: газовые, пересечённые и твёрдые ячейки, сумма Vsolid = 0,091125 |
| `checks/cut_wedge.json` | Смещённый клин: сумма Vsolid = 0,44 × 0,6² × sin(0,4) |
| `checks/none_cube_heating.json` | Полёт без межмолекулярных столкновений, нагрузки на куб и изменение температуры стенки |

В случае `none_cube_heating` также появляются `body_loads.csv` и
`final.dat.surface.csv`. Проверка равновесия и частот с надёжной статистикой
выполняется тестами, а не оценкой одного короткого снимка:

```bash
ctest --test-dir build --output-on-failure -R 'CellVolumeTest|HsScatteringTest|CollisionFrequencyTest|LowPopulationCellTest|VhsTest'
ctest --test-dir build --output-on-failure -R 'RunManifestIntegration|ConfigExamplesIntegration'
```

## Управление нагрузкой

`LOAD_MODE` принимает `full`, `balanced`, `quiet`, `pause`. Без параметра
управление выключено, если оно не задано в конфигурации. При указанном
параметре создаётся соседний файл `OUTPUT_DIR.control`; каталог результатов
остаётся пустым до запуска. Режим можно менять из другого терминала:

```bash
make load OUTPUT_DIR=output/check-hs LOAD_MODE=quiet
make load OUTPUT_DIR=output/check-hs LOAD_MODE=pause
make load OUTPUT_DIR=output/check-hs LOAD_MODE=balanced
make load OUTPUT_DIR=output/check-hs LOAD_MODE=full
```

Для собственного пути передавайте одинаковый `CONTROL_FILE=...` обеим
командам. `make load` только меняет запрос, не собирает и не запускает solver.
`pause` удерживает живой процесс; SIGINT/SIGTERM завершают его и на паузе.
Короткая проверка может закончиться раньше ручного переключения — для
наблюдения используйте более длинный сценарий. Инициализация и сборка
этими режимами не ограничиваются.

## Физический случай в SI

```bash
python3 tools/physical_config.py configs/physical/argon_cube.json --output-dir output/prepared-argon
make run CONFIG=output/prepared-argon/solver.json OUTPUT_DIR=output/check-si LOAD_MODE=balanced
python3 tools/replay_run.py output/check-si/manifest.json --output-dir output/check-si-repeat
```

`physical/argon_cube.json` задаёт аргон при 300 K, n=10²⁰ м⁻³, скорости
500 м/с в области 1 см × 1 см × 1 см. Куб имеет размер 2 мм, время расчёта
20 мкс; включены теплоёмкость стенки и коэффициенты силы. Это короткая
проверка полного пути SI → расчёт → паспорт → повтор, не установившийся
аэродинамический эталон. Газовые объёмы куба и клина проверяются аналитическими тестами этапа 5.

Таблица `output/prepared-argon/units.json` содержит множители перевода
результатов в SI; её копия также включена в паспорт. `ro` умножается на
числовую плотность, `E` — на энергию молекулы; сила и поглощённое тепло
используют отдельные макроскопические масштабы.

## Остальные конфигурации

В корневых JSON явно заданы HS и seed 123; размеры и длительность сохранены.
`body_flow.json`, `hypersonic_cube.json`, `hypersonic_wedge.json` и
`benchmark_body_check.json` предназначены для более длительных расчётов
или измерения производительности. Для быстрой проверки используйте `checks/`.

`argon/argon_vhs.json` — набор параметров газа, а не полный проверочный
сценарий. `argon_vhs_report.json` — отчёт аппроксимации, не вход решателя.
Сейчас разрешены только отсутствие тела, встроенный куб и встроенный клин.
STL, цилиндр и пирамида отклоняются, в том числе при `collision_model=none`.
Тело должно лежать строго внутри области; касание внешней границы также
отклоняется. Допустимы перенос, положительный масштаб и дробление сетки.

## Проверка газовых объёмов

Во всех запусках `cell_geometry.csv` содержит геометрию ячеек в порядке снимка:
`index;x;y;z;dx;dy;dz;region;Vcell;Vgas;Vsolid;gas_fraction`.
`region` принимает значения `gas`, `cut`, `solid`. Объёмы заданы в единицах
решателя; для SI умножайте их на куб масштаба длины из `units.json`.
Паспорт содержит метод `convex_clip_v1`, числа ячеек каждого класса,
суммарные объёмы и SHA-256 таблицы. Для двух новых конфигураций ожидаются
`Vcell` в сумме 1, относительная ошибка суммы `Vsolid` не более 10⁻¹⁰,
нулевые `N` и `ro` в твёрдых ячейках и конечные значения во всех снимках,
включая `data.dat` и `final.net`.

Инициализация сохраняет вес модельной частицы: число частиц в пересечённой
ячейке случайно и в среднем равно `N0 × Vgas/Vcell`; оно не перенормирует N0.
Одиночный снимок с четырьмя частицами на полную ячейку не является проверкой
сходимости или локального равновесия. Эти конфигурации проверяют интеграцию
геометрии, нормировки и вывода; статистику частот проверяет `CellVolumeTest`.
