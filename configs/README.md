# Готовые проверочные запуски

Из корня репозитория после сборки `cmake --build build -j4`:

```bash
build/MC3dSolver --config configs/checks/hs_low_population.json --output-dir output/check-hs
build/MC3dSolver --config configs/checks/vhs_hs_limit.json --output-dir output/check-vhs-hs
build/MC3dSolver --config configs/checks/vhs_argon_equilibrium.json --output-dir output/check-argon
build/MC3dSolver --config configs/checks/none_cube_heating.json --output-dir output/check-wall
```

Каждый каталог должен быть новым или пустым. Все случаи используют seed 123
и один поток. Ожидается `status: completed` в `manifest.json`, файлы
`initial.dat` и `final.dat` и конечные численные значения.

| Конфигурация | Что проверяет |
|---|---|
| `checks/hs_low_population.json` | HS с двумя начальными частицами на ячейку, периодическая область |
| `checks/vhs_hs_limit.json` | Предел VHS при ω=0,5; физические файлы должны совпасть с HS побайтово |
| `checks/vhs_argon_equilibrium.json` | Короткий равновесный VHS-запуск с профилем аргона и 20 частицами на ячейку |
| `checks/none_cube_heating.json` | Полёт без межмолекулярных столкновений, нагрузки на куб и изменение температуры стенки |

В последнем случае также появляются `body_loads.csv` и
`final.dat.surface.csv`. Проверка равновесия и частот с надёжной статистикой
выполняется тестами, а не оценкой одного короткого снимка:

```bash
ctest --test-dir build --output-on-failure -R 'HsScatteringTest|CollisionFrequencyTest|LowPopulationCellTest|VhsTest'
ctest --test-dir build --output-on-failure -R 'RunManifestIntegration|ConfigExamplesIntegration'
```

## Физический случай в SI

```bash
python3 tools/physical_config.py configs/physical/argon_cube.json --output-dir output/prepared-argon
build/MC3dSolver --config output/prepared-argon/solver.json --output-dir output/check-si
python3 tools/replay_run.py output/check-si/manifest.json --output-dir output/check-si-repeat
```

`physical/argon_cube.json` задаёт аргон при 300 K, n=10²⁰ м⁻³, скорости
500 м/с в области 1 см × 1 см × 1 см. Куб имеет размер 2 мм, время расчёта
20 мкс; включены теплоёмкость стенки и коэффициенты силы. Это короткая
проверка полного пути SI → расчёт → паспорт → повтор, не установившийся
аэродинамический эталон. Точность газовых объёмов пересечённых ячеек ещё
предстоит проверить на этапе 5.

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
STL в SI требует поля `file_unit_m`: размер единицы координат файла в метрах.
