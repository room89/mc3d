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

## Ограничение шага и разрешение (этап 6)

`collision_courant` по умолчанию равен 0.1: шаг ограничивается средней
частотой столкновений на частицу, оценённой сверху по RMS относительных
скоростей. Значение 0 отключает это ограничение. `cu` применяется только
к транспортному пределу. Старые конфиги получают новый ограничитель;
число шагов и траектории случайных столкновений могут измениться.

```bash
make run CONFIG=configs/checks/resolution_dense_10.json OUTPUT_DIR=output/resolution-dense-10 LOAD_MODE=balanced
make run CONFIG=configs/checks/resolution_dense_05.json OUTPUT_DIR=output/resolution-dense-05 LOAD_MODE=balanced
make run CONFIG=configs/checks/resolution_axes_10.json OUTPUT_DIR=output/resolution-axes LOAD_MODE=balanced
make run CONFIG=configs/checks/cut_cube.json OUTPUT_DIR=output/resolution-cut LOAD_MODE=balanced
make check-resolution OUTPUT_DIR=output/resolution-comparison
```

`resolution_diagnostics=false` отключает только вывод. `resolution_steps.csv`
содержит ограничения, фактический шаг, причину выбора, индекс ячейки/оси
(оси 0/1/2 = X/Y/Z), либо индекс граничного резервуара, и гистограммы.
`resolution_cells_<step>.csv` содержит подробности первого, последнего
и шагов со снимками; при остановке сохраняется последний завершённый шаг.
Данные относятся к состоянию до шага (`time`), `end_time` — его окончание.
При равенстве ограничений приоритет причины: конец расчёта, столкновения,
перенос; среди одинаковых ячеек выбирается первый индекс.

`nu_bound` — верхняя оценка средней частоты, а не максимум частоты каждой
частицы. `nu_eq` и `lambda_eq` — максвелловские оценки для диагностики;
они не предполагают, что текущий газ действительно равновесен. Для твёрдых,
малозаселённых, холодных ячеек и отключённых столкновений неприменимые
поля пусты, причина записана в `status`. Бесконечные ограничения также
записываются пустыми полями. В SI длины умножаются на `length_m`, времена —
на `time_s`, частоты делятся на `time_s`; отношения безразмерны.

Сравнение двух пучков использует 200 частиц, Kn=0.1, десять seed 123…132,
три момента времени и ограничения 0.1/0.05. Отчёт в `validation/stage6/`
содержит исходные моменты и допуски: max(5 SE парной разности, 0.06).
Прохождение сравнения подтверждает только отсутствие различимого эффекта
шага в этом опыте, не полную сходимость пространственно неоднородного течения.

## Оконная статистика и балансы (этап 7)

Новые измерения включены по умолчанию. `statistics_enabled=false` отключает
оконные таблицы, `balance_enabled=false` — файлы балансов. Физическое
обновление динамической стенки работает независимо от этих параметров.

`statistics_windows` задаёт список `[start,end]` в единицах времени решателя,
например `[[0,0.025],[0.025,0.05],[0.0125,0.0375]]`. Перекрытия разрешены;
границы конечны и удовлетворяют `0 <= start < end <= end_time`. Отсутствие
параметра или пустой список в JSON означают весь расчёт. В SI используется
`statistics_windows_s`; границы переводятся автоматически.

```bash
make run CONFIG=configs/checks/statistics_closed.json OUTPUT_DIR=output/statistics-closed LOAD_MODE=balanced
make run CONFIG=configs/checks/statistics_open.json OUTPUT_DIR=output/statistics-open LOAD_MODE=balanced
make run CONFIG=configs/checks/statistics_free.json OUTPUT_DIR=output/statistics-free LOAD_MODE=balanced
make run CONFIG=configs/checks/statistics_mirror.json OUTPUT_DIR=output/statistics-mirror LOAD_MODE=balanced
make run CONFIG=configs/checks/statistics_fixed_wall.json OUTPUT_DIR=output/statistics-fixed LOAD_MODE=balanced
make run CONFIG=configs/checks/statistics_dynamic_wall.json OUTPUT_DIR=output/statistics-dynamic LOAD_MODE=balanced
ctest --test-dir build --output-on-failure -R 'StatisticalMomentsTest|StatisticsIntegration'
```

`cell_statistics_<window>.csv` содержит моменты, среднюю популяцию,
плотность, скорость, `Tstat`, шесть компонент тензора давления и тепловой
поток. Вес наблюдения — пересечение шага с окном; используется состояние
до шага. Пустая ячейка увеличивает длительность наблюдения, но не A0.
A0 равен интегралу числа модельных частиц по времени. A1/A2/A3 накапливаются
в фиксированной системе скоростей с началом `origin`, указанным в таблице;
это сырые моменты в этой системе, а не сумма мгновенных центральных моментов.
Производные скорость и тепловой поток возвращаются в физические соглашения
решателя. `P` положителен при сжатии; напряжения Коши равны `-P`.

Новая температура получается из объединённых моментов. Она не усредняет
старое поле `T` и не меняет его смысл. Универсальная поправка Бесселя не
применяется: наблюдения одних частиц во времени коррелированы. Прохождение
теста независимых максвелловских выборок не делает произвольную конечную
коррелированную выборку строго несмещённой.

`surface_statistics_<window>.csv` содержит переданные телу импульс и тепло,
среднюю силу, мощность и тепловой поток на площадь треугольника. Знак плюс —
от газа к телу. Столбец `wall_temperature_at_output` — текущая температура
при записи таблицы; остальные интегралы относятся к указанному окну.
Для части шага обмен учитывается с весом длительность пересечения / dt.
Это среднее дискретных шаговых нагрузок, не реконструкция времени ударов.

`conservation.csv` содержит запасы газа до/после шага, суммарные вход/выход,
обмен с телом и зеркалами, изменение энергии стенки и невязки. Периодические
переносы не являются внешними потоками. Вход hyperfree учитывается до
возможного удара о тело. Выход учитывается при фактическом удалении частиц;
маршрутизация по отдельным граням не реконструируется. Необъяснимый выход
через замкнутые границы вызывает ошибку вместо маскировки потерянных частиц.

Rmass = ΔMgas - Min + Mout; Rp = ΔPgas - Pin + Pout + Jbody + Jmirror;
RE = ΔEgas - Ein + Eout + Qbody + Qmirror.
Для динамической стенки Rwall = ΔUwall - Qbody; для фиксированной Qthermostat
равно Qbody, а изменение энергии материала равно нулю. Положительный
Qthermostat означает отвод тепла термостатом, отрицательный — подвод.
`conservation_summary.json` хранит накопленные потоки и невязки, максимумы
шаговых энергетических/импульсных невязок и изменения энергии стенки.

Масса и импульс нормированы общим весом модельной частицы
`1/reference_particle_density`. Этот вес не перенормируется по текущей
популяции. В SI используйте масштабы `mass_kg`, `momentum_kg_m_s`,
`energy_j`, `power_w`, `pressure_pa`, `heat_flux_w_m2`, `velocity_m_s` и
`temperature_k` из паспорта. Отсутствующие значения выводятся пустыми,
состояния ячеек — `solid`, `unobserved`, `empty`, `valid`.

Оконные файлы атомарно обновляются при снимках, закрытии окна и завершении
или остановке. Флаг `complete` и `observed_time` отличают частичный результат.
При отсутствии шагов баланс содержит только заголовок, окна — нулевую
экспозицию. Существующие форматы снимков и body_loads.csv сохранены.
Дополнительная память растёт как число ячеек × число окон; на каждую пару
хранятся несколько сотен байт моментов. Обмены с телом буферизуются в ячейках
и сводятся в фиксированном порядке после завершения рабочих задач.

Отчёт `validation/stage7/report.json` содержит результаты шести запусков,
оценки температуры на независимых выборках, хеши исходников и бинарника.
В локальных запусках невязка числа частиц равна нулю, максимальная абсолютная
невязка энергии меньше 4.4e-15 в единицах решателя.
