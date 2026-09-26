# CCI Import Guide

OpenRdson extracts on-resistance from a **CCI database** — the set of files a
Calibre LVS run emits through the *Calibre Connectivity Interface* (CCI). This
guide walks a new user through producing that database from an LVS run and
wiring it into OpenRdson.

If you already have a CCI database, skip straight to
[Configuring OpenRdson](#configuring-openrdson).

---

## 1. What CCI is

CCI (Calibre Connectivity Interface) is Mentor/Siemens' standard export of an LVS
result: the extracted netlist, the port/terminal list, the device table, the
layer mapping, and the layout geometry, in a set of plain-text files that other
tools can read without talking to the Calibre database directly.

OpenRdson consumes these files, plus two technology inputs (the ICT process stack
and the logical→physical layer map) and an optional channel-model table.

## 2. Files OpenRdson needs

| File               | Source | Purpose |
|------|--------|---------|
| `*.agf` (or GDSII) | Calibre | The layout geometry (AGF is Calibre's ASCII layout format) |
| `*.gds.map`        | Calibre | GDS layer/datatype number → logical layer name |
| `*.ports`          | Calibre | Top-level ports (name, net, position, layer) |
| `*.devtab`         | Calibre | Device table (`seed` layer → model → ordered terminals) |
| `*.spi`            | Calibre | Golden LVS SPICE netlist (instances with `$X`/`$Y`, device params) |
| `*.lnn`            | Calibre | Layout net-name table (node id → net name) |
| `*.pin_xy_spi`     | Calibre | CCI netlist with the `.DEVTMPLT` terminal-order dictionary |
| `*.ict`            | Foundry | Process stack: conductors, dielectrics, vias, resistivities, z-heights |
| `*.map` (cci)      | Foundry | Logical → physical layer map used for extraction |
| `model.csv`        | Characterization | Bias-dependent channel model (`Id(T, Vgs, Vds)` table) |

The first seven are produced by Calibre; the last three come from the foundry's
process package or your own characterization and are not part of the LVS output.

### File formats at a glance

These are the conventions OpenRdson's parsers expect (matching Calibre CCI):

- **`gds.map`** — one line per layer: `LogicalName GDSLayer Datatype`, e.g.
  `Metal1 45 0`, `seed_device 45 0`.
- **`ports`** — SVDB Port Table; one row per top-level port with name, net,
  x/y, and layer.
- **`devtab`** — the Calibre `Device_Table`; a `Device Entry` per template with
  its `seed` layer, device kind, model/subcircuit name, ordered terminal list,
  and property layers.
- **`spi`** — one `X`-instance line per device, e.g.
  `XX0 1 3 2 1 <model> w=1e-05 weff=0.001 nf=100 nx=10 ny=10 m=1 $X=23680 $Y=10280 $D=202`.
  The `$X`/`$Y` (in database units) locate the device; `$D` references the
  devtab template id.
- **`lnn`** — `nodeId netName` lines.
- **`pin_xy_spi`** — contains the `.DEVTMPLT <id> <model>() <seed> <layer>(<terminal>) ...`
  dictionary that defines each model's terminal order (d/g/s/sub).

---

## 3. Prerequisites

- A Calibre LVS run that completes cleanly (the LVS report shows a match).
- CCI output enabled in the rule deck or the run command. (`CCI QUERY`)
- The foundry ICT (`*.ict`) and layer map (`cci.map`) for the process.
- (Optional) a `model.csv` channel-model table if you extract active devices.

## 4. Producing the CCI database from an LVS run

There are two equivalent ways to get the CCI files:

1. **Enable CCI output in the LVS rule deck** by adding
   `MASK SVDB DIRECTORY <SVDB_DIR> CCI QUERY` directives, or
2. **Run the Calibre Query Server** against the LVS results database (SVDB) and
   query out the CCI files with a script (see below).

The CCI output normally drops the seven Calibre files into the LVS output
directory, named after the run (e.g. `top.gds.map`, `top.ports`,
`top.pin_xy_spi`, `top.lnn`, `top.spi`, `top.devtab`).

## 5. Example Calibre Query Server script

The following query runs in the Calibre Query Server's native format — plain-text
commands, one per line.  Adjust command names if your Calibre version
differs.

Run the below script with the following command in your Unix environment:
`calibre -query_input CCI_script.txt -query <SVDB_DIR> <TOP_CELL> > query.log 2>&1`

```
#specifies GDS property values for attributes, corresponding to nets, instances and #devices

gds netprop number 5

gds placeprop number 6

gds devprop number 7

#Writes GDS map file

response file top.gds.map

gds map

response direct

#Writes Anntoated GDS file (AGF)

gds write top.agf

#Include trivial pins and empty cells in the layout netlist

#Use only node numbers for netlists

layout netlist trivial pins YES

layout netlist empty cells YES

layout netlist names NONE

#Writes the Node-tp-Net Name mapping

layout nametable write top.lnn

#Writes the Layout netlist

layout netlist hierarchy AGF

layout netlist write top.spi

#Writes a layout netlist with AGF hierarchy and $PIN_XY info

layout netlist pin locations YES

layout netlist write top.pin_xy_spi

#source and layout placement(cell instances) hierarchy(sph, lph)

source hierarchy write top.sph

layout hierarchy write top.lph

#Writes cross-refference file

net xref write top.nxf

instance xref write top.ixf

#Write top level and cell-level port tables

port table write top.ports

port table cells write top.ports_cells

# Write device table

response file top.devtab

device table

response direct
```

## 6. Verifying the output

Check that each file exists and looks right:

```sh
head top.gds.map      # LogicalName GDSLayer Datatype
head top.ports        # ports with net + position
head top.lnn          # nodeId netName
grep DEVTMPLT top.pin_xy_spi   # terminal-order dictionary present
grep '^X' top.spi     # device instances with $X/$Y and $D
```

Then copy the ICT, layer map, and model table alongside them.

## 7. Configuring OpenRdson

Point the config at the files:

```yaml
project:
  cci_dir: <path_to_cci_files>     # dir holding the 7 Calibre files
  layout: top.agf
  gds_map: top.gds.map
  ports: top.ports
  devtab: top.devtab
  spi: top.spi
  tech_ict: <path_to>.ict
  layer_map: <path_to>cci.map
  model_csv: <path_to>model.csv
```

Generate the template with `openrdson --print-default-config`, fill it in, then
run:

```sh
openrdson --config openrdson.yaml sheet-rds
```

If the device recognition or terminal resolution complains, the most common
causes are a missing/mismatched `devtab` (wrong `$D` template id) or a `cci.map`
that doesn't name every layer the extraction needs.
