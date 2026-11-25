# Input file essentials

MOSE reads its main configuration from `input.ini` in the case directory. The file uses INI sections; parameters not written in the file use their defaults.

## Basic structure

```ini
[MOSE-Physics]
equations = euler

[MOSE-Numerics]
time-scheme = RK3
cfl = 0.9
time-accurate = false
space-reconstruction = MUSCL
```

The sections group settings by purpose. For the complete list of sections, defaults, allowed values, and whether a parameter can change during a run, use the [Parameter Registry](../user/registry.md).

## Rules that catch common mistakes

| Input situation | What MOSE does |
|---|---|
| A parameter is omitted | Uses its default silently |
| A parameter name has the wrong case or is misspelled | Warns that the key was not read; the parameter keeps its default |
| A key is written without a value | Reports an input error and stops at startup |
| A value has the wrong type or is outside the parameter's allowed range | Reports an input error and stops at startup |
| The file is missing, unreadable, empty, or malformed | Stops at startup |

Parameter names and section names are case-sensitive. Copy their spelling from the registry.

Full-line comments may start with `!`, `;`, or `#`. After a parameter value, use `;` for an inline comment:

```ini
cfl = 0.9 ; convective CFL number
```

A `#` or `!` after a value is treated as part of that value, not as an inline comment.

## What happens when MOSE reads the file again?

MOSE normally re-reads `input.ini` every `ini-diter` iterations (10,000 by default). The registry classifies parameters into three groups:

| Registry label | Meaning |
|---|---|
| **taken** | The new value is applied at the next read |
| **fixed** | The value chosen at startup must remain unchanged |
| **start only** | The setting is used only during startup |

For example, `cfl` and output intervals are read again; physics choices and `sol-format` are fixed; `newrun` is start-only. Changing `newrun` in a file while MOSE is already running does not start a restart.

If a runtime read contains an error, `ini-error = warn` (the default) rejects that reading and keeps the values from the previous successful read. Set `ini-error = stop` to stop the run instead. Errors during the initial read always stop setup.

## Where to look next

- [Parameter Registry](../user/registry.md): full parameter list and exact allowed values.
- [Full input-file guide](../user/input.md): multigrid, multirate, time-step limits, model interactions, and runtime-read details.
- [Restart a simulation](restart.md): how the startup-only `newrun` setting selects an initial condition or a saved solution.
