# Input validation and runtime changes

MOSE reads `input.ini` at startup and normally reads it again every `ini-diter` iterations (10,000 by default). This page explains what the reader accepts and what happens when a later read finds a problem.

For the complete list of parameters and allowed values, use the automatically generated [Parameter Registry](registry.md).

## INI syntax

Each section starts with a name in square brackets. Parameters are written as `key = value`:

```ini
[MOSE-Numerics]
time-scheme = RK3
cfl = 0.9

[MOSE-IO]
sol-diter = 500
```

Section names, parameter names, and values that select named options must use the spelling shown in the registry.

| Syntax | Meaning |
|---|---|
| A line starting with `!`, `;`, or `#` | Full-line comment. |
| `key = value ; comment` | `;` starts an inline comment. |
| `key = value # comment` | `#` is part of the value, not an inline comment. |
| Parameter omitted | MOSE uses the documented default. |

Avoid `#` and `!` after a value. Use `;` for inline comments.

## What is checked at startup?

| Input problem | Startup behaviour |
|---|---|
| File missing, unreadable, empty, or malformed | MOSE stops. |
| Required parameter missing | MOSE stops and names the missing parameter. |
| Unknown or misspelled key | MOSE warns, ignores that key, and uses the parameter default. |
| Parameter value has the wrong type | MOSE stops. |
| Value is outside the allowed rule in the registry | MOSE stops. |
| Real value is NaN or infinity | MOSE stops. |

Examples of values that do not match the expected type include `0,5` for a real, `2.0` for an integer, and `yes` for a logical. Use the parameter's documented syntax rather than relying on permissive parsing.

A key from the registry written without a value is an input error. A key that MOSE does not recognize is different: the reader warns about it and continues with the default for the parameter it was intended to set.

## What happens on a later read?

The runtime reader applies the classification in the registry's **During the run** column:

| Classification | Runtime behaviour | Examples |
|---|---|---|
| `taken` | New value is applied. | `cfl`, `vnn`, output intervals. |
| `fixed` | Startup value must remain unchanged. | `equations`, `turbulence`, `time-scheme`, `time-accurate`, `space-reconstruction`, `sol-format`. |
| `start only` | Later reads do not change the setting. | `newrun`, `ic-format`, `decomposition-map`. |

Keys left out of a later file keep their current values. A fixed key written with a different value is an error of that read. Some setup choices also depend on other settings; MOSE rechecks these combinations when it reads the file again.

**Important:** `newrun = false` written while a simulation is running does not restart it. `newrun` is read only at startup. To continue from a saved solution, stop the current process and launch MOSE again with `newrun = false`.

## `ini-error`: reject a bad runtime read or stop

Set `ini-error` in `[MOSE-IO]`:

```ini
[MOSE-IO]
ini-error = warn
```

| Value | Result when a runtime read contains an error |
|---|---|
| `warn` (default) | Reject the entire read, print a warning, and continue with the values from the previous successful read. |
| `stop` | Stop the run with an error. |

The `ini-error` value in force when the read starts decides what happens. Startup errors always stop the run, regardless of `ini-error`.

## MPI behaviour

In an MPI run, rank 0 reads `input.ini` and distributes its contents to the other ranks. The other ranks do not independently choose different settings from their own local copy of the file.

## Where to look next

- [Parameter Registry](registry.md): all keys, defaults, allowed values, and runtime classifications.
- [Input file guide](input.md): sections and model/numerics interactions.
- [Multigrid](multigrid.md): grid levels and FAS cycle settings.
- [Multirate time stepping](multirate.md): local time-stepping parameters and restrictions.
